#!/usr/bin/env python3

import os
import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


def _detect_viphog_col(df):
    """Return (viphog_col, query_col) by checking which ID column contains 'ViPhOG'."""
    for col in ("target_name", "query_name"):
        if col in df.columns and df[col].str.contains("ViPhOG", na=False).any():
            other = "query_name" if col == "target_name" else "target_name"
            log.debug("ViPhOG column detected: %s  (sequence column: %s)", col, other)
            return col, other
    raise ValueError(
        "Neither 'target_name' nor 'query_name' contains 'ViPhOG' entries. "
        "Check that the input was produced by hmmer_format_table.py."
    )


def ratio_evalue(vphmm_df, taxa_dict, evalue):
    """This function takes a dataframe containing the result of the comparison
       between predicted viral proteins and the ViPhOG database, and outputs a
       table storing the profile hit length ratio and total sequence Evalue for
       each profile-protein hit
    """

    log.debug("ratio_evalue: input df shape=%s, taxa_dict size=%d, evalue cutoff=%s",
              vphmm_df.shape, len(taxa_dict), evalue)
    log.debug("ratio_evalue: input df columns: %s", list(vphmm_df.columns))

    viphog_col, query_col = _detect_viphog_col(vphmm_df)

    informative_df = vphmm_df[
        (vphmm_df[viphog_col].isin(taxa_dict.keys())) &
        (pd.to_numeric(vphmm_df["domain_E-value"], errors="coerce") <= evalue)
    ]

    log.debug("ratio_evalue: rows after filtering=%d (before=%d)",
              len(informative_df), len(vphmm_df))

    if len(informative_df) < 1:
        log.warning("ratio_evalue: no informative hits remain after filtering")
        return None

    informative_df = informative_df.reset_index(drop=True)
    vphmm_hits = list(informative_df[viphog_col].value_counts().index)
    log.debug("ratio_evalue: unique ViPhOG targets after filter: %d", len(vphmm_hits))
    final_hit_list = []

    for vphmm in vphmm_hits:
        vphmm_specific_df = informative_df[informative_df[viphog_col] == vphmm] \
            .reset_index(drop=True)
        query_vcounts = vphmm_specific_df[query_col].value_counts()
        more_than_one = list(query_vcounts[query_vcounts > 1].index)

        query_list = []

        for i in range(len(vphmm_specific_df)):

            query_name = vphmm_specific_df[query_col][i]
            coord_to = vphmm_specific_df["hmm_coord_to"][i]
            coord_from = vphmm_specific_df["hmm_coord_from"][i]
            t_len = vphmm_specific_df["tlen"][i]

            if query_name in query_list:
                continue

            query_list.append(query_name)

            if query_name in more_than_one:
                number = query_vcounts[query_name]
                coords_list = sorted([(
                    vphmm_specific_df["hmm_coord_from"][i + j],
                    vphmm_specific_df["hmm_coord_to"][i + j]) for j in range(number)
                ])
                reference_pair = list(coords_list[0])

                final_coords_list = []
                for elem in range(1, len(coords_list)):
                    if coords_list[elem][0] <= reference_pair[1]:
                        reference_pair[1] = max(
                            reference_pair[1], coords_list[elem][1])
                    else:
                        final_coords_list.append(tuple(reference_pair))
                        reference_pair[0] = coords_list[elem][0]
                        reference_pair[1] = coords_list[elem][1]

                final_coords_list.append(tuple(reference_pair))
                total_hmm_align = sum(
                    [final - initial + 1 for initial, final in final_coords_list]
                )
                hmm_ratio = total_hmm_align / t_len
            else:
                hmm_ratio = (coord_to - coord_from + 1) / t_len

            raw_evalue = vphmm_specific_df["full_E-value"][i]
            log.debug("vphmm=%s query=%s raw E-value=%r", vphmm, query_name, raw_evalue)
            try:
                fs_e_value = float(raw_evalue)
                e_value_exponential = abs(int(("%E" % fs_e_value).split("E")[-1]))
            except Exception as exc:
                log.error("Failed parsing E-value for vphmm=%s query=%s value=%r: %s",
                          vphmm, query_name, raw_evalue, exc)
                raise

            final_hit_list.append((
                vphmm,
                vphmm_specific_df[query_col][i],
                hmm_ratio,
                e_value_exponential
            ))

        final_df = pd.DataFrame(final_hit_list, columns=[
                                "ViPhOG", "query", "Ratio", "Abs_Evalue_exp"])

        # TODO: explain
        final_df["Abs_Evalue_exp"] = final_df["Abs_Evalue_exp"].apply(
            lambda x: 277 if x == 0 else x
        )

        final_df["Taxon"] = final_df["ViPhOG"].apply(lambda x: taxa_dict[x])

    return final_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate dataframe that stores the profile alignment ratio "
                    " and total e-value for each ViPhOG-query pair")
    parser.add_argument("-i", "--input", dest="input_file",
                        help="domtbl generated with Generate_vphmm_hmmer_matrix.py",
                        required=True)
    parser.add_argument("-t", "--taxa", dest="taxa_tsv",
                        help="TSV file: additional_data_vpHMMs_v{1,2,3,4}.tsv", required=True)
    parser.add_argument("-o", "--outfile", dest="out_file",
                        help="Output table name (default: cwd)",
                        default=".")
    parser.add_argument("-e", "--evalue", dest="evalue",
                        help="E-value cutoff for each HMM hit",
                        default=0.01)
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.out_file
    evalue = args.evalue

    log.info("Reading input domtbl: %s", input_file)
    try:
        input_df = pd.read_csv(input_file, sep="\t")
        log.info("Input df loaded: shape=%s, columns=%s", input_df.shape, list(input_df.columns))
        log.debug("Input df head:\n%s", input_df.head(3).to_string())
    except Exception as exc:
        log.error("Failed to read input file %s: %s", input_file, exc)
        raise

    taxa_dict = {}

    log.info("Reading taxa TSV: %s", args.taxa_tsv)
    try:
        tsv_df = pd.read_csv(args.taxa_tsv, sep="\t")
        log.info("Taxa TSV loaded: shape=%s, columns=%s", tsv_df.shape, list(tsv_df.columns))
    except Exception as exc:
        log.error("Failed to read taxa file %s: %s", args.taxa_tsv, exc)
        raise

    for i in range(len(tsv_df)):
        taxa_dict["ViPhOG" + str(tsv_df["Number"][i]) + ".faa"] = tsv_df["Associated"][i]
    log.info("taxa_dict built: %d entries", len(taxa_dict))

    log.info("Running ratio_evalue with evalue cutoff=%s", evalue)
    output_df = ratio_evalue(input_df, taxa_dict, float(evalue))

    log.info("Writing output to: %s", output_file)
    with open(output_file, "w") as of_handle:
        if output_df is None or output_df.empty:
            log.warning("No informative hits — writing empty output")
            print("No informative hits against the ViPhOG database "
                  "were obtained for the contigs provided")
        else:
            log.info("Output df shape=%s", output_df.shape)
            output_df.to_csv(of_handle, sep="\t", index=False)
    log.info("Done")