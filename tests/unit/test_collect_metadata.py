import csv
import gzip
import os
import sys

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "bin")
sys.path.insert(0, BIN_DIR)

import collect_metadata  # noqa: E402

META_HEADER = ["sequence_id", "original_name", "description", "type", "source_of_prediction", "biomes",
               "sequence_length", "rrna", "sequence_sha256", "checkv_quality", "definition"]

# (catalogue ID, original name in the rename map, description as choose_sequences.py writes it)
CATALOGUE_VIRUS = ("MGYV1", "MGYG000000001_9 viral_sequence|1:3862", "MGYG000000001_9|viral_sequence-1:3862")
ANALYSIS_VIRUS = ("MGYV2", "ERZ2772590_19|viral_sequence-1:147360", "ERZ2772590_19-viral_sequence-1:147360")
ANALYSIS_PROPHAGE = ("MGYV3", "ERZ2772590_1|prophage-10235:62347", "ERZ2772590_1-prophage-10235:62347")
ANALYSIS_PLASMID = ("MGYV4", "ERZ2772590_7|plasmid-1:5000", "ERZ2772590_7-plasmid-1:5000")
RECORDS = [CATALOGUE_VIRUS, ANALYSIS_VIRUS, ANALYSIS_PROPHAGE, ANALYSIS_PLASMID]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return str(path)


def read_gz_tsv(path):
    with gzip.open(path, "rt", newline="") as f:
        return {row["ID"] if "ID" in row else row["Rep_ID"]: row for row in csv.DictReader(f, delimiter="\t")}


def run(tmp_path, monkeypatch):
    combined = write_tsv(tmp_path / "combined.tsv", META_HEADER, [
        [mgyv, orig, desc, "metagenome", "geNomad", "marine", "100", "No", f"sha{mgyv}", "High-quality",
         "plasmid" if "plasmid" in orig else "virus"]
        for mgyv, orig, desc in RECORDS
    ])
    rename_map = write_tsv(tmp_path / "map.tsv", ["original", "temporary", "short", "biome", "type", "definition"],
                           [[orig, mgyv, orig.split(" ")[0], "marine", "metagenome", "virus"]
                            for mgyv, orig, _ in RECORDS])
    viruses_cluster = write_tsv(tmp_path / "vc.tsv", ["object", "cluster"],
                                [["MGYV1", "MGYV1"], ["MGYV2", "MGYV2"], ["MGYV3", "MGYV2"]])
    plasmids_cluster = write_tsv(tmp_path / "pc.tsv", ["object", "cluster"], [["MGYV4", "MGYV4"]])
    genomes = write_tsv(tmp_path / "genomes.tsv",
                        ["Genome", "Genome_accession", "Lineage", "Sample_accession", "Study_accession"],
                        [["MGYG000000001", "GCA_1", "d__Bacteria", "SAMEA1", "ERP1"]])
    analyses = write_tsv(tmp_path / "analyses.tsv", ["accession", "study_id", "sample_id", "biome_id"],
                         [["ERZ2772590", "ERP2", "SAMEA2", "root:Environmental"]])
    genomad = write_tsv(tmp_path / "genomad.tsv", ["seq_name", "taxonomy"],
                        [["MGYV1", "Viruses;A"], ["MGYV2", "Viruses;B"]])

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "collect_metadata.py", "--combined_meta", combined, "--viruses_cluster", viruses_cluster,
        "--plasmids_cluster", plasmids_cluster, "--additional_metadata", genomes,
        "--analyses_metadata", analyses, "--genomad", genomad, "--map", rename_map,
        "--output_viruses", "viruses.tsv", "--output_plasmids", "plasmids.tsv",
        "--output_reps", "reps.tsv", "--compress",
    ])
    collect_metadata.main()
    return (read_gz_tsv(tmp_path / "viruses.tsv.gz"), read_gz_tsv(tmp_path / "plasmids.tsv.gz"),
            read_gz_tsv(tmp_path / "reps.tsv.gz"))


def test_catalogue_and_analysis_names_are_mapped(tmp_path, monkeypatch):
    viruses, plasmids, _ = run(tmp_path, monkeypatch)
    assert set(viruses) == {"MGYV1", "MGYV2", "MGYV3"}
    assert set(plasmids) == {"MGYV4"}

    assert viruses["MGYV1"]["Original_ID"] == "MGYG000000001_9|viral_sequence-1:3862"
    assert viruses["MGYV2"]["Original_ID"] == "ERZ2772590_19|viral_sequence-1:147360"
    assert viruses["MGYV3"]["Original_ID"] == "ERZ2772590_1|prophage-10235:62347"
    assert plasmids["MGYV4"]["Original_ID"] == "ERZ2772590_7|plasmid-1:5000"

    # values joined on the catalogue ID are filled for analysis sequences too
    assert viruses["MGYV2"]["Cluster_rep"] == "MGYV2"
    assert viruses["MGYV3"]["Cluster_rep"] == "MGYV2"
    assert viruses["MGYV2"]["genomad_lineage"] == "Viruses;B"
    assert plasmids["MGYV4"]["Cluster_rep"] == "MGYV4"


def test_analysis_metadata_joined_on_erz_without_contig_number(tmp_path, monkeypatch):
    viruses, plasmids, _ = run(tmp_path, monkeypatch)
    assert viruses["MGYV2"]["Source_accession"] == "ERZ2772590"
    assert viruses["MGYV2"]["Source_project"] == "ERP2"
    assert plasmids["MGYV4"]["Source_sample"] == "SAMEA2"
    # catalogue genomes still come from the genome metadata
    assert viruses["MGYV1"]["Source_accession"] == "GCA_1"


def test_cluster_stats_include_analysis_members(tmp_path, monkeypatch):
    _, _, reps = run(tmp_path, monkeypatch)
    assert reps["MGYV2"]["cluster_size"] == "2"
    assert reps["MGYV2"]["Rep_genomad_lineage"] == "Viruses;B"
    assert reps["MGYV2"]["Rep_source"] == "metagenome"
