import gzip
import json
import os
import sys

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import check_catalogue_results as ccr  # noqa: E402

VIRUS_HEADER = ["ID", "Original_ID", "Cluster_rep", "vitap_lineage", "Sequence_length"]


def write(path, text, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)


def tsv(header, rows):
    return "\n".join("\t".join(r) for r in [header, *rows]) + "\n"


def make_results(tmp_path, params=None):
    """A results dir with every file expected under default options."""
    res = tmp_path / "results"
    params = {**ccr.DEFAULT_PARAMS, **(params or {})}
    for _group, patterns, condition in ccr.EXPECTED_FILES:
        if condition is not None and not condition(params):
            continue
        pattern = patterns if isinstance(patterns, str) else patterns[0]
        rel = pattern.replace("*", "x")
        if rel.startswith("pipeline_info/params_"):
            continue  # written below with the real params
        write(res / rel, "{}\n" if rel.endswith(".json") else "content\n", gz=rel.endswith(".gz"))

    write(res / "pipeline_info" / "params_2026-01-01_00-00-00.json", json.dumps(params))
    write(res / "catalogue.json", json.dumps({
        "total_sequences": 6, "duplicate_sequences": 1, "deduplicated_sequences": 5,
        "rrna_excluded_sequences": 1, "not_determined_excluded_sequences": 0, "qc_excluded_sequences": 1,
        "unique_sequences": 4,
        "viral_sequences": 2, "prophages": 1, "viral_clusters": 2, "plasmids": 1,
    }))
    write(res / "additional_data" / "rename_mapfile.tsv",
          "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"
          + "".join(f"s{i}\tMGYV{i}\ts{i}\tsoil\tmetagenome\tvirus\n" for i in range(1, 7)))
    write(res / "viruses" / "viruses-all-metadata.tsv.gz", tsv(VIRUS_HEADER, [
        ["MGYV1", "a|viral_sequence", "MGYV1", "Caudoviricetes", "100"],
        ["MGYV2", "b|prophage", "MGYV1", "missing", "200"],
        ["MGYV3", "c|viral_sequence", "MGYV3", "NA", ""],  # Sequence_length must not be missing
    ]), gz=True)
    write(res / "viruses" / "viruses-cluster-stats.tsv.gz",
          tsv(["Rep_ID", "Rep_biome"], [["MGYV1", "missing"], ["MGYV3", "missing"]]), gz=True)
    write(res / "plasmids" / "plasmids-all-metadata.tsv.gz",
          tsv(["ID", "Biome"], [["MGYV4", "marine"]]), gz=True)
    header = ["sequence_id", "biome"]
    # every unique sequence (5), and the excluded ones: 1 duplicate (not in the table above) and 1 rRNA
    write(res / "additional_data" / "catalogue_filtered_metadata.tsv",
          tsv(header, [["MGYV1", "m"], ["MGYV2", "m"], ["MGYV3", "m"], ["MGYV4", "m"], ["MGYV5", "m"]]))
    write(res / "additional_data" / "catalogue_excluded_records_metadata.tsv",
          tsv(["filter_reason"] + header, [["rrna", "MGYV5", "m"], ["duplicate", "MGYV6", "m"]]))
    write_sequence_set(res, "viruses/viral_sequences", ["MGYV1", "MGYV3"])
    write_sequence_set(res, "viruses/prophages", ["MGYV2"])
    write_sequence_set(res, "plasmids/plasmids", ["MGYV4"])
    write_sequence_set(res, f"{ccr.VR}/viruses_reps", ["MGYV1", "MGYV3"], fna_ext="fasta")
    write_sequence_set(res, f"{ccr.PR}/plasmids_reps", ["MGYV4"], fna_ext="fasta")
    for gff in ccr.ANNOTATED_GFFS:
        write(res / gff, gff_text(["MGYV1"]), gz=True)
    return res


def gff_text(accessions, region_id=None):
    lines = ["##gff-version 3"]
    for acc in accessions:
        lines.append(f"{acc}\tVIRify\tviral_sequence\t1\t100\t.\t.\t.\tID={region_id or acc};quality=HQ")
        lines.append(f"{acc}\tProdigal\tCDS\t1\t60\t.\t+\t0\tID=orig_{acc}_1;product=hypothetical")
    return "\n".join(lines) + "\n"


def write_sequence_set(res, stem, accessions, fna_ext="fna", fna_ids=None, extra_proteins=()):
    fna_ids = fna_ids or accessions
    write(res / f"{stem}.{fna_ext}.gz", "".join(f">{a} description\nACGT\n" for a in fna_ids), gz=True)
    write(res / f"{stem}.gff.gz", gff_text(accessions), gz=True)
    proteins = [f"orig_{a}_1" for a in accessions] + list(extra_proteins)
    write(res / f"{stem}.faa.gz", "".join(f">{p} # 1 # 60\nMK\n" for p in proteins), gz=True)


def run(tmp_path, res):
    out = tmp_path / "report.md"
    code = ccr.main(["-r", str(res), "-o", str(out)])
    return code, out.read_text()


def test_complete_results_files_ok(tmp_path):
    res = make_results(tmp_path)
    rows = ccr.check_files(str(res), ccr.DEFAULT_PARAMS, params_known=True)
    assert [r for r in rows if r["status"] in (ccr.FAIL, ccr.WARN)] == []


def test_missing_and_empty_files_fail(tmp_path):
    res = make_results(tmp_path)
    os.remove(res / "plasmids" / "plasmids.faa.gz")
    write(res / "viruses" / "prophages.fna.gz", "", gz=True)
    rows = {r["pattern"]: r for r in ccr.check_files(str(res), ccr.DEFAULT_PARAMS, params_known=True)}
    assert rows["plasmids/plasmids.faa.gz"]["status"] == ccr.FAIL
    assert rows["viruses/prophages.fna.gz"]["status"] == ccr.FAIL
    assert "empty" in rows["viruses/prophages.fna.gz"]["note"]


def test_execution_reports_one_per_run(tmp_path):
    res = make_results(tmp_path)
    pattern = "pipeline_info/execution_report*.html"
    for f in (res / "pipeline_info").glob("execution_report*.html"):
        f.unlink()

    def status():
        rows = {r["pattern"]: r for r in ccr.check_files(str(res), ccr.DEFAULT_PARAMS, params_known=True)}
        return rows[pattern]["status"], rows[pattern]["note"]

    assert status() == (ccr.FAIL, "missing")
    write(res / "pipeline_info" / "execution_report_2026-01-01_00-00-00.html", "")
    assert status()[0] == ccr.FAIL  # only an empty report
    write(res / "pipeline_info" / "execution_report_2026-01-02_00-00-00.html", "<html></html>")
    write(res / "pipeline_info" / "execution_report_2026-01-03_00-00-00.html", "<html></html>")
    st, note = status()
    assert st == ccr.OK and note.startswith("3 file(s)") and "1 empty" in note


def test_skipped_steps_are_not_expected(tmp_path):
    res = make_results(tmp_path, {"skip_iphop": True})
    rows = {r["pattern"]: r for r in ccr.check_files(str(res), {**ccr.DEFAULT_PARAMS, "skip_iphop": True},
                                                     params_known=True)}
    assert rows[f"{ccr.VR}/host_detection/iphop/iphop_genome.csv"]["status"] == "SKIP"


def test_missing_values_report(tmp_path):
    res = make_results(tmp_path)
    code, text = run(tmp_path, res)
    assert code == 1  # Sequence_length is documented as never missing
    assert "| ❌ | `Sequence_length` | 1 | 33.3% | 0 | 0 | 0 | 1 | documented as never missing |" in text
    assert "| ✅ | `vitap_lineage` | 2 | 66.7% | 1 | 1 | 0 | 0 |  |" in text
    assert "| ⚠️ | `Rep_biome` | 2 | 100.0% | 2 | 0 | 0 | 0 | all values missing |" in text
    # row counts against catalogue.json
    assert "| ✅ | `viruses/viruses-all-metadata.tsv.gz` | `viral_sequences` + `prophages` | 3 | 3 |" in text
    assert "| ✅ | `additional_data/catalogue_excluded_records_metadata.tsv` |" in text


def test_sequence_counts_add_up(tmp_path):
    res = make_results(tmp_path)
    code, text = run(tmp_path, res)
    section = text.split("## 3.")[1].split("## 4.")[0]
    assert "| ✅ | `additional_data/rename_mapfile.tsv` | `total_sequences` | 6 | 6 |" in section
    assert ("| ✅ | `additional_data/catalogue_excluded_records_metadata.tsv` | "
            "`duplicate_sequences` + `qc_excluded_sequences` | 2 | 2 |") in section
    assert ("| ✅ | `additional_data/catalogue_excluded_records_metadata.tsv` filter_reason `rrna` | "
            "`rrna_excluded_sequences` | 1 | 1 |") in section
    assert "| ✅ | `additional_data/catalogue_filtered_metadata.tsv` | `deduplicated_sequences` | 5 | 5 |" in section
    assert "| ✅ | `total_sequences` | 6 | `duplicate_sequences` + `deduplicated_sequences` | 6 |" in section
    assert "| ✅ | `unique_sequences` | 4 | `viral_sequences` + `prophages` + `plasmids` | 4 |" in section

    # one input sequence too many in the rename map
    with open(res / "additional_data" / "rename_mapfile.tsv", "a") as f:
        f.write("s7\tMGYV7\ts7\tsoil\tmetagenome\tvirus\n")
    stats = json.loads((res / "catalogue.json").read_text())
    stats["qc_excluded_sequences"] = 5
    write(res / "catalogue.json", json.dumps(stats))
    code, text = run(tmp_path, res)
    section = text.split("## 3.")[1].split("## 4.")[0]
    assert code == 1
    assert "| ❌ | `additional_data/rename_mapfile.tsv` | `total_sequences` | 6 | 7 |" in section
    assert "| ❌ | `qc_excluded_sequences` | 5 | `rrna_excluded_sequences` + `not_determined_excluded_sequences` | 1 |" in section


def test_row_count_mismatch_fails(tmp_path):
    res = make_results(tmp_path)
    write(res / "catalogue.json", json.dumps({"plasmids": 5}))
    _code, text = run(tmp_path, res)
    assert "| ❌ | `plasmids/plasmids-all-metadata.tsv.gz` | `plasmids` | 5 | 1 |" in text


def test_documented_columns_parsed():
    docs = ccr.parse_documented_columns(ccr.DEFAULT_OUTPUT_DOCS)
    assert docs["viruses-all-metadata.tsv.gz"]["ID"] is False
    assert docs["viruses-all-metadata.tsv.gz"]["vitap_lineage"] is True
    assert "Rep_ID" in docs["viruses-cluster-stats.tsv.gz"]


def id_checks(res):
    stats = json.loads((res / "catalogue.json").read_text())
    return {r["name"] if not r.get("annotated") else r["gff"]: r
            for r in ccr.check_sequence_ids(str(res), "MGYV", stats)}


def test_sequence_ids_ok(tmp_path):
    checks = id_checks(make_results(tmp_path))
    assert all(r["status"] == ccr.OK for r in checks.values()), [r["problems"] for r in checks.values()]
    vs = checks["Viral sequences"]
    assert (vs["n_seqs"], vs["exp_seqs"], vs["n_cds"], vs["n_prot"]) == (2, 2, 2, 2)


def test_unrenamed_fasta_header_fails(tmp_path):
    res = make_results(tmp_path)
    write_sequence_set(res, "plasmids/plasmids", ["MGYV4"], fna_ids=["MGYG000519004_12|plasmid-1:500"])
    r = id_checks(res)["Plasmids"]
    assert r["status"] == ccr.FAIL
    messages = " ".join(m for _, m, _ in r["problems"])
    assert "1 FASTA header(s) are not `MGYV` accessions" in messages
    assert "GFF sequence(s) not in the FASTA" in messages
    assert ["MGYG000519004_12|plasmid-1:500"] in [ex for _, _, ex in r["problems"]]


def test_unrenamed_gff_id_fails_but_cds_ids_are_ignored(tmp_path):
    res = make_results(tmp_path)
    write(res / "viruses" / "prophages.gff.gz", gff_text(["MGYV2"], region_id="MGYG1_5|prophage-1:100"), gz=True)
    r = id_checks(res)["Prophages"]
    messages = [m for _, m, _ in r["problems"]]
    assert messages == ["1 non-CDS GFF record ID(s) are not `MGYV` accessions"]


def test_cds_and_protein_counts_must_match(tmp_path):
    res = make_results(tmp_path)
    write_sequence_set(res, "viruses/viral_sequences", ["MGYV1", "MGYV3"], extra_proteins=["stray_1"])
    r = id_checks(res)["Viral sequences"]
    messages = " ".join(m for _, m, _ in r["problems"])
    assert "2 CDS records in the GFF but 3 proteins in the FAA" in messages
    assert "1 FAA protein(s) without a CDS in the GFF" in messages


def test_report_has_identifier_section(tmp_path):
    _code, text = run(tmp_path, make_results(tmp_path))
    assert "## 4. Sequence identifiers and protein counts" in text
    assert "| ✅ | Viral sequences | 2 | 2 | 2 | 2 | 2 | - |" in text


def read_mismatches(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    assert rows[0] == ["issue", "id"]
    return sorted(tuple(r) for r in rows[1:])


def test_mismatches_written_only_when_asked_and_found(tmp_path):
    res = make_results(tmp_path)
    # FASTA: MGYV1 and an extra MGYV9 (twice); GFF: MGYV1, MGYV3; FAA: CDS proteins plus a stray one
    write(res / "viruses" / "viral_sequences.fna.gz", ">MGYV1\nA\n>MGYV9\nA\n>MGYV9\nA\n", gz=True)
    write(res / "viruses" / "viral_sequences.faa.gz",
          ">orig_MGYV1_1\nM\n>stray_1\nM\n>stray_1\nM\n", gz=True)
    stats = json.loads((res / "catalogue.json").read_text())

    out_dir = tmp_path / "mismatches"
    assert all(r.get("mismatch_file") is None for r in ccr.check_sequence_ids(str(res), "MGYV", stats))
    assert not out_dir.exists()

    checks = {r.get("name"): r for r in ccr.check_sequence_ids(str(res), "MGYV", stats, str(out_dir))}
    path = checks["Viral sequences"]["mismatch_file"]
    assert path == str(out_dir / "viral_sequences.mismatches.tsv")
    assert read_mismatches(path) == sorted([
        ("in_gff_not_in_fasta", "MGYV3"),
        ("in_fasta_not_in_gff", "MGYV9"),
        ("duplicated_in_fasta", "MGYV9"),
        ("in_faa_not_in_gff", "stray_1"),
        ("in_gff_not_in_faa", "orig_MGYV3_1"),
        ("duplicated_in_faa", "stray_1"),
    ])
    # consistent sets get no file
    assert checks["Prophages"]["mismatch_file"] is None
    assert sorted(os.listdir(out_dir)) == ["viral_sequences.mismatches.tsv"]

    messages = " ".join(m for _, m, _ in checks["Viral sequences"]["problems"])
    assert "1 FASTA sequence(s) not in the GFF" in messages


def test_report_lists_mismatch_files(tmp_path):
    res = make_results(tmp_path)
    write_sequence_set(res, "viruses/viral_sequences", ["MGYV1", "MGYV3"], extra_proteins=["stray_1"])
    out = tmp_path / "report.md"
    ccr.main(["-r", str(res), "-o", str(out), "--mismatches-dir", str(tmp_path / "mm")])
    text = out.read_text()
    assert "### Mismatched identifiers" in text
    assert "viral_sequences.mismatches.tsv" in text
