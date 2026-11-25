import gzip
import os
import sys

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "bin")
sys.path.insert(0, BIN_DIR)

import collect_catalogue_stats as ccs  # noqa: E402


def test_reps_list_without_header_counts_every_line(tmp_path):
    # extract_reps_stats.py --output-reps-list: one representative per line, no header
    reps = tmp_path / "viruses_reps.tsv"
    reps.write_text("MGYG1_1|viral_sequence-1:100\nMGYG2_3|prophage-5:900\nERZ1_2|viral_sequence-1:50\n")
    assert ccs.count_clusters(str(reps)) == 3


def test_header_is_skipped_when_present(tmp_path):
    table = tmp_path / "clusters.tsv.gz"
    with gzip.open(table, "wt") as f:
        f.write("object\tcluster\nA\t0\nB\t0\nC\t1\n")
    assert ccs.count_clusters(str(table)) == 3  # unique values of column 1


def test_single_cluster(tmp_path):
    reps = tmp_path / "plasmids_reps.tsv"
    reps.write_text("MGYG1_7|plasmid-1:2000\n")
    assert ccs.count_clusters(str(reps)) == 1


def test_main_sequence_counts(tmp_path, monkeypatch, capsys):
    import json
    (tmp_path / "stats.json").write_text(json.dumps({
        "number_of_viral_sequences": 2, "number_of_viral_sequence_proteins": 4,
        "number_of_prophages": 1, "number_of_prophage_proteins": 2,
        "number_of_plasmids": 1, "number_of_plasmid_proteins": 3,
    }))
    # 8 input sequences
    (tmp_path / "map.tsv").write_text(
        "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"
        + "".join(f"s{i}\tMGYV{i}\ts{i}\tsoil\tmetagenome\tvirus\n" for i in range(8)))
    # 6 unique after deduplication
    head = "sequence_id\toriginal_name\tbiomes\n"
    (tmp_path / "all.tsv").write_text(head + "".join(f"MGYV{i}\ts{i}\tsoil\n" for i in range(6)))
    # 4 kept after QC
    (tmp_path / "kept.tsv").write_text(head + "".join(f"MGYV{i}\ts{i}\tsoil\n" for i in range(4)))
    (tmp_path / "excl.tsv").write_text(
        "filter_reason\t" + head
        + "rrna\tMGYV4\ts4\tsoil\nnot_determined\tMGYV5\ts5\tsoil\n"
        + "duplicate\tMGYV6\ts6\tsoil\nduplicate\tMGYV7\ts7\tsoil\n")
    (tmp_path / "reps.tsv").write_text("MGYV0\nMGYV2\n")
    (tmp_path / "p.faa").write_text(">a\nM\n")
    out = tmp_path / "catalogue.json"
    monkeypatch.setattr(sys, "argv", [
        "collect_catalogue_stats.py", "--initial-json", str(tmp_path / "stats.json"),
        "--map-file", str(tmp_path / "map.tsv"), "--initial-metadata", str(tmp_path / "all.tsv"),
        "--metadata", str(tmp_path / "kept.tsv"), "--excluded-metadata", str(tmp_path / "excl.tsv"),
        "--clusters-viruses", str(tmp_path / "reps.tsv"), "--clusters-plasmids", str(tmp_path / "reps.tsv"),
        "--viruses-reps-proteins", str(tmp_path / "p.faa"), "--plasmids-reps-proteins", str(tmp_path / "p.faa"),
        "-o", str(out)])
    ccs.main()
    stats = json.loads(out.read_text())
    assert stats["total_sequences"] == 8
    assert stats["duplicate_sequences"] == 2
    assert stats["deduplicated_sequences"] == 6
    assert stats["rrna_excluded_sequences"] == 1
    assert stats["not_determined_excluded_sequences"] == 1
    assert stats["qc_excluded_sequences"] == 2
    assert stats["unique_sequences"] == 4
    assert "Warning" not in capsys.readouterr().err


def test_repeated_header_is_not_counted(tmp_path):
    table = tmp_path / "clusters.tsv"
    table.write_text("object\tcluster\nA\t0\nB\t0\nobject\tcluster\nC\t1\n")
    assert ccs.count_clusters(str(table)) == 3


def test_no_header_counts_first_line(tmp_path):
    reps = tmp_path / "reps.tsv"
    reps.write_text("MGYV1\nMGYV2\nMGYV2\n")
    assert ccs.count_clusters(str(reps)) == 2
