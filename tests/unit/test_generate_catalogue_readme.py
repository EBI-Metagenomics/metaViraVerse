import gzip
import json
import os
import sys

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "bin")
sys.path.insert(0, BIN_DIR)

import generate_catalogue_readme as gcr  # noqa: E402

REPO_DOCS = os.path.join(os.path.dirname(__file__), "..", "..", "docs")

STATS = {
    "total_sequences": 22572,
    "duplicate_sequences": 1000,
    "deduplicated_sequences": 21572,
    "rrna_excluded_sequences": 1376,
    "not_determined_excluded_sequences": 0,
    "qc_excluded_sequences": 1376,
    "unique_sequences": 20196,
    "viral_sequences": 5426,
    "plasmids": 13528,
    "prophages": 1242,
    "number_of_biomes": 2,
    "viral_clusters": 6616,
    "plasmid_clusters": 12440,
    "total_proteins_viruses": 11000,
    "total_proteins_plasmids": 20000,
    "total_proteins_viruses_reps": 5000,
    "total_proteins_plasmids_reps": 2500,
}


def make_results(tmp_path, version="v1.2.0"):
    res = tmp_path / "results"
    (res / "pipeline_info").mkdir(parents=True, exist_ok=True)
    (res / "viruses").mkdir(exist_ok=True)
    (res / "catalogue.json").write_text(json.dumps(STATS))
    (res / "pipeline_info" / "metaviraverse_software_mqc_versions.yml").write_text(
        f"CHECKV:\n  checkv: 1.0.3\nWorkflow:\n  EBI-Metagenomics/metaviraverse: {version}\n  Nextflow: 25.04.6\n"
    )
    with gzip.open(res / "viruses" / "viruses-all-metadata.tsv.gz", "wt") as f:
        f.write("ID\tBiome\tSource_project\nMGYV1\tmarine,soil\tERP1\nMGYV2\tNA\tmissing\n")
    return res


def run(tmp_path, *extra):
    out = tmp_path / "README.md"
    gcr.main(["-r", str(make_results(tmp_path)), "-o", str(out), *extra])
    return out.read_text()


def test_readme_contains_version_stats_and_methods(tmp_path):
    text = run(tmp_path, "--catalogue-version", "v1.0")
    assert "metaViraVerse pipeline v1.2.0: https://github.com/EBI-Metagenomics/metaViraVerse/releases/tag/v1.2.0" in text
    assert "catalogue.json" not in text and "| Statistic |" not in text
    assert "6,616 viral clusters" in text
    assert "Sequences come from genomes of 1 study." in text
    assert "ERP1" not in text
    assert "2 biome(s): marine, soil" in text
    assert "## How the catalogue was built" in text
    # output descriptions are not copied, only linked
    assert "\n## Viruses" not in text and "Columns of" not in text
    assert "https://github.com/EBI-Metagenomics/metaViraVerse/blob/v1.2.0/docs/output.md" in text
    for page in ("methods.md", "databases.md", "citations.md", "faq.md"):
        assert f"https://github.com/EBI-Metagenomics/metaViraVerse/blob/v1.2.0/docs/{page}" in text
        assert os.path.exists(os.path.join(REPO_DOCS, page)), page


def test_tools_follow_params_and_versions(tmp_path):
    res = make_results(tmp_path)
    (res / "pipeline_info" / "params_2026-01-01_00-00-00.json").write_text(json.dumps(
        {"skip_iphop": True, "predict_host_from_custom_spacers": True, "skip_vitap": True,
         "skip_amrfinderplus_for_viruses": False}))
    out = tmp_path / "README.md"
    gcr.main(["-r", str(res), "-o", str(out)])
    text = out.read_text()
    tools = text.split("### Tools and annotations")[1].split("###")[0]
    assert "[CheckV](https://bitbucket.org/berkeleylab/checkv) v1.0.3" in tools
    assert "iPHoP" not in tools and "SpacePHARER" in tools
    assert "VITAP" not in tools and "geNomad" in tools
    assert "Antimicrobial resistance genes on viruses: [AMRFinderPlus]" in tools
    assert "Antimicrobial resistance genes on plasmids: [DeepARG]" in tools
    assert "linclust at 90% identity and 80% coverage" in tools
    assert "iPHoP" not in text.split("## How the catalogue was built")[1]


def test_exclusion_sentence_follows_skipped_qc(tmp_path):
    res = make_results(tmp_path)
    (res / "pipeline_info" / "params_2026-01-01_00-00-00.json").write_text(json.dumps(
        {"skip_checkv": True, "skip_rrna_detection": True}))
    out = tmp_path / "README.md"
    gcr.main(["-r", str(res), "-o", str(out)])
    text = out.read_text()
    assert "were excluded" not in text and "CheckV" not in text and "barrnap" not in text


def test_dev_version_has_no_release_link(tmp_path):
    text = run(tmp_path, "--pipeline-version", "1.0.0dev")
    assert "pipeline v1.0.0dev (" in text
    assert "releases/tag" not in text
    assert "https://github.com/EBI-Metagenomics/metaViraVerse/blob/dev/docs/output.md" in text



def test_release_tag_ignores_commit_suffix():
    assert gcr.release_tag("v1.0.0") == "v1.0.0"
    assert gcr.release_tag("v1.0.0-gabc1234") == "v1.0.0"
    assert gcr.release_tag("v1.0.0dev") is None
    assert gcr.release_tag("v1.0.0dev-gabc1234") is None
    assert gcr.docs_url("v1.2.0-g0123456789abcdef") == (
        "https://github.com/EBI-Metagenomics/metaViraVerse/blob/v1.2.0/docs/output.md")


def test_metadata_stats(tmp_path):
    viruses = tmp_path / "viruses-all-metadata.tsv.gz"
    with gzip.open(viruses, "wt") as f:
        f.write(
            "ID\tSource_of_prediction\tSource\tcheckv_quality\tcompleteness\tcontamination\tSource_project\n"
            "MGYV1\tgeNomad\tmetagenome\tComplete\t100.0\t0.0\tERP1\n"
            "MGYV2\tgeNomad\tgenome\tLow-quality\t10.0\t5.0\tmissing\n"
            "ID\tSource_of_prediction\tSource\tcheckv_quality\tcompleteness\tcontamination\tSource_project\n"
            "MGYV3\tthird_party\tNA\tNot-determined\tNA\tNA\tERP2\n"
        )
    plasmids = tmp_path / "plasmids-all-metadata.tsv"
    plasmids.write_text("ID\tSource\tSource_lineage\nMGYV4\tgenome\tmissing\nMGYV5\tmetagenome\t\n")

    st = gcr.collect_metadata_stats([str(viruses), str(plasmids)])
    v, p = st["viruses"], st["plasmids"]
    assert v["rows"] == 3
    assert v["categories"]["Source_of_prediction"] == {"geNomad": 2, "third_party": 1}
    assert v["categories"]["Source"] == {"metagenome": 1, "genome": 1, "NA": 1}
    assert v["categories"]["checkv_quality"]["Not-determined"] == 1
    assert v["numeric"]["completeness"] == {"n": 2, "min": 10.0, "max": 100.0, "median": 55.0, "mean": 55.0}
    assert v["missing"]["completeness"] == {"NA": 1}
    assert v["missing"]["Source_project"] == {"missing": 1}
    assert "checkv_quality" not in p["categories"] and p["numeric"] == {}
    assert p["missing"] == {"Source_lineage": {"missing": 1, "": 1}}

    text = "\n".join(gcr.render_metadata_stats(st))
    assert "### Viruses" in text and "### Plasmids" in text
    assert "| geNomad | 2 | 66.7% |" in text
    assert "| `completeness` | 2 | 10.00 | 100.00 | 55.00 | 55.00 |" in text
    assert "| `Source_lineage` | 1 | 1 | 2 | 100.0% |" in text


def test_readme_has_metadata_stats_section(tmp_path):
    text = run(tmp_path)
    section = text.split("## Metadata statistics")[1].split("## Files")[0]
    assert "### Viruses" in section and "2 records." in section
    assert "| `Source_project` | 1 | 0 | 1 | 50.0% |" in section
    assert "| `Biome` | 0 | 1 | 1 | 50.0% |" in section


def test_summary_has_duplicate_and_qc_counts(tmp_path):
    text = run(tmp_path)
    assert ("A total of 22,572 viral sequences, prophages and plasmids were collected. 1,000 identical copies "
            "were merged, giving 21,572 unique sequences; 1,376 of them were removed at the quality-control step "
            "(1,376 carrying rRNA genes), leaving 20,196 sequences in the catalogue.") in text
    assert "checksum (1,000 copies)" in text


def test_summary_of_older_catalogue_json():
    old = {"total_sequences": 10, "unique_sequences": 8, "qc_excluded_sequences": 2}
    assert gcr.summary_bullets(old)[0].startswith("A total of 10 ")
