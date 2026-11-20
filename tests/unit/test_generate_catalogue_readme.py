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
    "unique_sequences": 20196,
    "qc_excluded_sequences": 2380,
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
    assert "## Viruses" not in text and "Columns of" not in text
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
