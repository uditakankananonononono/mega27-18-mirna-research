"""Hermetic checks on the committed ENCODE miRNA-abundance audit."""
import json, os, sys
import pytest

HERE = os.path.dirname(__file__)
J = json.load(open(os.path.join(HERE, "..", "results", "encode_mirna_abundance.json")))


def test_files_and_lines():
    assert J["n_files"] == 21 and J["n_lines"] == 12
    assert len({f["acc"] for f in J["files"]}) == 21 and all(f["acc"].startswith("ENCFF") for f in J["files"])


def test_gates_recorded_honestly():
    g = J["gates"]
    assert g["G0_rows_ok"] and g["G1_miR21_top_decile"] and g["G3_pass"]
    assert g["G2_pass"] is False and all(v == 0.0 for v in g["G2_miR122_per_line"].values())
    ph = J["posthoc_G2b_tissue_controls"]
    assert "not pre-registered" in ph["note"] and ph["all_pass"]
    assert set(ph["MIR142_top3"]) == {"GM12878", "K562", "HL-60"}


def test_hypotheses_falsified():
    for h in ["H1_abundance_vs_cnn_auroc", "H2_abundance_vs_delta", "H3_partial_given_log_npos"]:
        assert J[h]["verdict"] == "FALSIFIED" and J[h]["p_one_sided"] >= 0.05 and J[h]["n"] == 131
        assert abs(J[h]["rho"]) < 0.2


def test_precursor_mapping_rules():
    sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
    pytest.importorskip("requests")
    from encode_mirna_abundance import precursors
    S = ["MIR1-1", "MIR1-2", "MIR10A", "MIR100", "MIRLET7A1", "MIRLET7A2", "MIRLET7A3", "MIR125B1", "MIR125B2", "MIR21", "MIR210"]
    assert precursors("hsa-miR-1-3p", S) == ["MIR1-1", "MIR1-2"]
    assert precursors("hsa-let-7a-5p", S) == ["MIRLET7A1", "MIRLET7A2", "MIRLET7A3"]
    assert precursors("hsa-miR-125b-5p", S) == ["MIR125B1", "MIR125B2"]
    assert precursors("hsa-miR-21-5p", S) == ["MIR21"]


def test_partial_spearman_sanity():
    sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
    pytest.importorskip("requests")
    import numpy as np
    from encode_mirna_abundance import partial_spearman
    rng = np.random.default_rng(0); z = rng.normal(size=200)
    x = z + 0.1 * rng.normal(size=200); y = z + 0.1 * rng.normal(size=200)
    r, p = partial_spearman(x, y, z)
    assert abs(r) < 0.3 and p > 0.001
