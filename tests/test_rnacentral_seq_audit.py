"""Hermetic checks on the committed RNAcentral sequence audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "rnacentral_seq_audit.json")))


def test_gates():
    g = J["gates"]
    assert g["G1_let7a_exact"] and g["G2_pass"] and g["G3_pass"] and g["G3_deranged_seed_identity"] < 0.2


def test_identity_confirmed_and_consistent():
    R = [r for r in J["rows"] if r["resolved"]]
    assert len(R) == J["n_resolved"] == 132
    assert all(r["seed_match"] and r["full_match"] for r in R)
    assert J["H1_seed_identity"]["verdict"] == "CONFIRMED" and J["H2_full_identity"]["verdict"] == "CONFIRMED"


def test_seeds_are_rna_and_length_ok():
    for r in J["rows"]:
        s = r["targetscan_seq"]
        assert set(s) <= set("ACGU") and 17 <= len(s) <= 28 and r["mimat"].startswith("MIMAT")


def test_nothing_unresolved():
    assert J["unresolved"] == [] and J["n_resolved"] == J["n_with_targetscan_accession"] == 132
