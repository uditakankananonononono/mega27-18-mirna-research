"""Hermetic checks on the committed Complex Portal subunit test."""
import json, os, re, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "complexportal_subunits.json")))
P = J["per_mirna"]


def test_accession_regex():
    from complexportal_subunits import ACC
    assert ACC.findall("P84022(1)|Q13485(1)|CPX-1(2)|A0A024R161(1)") == ["P84022", "Q13485", "A0A024R161"]


def test_gates_and_verdicts():
    g = J["gates"]
    assert g["G1_pass"] and g["G2_pass"] and g["G3_pass"]
    assert J["X1_vs_uniform"]["verdict"] == "CONFIRMED" and J["X2_vs_expr"]["verdict"] == "FALSIFIED" and J["X3_vs_utrlen"]["verdict"] == "FALSIFIED"
    assert sum(v["cnn"] < v["expr_mean"] for v in P.values()) == J["X2_vs_expr"]["wins"]
