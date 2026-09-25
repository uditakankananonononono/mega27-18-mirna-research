"""Hermetic checks on the committed BioMart paralog test and its parser."""
import json, os, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "biomart_paralogs.json")))


def test_paralog_counter():
    from biomart_paralogs import paralog_counts
    p = os.path.join(tempfile.mkdtemp(), "p.tsv")
    open(p, "w").write("Gene stable ID\tGene name\tHuman paralogue gene stable ID\nE1\tA\tE2\nE1\tA\tE3\nE1\tA\tE3\nE4\tB\t\n")
    c = paralog_counts(p)
    assert c["A"] == 2 and c["B"] == 0


def test_gates_failed_recorded_and_verdicts():
    assert J["gates"]["G1_pass"] is False and J["gates"]["G2_pass"] is False
    assert J["B1_more_paralogs_vs_uniform"]["verdict"] == "FALSIFIED" and J["B3_loeuf_vs_paralog_matched"]["verdict"] == "CONFIRMED"
