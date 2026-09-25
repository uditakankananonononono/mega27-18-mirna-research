"""Hermetic checks on the committed HGNC family test and its parser."""
import json, os, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "hgnc_families.json")))


def test_hgnc_parser():
    from hgnc_families import load_hgnc, g2_ratio
    p = os.path.join(tempfile.mkdtemp(), "h.txt")
    rows = ["symbol\tlocus_group\tstatus\tlocation\tgene_group\tgene_group_id",
            "ZNF1\tprotein-coding gene\tApproved\t19q13.4\tZinc fingers C2H2-type\t28",
            "ZNF2\tprotein-coding gene\tApproved\t7q11\tZinc fingers C2H2-type|KRAB domain containing\t28|1",
            "GENEA\tprotein-coding gene\tApproved\t19p13\t\t",
            "GENEB\tprotein-coding gene\tApproved\t1p36\tX family\t5",
            "PSEUD\tpseudogene\tApproved\t19q13\tZinc fingers C2H2-type\t28"]
    open(p, "w").write("\n".join(rows) + "\n")
    h = load_hgnc(p)
    assert "PSEUD" not in h.index and bool(h.loc["ZNF2", "znf"]) and not h.loc["GENEB", "znf"]
    assert abs(h.loc["ZNF1", "fam"] - 1.584962500721156) < 1e-9 and h.loc["GENEA", "fam"] == 0
    assert abs(g2_ratio(h) - 1.0) < 1e-9


def test_verdicts_and_failed_gate_recorded():
    assert J["gates"]["G1_pass"] is False and J["gates"]["G2_pass"] is True
    assert J["F1_znf_vs_expr"]["verdict"] == "FALSIFIED" and J["F3_loeuf_noznf_vs_expr"]["verdict"] == "CONFIRMED"
