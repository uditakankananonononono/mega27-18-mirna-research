"""Hermetic checks on the committed MGI test and its report parser."""
import json, os, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "mgi_mouse_ko.json")))


def test_hmd_parser():
    from mgi_mouse_ko import load_hmd
    p = os.path.join(tempfile.mkdtemp(), "h.rpt")
    open(p, "w").write("AAA\t1\tAaa\tMGI:1\tMP:0010768, MP:0005380\t\nBBB\t2\tBbb\tMGI:2\t\t\nAAA\t1\tAaa2\tMGI:3\tMP:0005376\t\n")
    rows, ph = load_hmd(p)
    assert rows == {"AAA", "BBB"} and ph == {"AAA": {"MP:0010768", "MP:0005380", "MP:0005376"}}


def test_verdicts_mixed_recorded():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"]
    v = [J[k]["verdict"] for k in ("M1_mort_vs_expr", "M2_mort_vs_pub", "M3_mort_vs_utrlen", "M4_emb_vs_pub")]
    assert v == ["CONFIRMED", "FALSIFIED", "FALSIFIED", "CONFIRMED"]
    assert all(J[k]["p_one_sided"] > 0.0125 for k in ("M1_mort_vs_expr", "M4_emb_vs_pub"))
