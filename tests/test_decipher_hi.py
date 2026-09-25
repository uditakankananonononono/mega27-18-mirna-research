"""Hermetic checks on the committed DECIPHER HI test and its parser."""
import gzip, json, os, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "decipher_hi.json")))


def test_parser_takes_max_per_gene():
    from decipher_hi import load_hi
    p = os.path.join(tempfile.mkdtemp(), "t.bed.gz")
    with gzip.open(p, "wt") as f:
        f.write("track name='x'\nchr1\t1\t2\tAAA|0.1|50%\t0.1\t.\t1\t2\t0\nchr1\t3\t4\tAAA|0.7|10%\t0.7\t.\t3\t4\t0\n")
    assert load_hi(p)["AAA"] == 0.7


def test_gates_and_verdicts():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"]
    assert J["H3_vs_utrlen"]["verdict"] == "FALSIFIED" and J["H2_vs_expr"]["verdict"] == "CONFIRMED"
    P = J["per_mirna"]
    assert sum(v["cnn"] < v["len_mean"] for v in P.values()) == J["H3_vs_utrlen"]["wins"]
