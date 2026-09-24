import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from rnahybrid_benchmark import parse_compact


def test_parse_compact_good_line():
    assert parse_compact("GENE1:2500:q:22:-18.4:0.031:120:rest") == ("GENE1", -18.4)


def test_parse_compact_rejects_malformed():
    assert parse_compact("too:short") is None
    assert parse_compact("a:b:c:d:notanumber:f") is None


def test_parse_compact_keeps_negative_energy():
    t, e = parse_compact("X:100:q:21:-0.5:1.0:3:z")
    assert t == "X" and e == -0.5
