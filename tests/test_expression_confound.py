import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from expression_confound import attach_expression, per_mirna, standardize


def test_attach_expression_joins_and_counts_unmatched():
    df = pd.DataFrame({"gene": ["A", "B", "C", "A"], "label": [1, 0, 0, 1]})
    expr = pd.DataFrame({"Gene name": ["A", "B", "B"], "nTPM": [10.0, 1.0, 3.0]})
    out, unmatched = attach_expression(df, expr)
    assert unmatched == 1 and len(out) == 3
    assert np.isclose(out.loc[out.gene == "B", "lexpr"].iloc[0], np.log1p(3.0))


def test_standardize_constant_column_safe():
    X = np.array([[1.0, 5.0], [3.0, 5.0]])
    Z = standardize(X)
    assert np.allclose(Z[:, 1], 0) and np.allclose(Z[:, 0], [-1, 1])


def test_per_mirna_expression_signal_detected_and_cnn_noise_adds_nothing():
    rng = np.random.default_rng(0)
    n = 600
    lexpr = rng.normal(size=n)
    y = (lexpr + 0.3 * rng.normal(size=n) > 0).astype(int)
    g = pd.DataFrame({"n8": rng.integers(0, 3, n), "n7": rng.integers(0, 3, n),
                      "len": rng.integers(500, 3000, n), "lexpr": lexpr,
                      "min": rng.normal(size=n), "sum": rng.normal(size=n), "label": y})
    a, c = per_mirna(g)
    assert a > 0.9 and abs(c - a) < 0.02


def test_per_mirna_returns_none_when_unbalanced():
    g = pd.DataFrame({c: np.zeros(60) for c in ["n8", "n7", "len", "lexpr", "min", "sum"]})
    g["label"] = [1] * 10 + [0] * 50
    assert per_mirna(g) is None
