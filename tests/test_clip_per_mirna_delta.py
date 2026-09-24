import os, sys, importlib.util
import numpy as np, pandas as pd
HERE = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('cpd', os.path.join(HERE, '..', 'scripts', 'clip_per_mirna_delta.py'))
cpd = importlib.util.module_from_spec(spec); spec.loader.exec_module(cpd)


def _synth(n=400, seed=0, signal=True):
    rng = np.random.default_rng(seed)
    n8 = rng.integers(0, 3, n); n7 = rng.integers(0, 4, n); ln = rng.integers(300, 3000, n)
    mn = rng.normal(0, 1, n)
    p = 1 / (1 + np.exp(-(0.4 * n8 + 0.2 * n7 + (1.2 * mn if signal else 0))))
    y = rng.binomial(1, p)
    return pd.DataFrame({'n8': n8, 'n7': n7, 'len': ln, 'min': mn, 'sum': mn * 2, 'label': y})


def test_delta_positive_when_cnn_feature_carries_signal():
    r = cpd.per_mirna_delta(_synth(signal=True), min_pos=20, min_neg=20)
    assert r is not None and r[2] > 0.02


def test_delta_near_zero_when_cnn_features_are_noise():
    r = cpd.per_mirna_delta(_synth(signal=False), min_pos=20, min_neg=20)
    assert r is not None and abs(r[2]) < 0.05


def test_class_imbalance_returns_none():
    df = _synth(n=100); df['label'] = 0; df.loc[df.index[:5], 'label'] = 1
    assert cpd.per_mirna_delta(df) is None
