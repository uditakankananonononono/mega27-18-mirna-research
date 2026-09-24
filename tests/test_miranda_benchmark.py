"""Hermetic checks on the committed miRanda fourth-arm benchmark."""
import csv, json, os

R = os.path.join(os.path.dirname(__file__), "..", "results")
M = json.load(open(os.path.join(R, "miranda_benchmark.json")))
RH = json.load(open(os.path.join(R, "rnahybrid_benchmark.json")))
DB = json.load(open(os.path.join(R, "mirdb_benchmark.json")))


def test_coverage_and_shape():
    assert M["n_mirnas"] == 80 and M["n_rows"] == 7914
    assert 0.0 < M["median_coverage"] <= 1.0
    rows = list(csv.DictReader(open(os.path.join(R, "miranda_rows.csv"))))
    assert len(rows) == M["n_rows"]
    assert all(r["miranda_score"] == "" or float(r["miranda_score"]) > -1e8 for r in rows)


def test_cnn_beats_miranda_significantly():
    assert M["median_per_mirna_cnn"] > M["median_per_mirna_miranda"]
    assert M["wilcoxon_p_miranda_minus_cnn"] < 0.001
    assert M["miranda_wins"] <= 20  # CNN wins the large majority


def test_four_way_ordering():
    tw = M["identical_row_three_way"]
    assert tw["n_rows"] == 7914
    # miRDB (conservation+ML) > CNN > miRanda > RNAhybrid
    assert DB["pooled_mirdb_auroc"] > tw["pooled_auroc_cnn"]
    assert tw["pooled_auroc_cnn"] > tw["pooled_auroc_miranda"]
    assert tw["pooled_auroc_miranda"] > tw["pooled_auroc_rnahybrid"]
    assert abs(tw["pooled_auroc_rnahybrid"] - RH["pooled_auroc_rnahybrid"]) < 1e-9


def test_biophysics_arms_above_chance_pooled():
    assert M["pooled_auroc_miranda"] > 0.55
    assert M["pooled_auroc_miranda"] < M["pooled_auroc_cnn"]
