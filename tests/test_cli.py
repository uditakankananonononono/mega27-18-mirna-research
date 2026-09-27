"""Hermetic tests for mirtarget-score (bundled model + data, tiny FASTA)."""
import json, pathlib, subprocess, sys, os

REPO = pathlib.Path(__file__).resolve().parents[1]
ENV = {**os.environ, "PYTHONPATH": str(REPO / "src")}


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "mirna.cli", *args],
                          capture_output=True, text=True, env=ENV, timeout=600)


def test_named_target_mode_reports_missing_original_inputs():
    """Named target mode cannot run without its unbundled TargetScan inputs."""
    import tempfile
    with tempfile.TemporaryDirectory() as empty_data:
        r = run_cli("--mirna", "hsa-miR-7-5p", "--gene", "CDR1as", "--top", "5",
                    "--data-dir", empty_data)
    assert r.returncode != 0
    assert "TargetScan miR_Family_Info.txt is absent" in r.stderr

def test_known_sponge_hit_with_original_inputs():
    """Integration test, only if original TargetScan inputs are mounted."""
    data = REPO / "data"
    if not all((data / name).is_file() for name in ("miR_Family_Info.txt", "human_utrs.tsv")):
        import pytest
        pytest.skip("original TargetScan files not bundled; no scientific fixture substituted")
    r = run_cli("--mirna", "hsa-miR-7-5p", "--gene", "CDR1as", "--top", "5")
    assert r.returncode == 0, r.stderr
    out = json.loads(r.stdout)
    assert out["sites"] and all(s["target"].startswith("CDR1as") for s in out["sites"])
    assert out["sites"][0]["score"] < -1.0


def test_custom_seq_fasta(tmp_path):
    # synthetic UTR containing the let-7a 8mer seed complement
    sys.path.insert(0, str(REPO / "src"))
    from mirna.cli import rc
    seed8 = rc("UGAGGUAGUAGGUUGUAUAGUU"[1:8])
    utr = "ACGU" * 40 + seed8 + "A" + "UUUU" * 20
    f = tmp_path / "u.fa"; f.write_text(f">syn1\n{utr}\n")
    r = run_cli("--seq", "UGAGGUAGUAGGUUGUAUAGUU", "--utr-fasta", str(f))
    assert r.returncode == 0, r.stderr
    out = json.loads(r.stdout)
    assert out["sites"] and out["sites"][0]["target"] == "syn1"
    assert out["sites"][0]["site_type"] == "8mer"


def test_no_match_is_empty_not_error(tmp_path):
    f = tmp_path / "u.fa"; f.write_text(">syn2\n" + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    r = run_cli("--seq", "UGGAAGACUAGUGAUUUUGUUGU", "--utr-fasta", str(f))
    assert r.returncode == 0 and json.loads(r.stdout)["sites"] == []
