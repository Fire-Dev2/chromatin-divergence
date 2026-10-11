"""Smoke tests for the Aim 3 generator, runner and scoring."""
import numpy as np
import pandas as pd

from cdtest import crt, simulate as S


def test_generator_targets_partial_r2_and_structure():
    s = S.Setting(0.3, 0.6, 0.3, "hetero", "large", G=2000)
    r2 = []
    for i in range(20):
        sim = S.simulate_replicate(s, np.random.default_rng(i), n_decoys=3)
        Z = crt.build_design(sim["C"], S.SPLINE_COLS)
        r2.append(crt.partial_r2(sim["M"], sim["x"], Z))
        assert sim["decoys"].shape == (2000, 3)
        assert abs(np.corrcoef(sim["x"], sim["C"]["d"])[0, 1] - 0.3) < 0.1
    assert abs(np.mean(r2) - 0.05) < 0.01
    assert len(S.grid_settings()) == 192


def test_runner_checkpoint_resume_and_summary(tmp_path):
    st = [S.Setting(0.5, 0.6, 0.0, "homo", e, G=300) for e in ("zero", "large")]
    kw = dict(reps=4, B=19, n_jobs=1, n_decoys=9, batch=2, settings=st, verbose=0,
              part_dir=str(tmp_path))
    a = S.run_grid(**kw)
    assert len(a) == 8 and len(list(tmp_path.glob("*.parquet"))) == 4
    b = S.run_grid(**kw)  # resumes from parts
    pd.testing.assert_frame_equal(a, b, check_dtype=False)
    sm = S.summarize(a)
    assert set(sm.test) == set(S.ALL_TESTS)
    assert {"rej05", "rej01", "coverage", "bias", "power_matched"} <= set(sm.columns)


def test_module_designs():
    rng = np.random.default_rng(0)
    raw = np.array([S.module_size(rng, "mixed") for _ in range(20000)])
    assert raw.min() >= 1 and raw.max() <= 150
    assert abs(np.median(raw) - 8) <= 1 and raw.mean() > 11  # long right tail
    g5 = np.bincount(S.make_genome(2000, np.random.default_rng(1), 5.0).module).mean()
    g25 = np.bincount(S.make_genome(2000, np.random.default_rng(1), 25.0).module).mean()
    assert g5 < g25
    s = S.Setting(0.3, 0.3, 0.0, "homo", "zero", module_design="mixed")
    assert S.simulate_replicate(s, np.random.default_rng(2), 0)["M"].shape == (2000,)
