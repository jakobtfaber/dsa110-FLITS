"""Regression tests for `_gain_marginal_multi_band` (per-channel gain-marginal
evidence for one band).

The function integrates per-component gains g ~ N(0, s2 I_N) analytically. The
key correctness statement (asserted by `test_brute_force_woodbury`) is the
Woodbury identity it is built on: for each channel f with kernels K_f (N, T),
data d_f (T,), noise variance sigma_f^2, the analytic per-channel evidence equals
the direct Gaussian evidence of d_f under the marginal covariance

    Sigma_f = sigma_f^2 I_T + s2 * K_f^T K_f          (T, T)
    lnZ_f   = -0.5 d_f^T Sigma_f^-1 d_f - 0.5 ln det(2 pi Sigma_f)

summed over valid channels. We use a DUCK-TYPED fake model that returns
PRESCRIBED component kernels so the brute side is exactly computable and decoupled
from the real forward model. s2 is passed as a fixed float so the s2-profiling
path is never exercised.
"""

import numpy as np
import pytest

from scattering.scat_analysis.burstfit import FRBParams
from scattering.scat_analysis.burstfit_joint import _gain_marginal_multi_band


class _FakeModel:
    """Minimal duck-typed stand-in for FRBModel.

    `_gain_marginal_multi_band` only touches `.data`, `.noise_std`, `.valid`, and
    `model(p, mk, freq_subset=valid)`. We ignore the params/model_key entirely and
    return a pre-baked per-component kernel cube indexed by call order, so the
    kernels in the evidence math are exactly what the test prescribes.

    kernels: (N, F, T) -- component, full-frequency, time. data: (F, T).
    noise_std: (F,). valid: (F,) bool.
    """

    def __init__(self, kernels, data, noise_std, valid):
        self._kernels = np.asarray(kernels, dtype=float)  # (N, F, T)
        self.data = np.asarray(data, dtype=float)  # (F, T)
        self.noise_std = np.asarray(noise_std, dtype=float)  # (F,)
        self.valid = np.asarray(valid, dtype=bool)  # (F,)
        self._i = 0

    def __call__(self, p, model_key, freq_subset=None):
        # The function calls us once per component, in params_list order. Hand back
        # that component's kernel rows for the requested frequency subset.
        k = self._kernels[self._i]
        self._i += 1
        if freq_subset is not None:
            k = k[freq_subset]
        return k


def _brute_lnZ(kernels, data, noise_std, valid, s2):
    """Direct (non-Woodbury) Gaussian evidence summed over valid channels.

    Per channel: Sigma = sigma^2 I_T + s2 * K^T K  (K is (N, T)); the marginal of
    d ~ N(0, Sigma). Uses slogdet + solve so it is an INDEPENDENT implementation of
    the same quantity the function computes via the Woodbury form.
    """
    K_all = np.asarray(kernels, dtype=float)  # (N, F, T)
    N, F, T = K_all.shape
    total = 0.0
    for f in range(F):
        if not valid[f]:
            continue
        K = K_all[:, f, :]  # (N, T)
        d = np.asarray(data)[f]  # (T,)
        sig2 = float(noise_std[f]) ** 2
        Sigma = sig2 * np.eye(T) + s2 * (K.T @ K)  # (T, T)
        sign, logdet = np.linalg.slogdet(Sigma)
        assert sign > 0, "Sigma must be SPD"
        quad = d @ np.linalg.solve(Sigma, d)
        lnZ_f = -0.5 * quad - 0.5 * (logdet + T * np.log(2.0 * np.pi))
        total += lnZ_f
    return total


def _make_well_conditioned_case(seed, N, F, T, scale=1.0):
    """Seeded random kernels/data/noise with full-rank M_f on every channel.

    Random (N, T) Gaussian kernels are generically full row-rank for T > N, so the
    eigenvalue guard does not fire and the well-conditioned path is exercised.
    """
    rng = np.random.default_rng(seed)
    kernels = rng.standard_normal((N, F, T)) * scale
    data = rng.standard_normal((F, T)) * scale
    noise_std = 0.5 + rng.random(F)  # in [0.5, 1.5), strictly > 1e-9
    valid = np.ones(F, dtype=bool)
    return kernels, data, noise_std, valid


# Dummy params/keys: the fake model ignores their CONTENT, but the function calls
# `replace(p, c0=1.0, gamma=0.0)` on each, so they must be real FRBParams dataclass
# instances (and there must be N of them to drive N component calls).
def _dummy_params_keys(N):
    p = FRBParams(c0=1.0, t0=0.0, gamma=0.0, zeta=0.1, tau_1ghz=0.1, alpha=4.0)
    return [p] * N, ["M3"] * N


def test_brute_force_woodbury():
    """(a) Total lnZ equals the direct brute-force Gaussian evidence.

    This is the load-bearing correctness check: the function's Woodbury-collapsed
    per-channel evidence must equal -0.5 d^T Sigma^-1 d - 0.5 ln det(2 pi Sigma)
    summed over valid channels, with Sigma = sigma^2 I + s2 K^T K. If this fails,
    investigate the kernel/construction or a real discrepancy -- do NOT loosen tol.
    """
    N, F, T = 2, 3, 8
    s2 = 0.7
    kernels, data, noise_std, valid = _make_well_conditioned_case(seed=1234, N=N, F=F, T=T)
    model = _FakeModel(kernels, data, noise_std, valid)
    params, keys = _dummy_params_keys(N)

    lnZ, diag = _gain_marginal_multi_band(model, params, keys, s2=s2)
    brute = _brute_lnZ(kernels, data, noise_std, valid, s2)

    assert np.isfinite(lnZ)
    assert diag["frac_culled"] == 0.0
    assert diag["n_supported"] == F
    assert diag["s2"] == pytest.approx(s2)
    np.testing.assert_allclose(lnZ, brute, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("seed", [0, 7, 42])
def test_brute_force_woodbury_varied(seed):
    """(a, extended) Woodbury cross-check holds across several seeds and a larger
    component count -- guards against an accidental construction that only matches
    for one particular draw."""
    N, F, T = 3, 4, 12
    s2 = 1.3
    kernels, data, noise_std, valid = _make_well_conditioned_case(seed=seed, N=N, F=F, T=T)
    model = _FakeModel(kernels, data, noise_std, valid)
    params, keys = _dummy_params_keys(N)

    lnZ, _ = _gain_marginal_multi_band(model, params, keys, s2=s2)
    brute = _brute_lnZ(kernels, data, noise_std, valid, s2)
    np.testing.assert_allclose(lnZ, brute, rtol=1e-9, atol=1e-9)


def test_label_swap_invariance():
    """(b) Permuting the component order leaves total lnZ unchanged.

    The marginal evidence depends on the components only through K^T K (and b, M),
    which are symmetric under relabelling. Swapping the component axis of the kernel
    cube (and the params/keys we feed in the same order) must give the same lnZ.
    """
    N, F, T = 3, 4, 10
    s2 = 0.9
    kernels, data, noise_std, valid = _make_well_conditioned_case(seed=2024, N=N, F=F, T=T)
    perm = [2, 0, 1]
    kernels_perm = kernels[perm]

    lnZ_a, _ = _gain_marginal_multi_band(
        _FakeModel(kernels, data, noise_std, valid), *_dummy_params_keys(N), s2=s2
    )
    lnZ_b, _ = _gain_marginal_multi_band(
        _FakeModel(kernels_perm, data, noise_std, valid),
        *_dummy_params_keys(N),
        s2=s2,
    )
    np.testing.assert_allclose(lnZ_a, lnZ_b, rtol=1e-12, atol=1e-12)


def test_rank1_fallback_on_collinear_channel():
    """(c) A collinear channel triggers the eigenvalue guard's rank-1 fallback.

    We make the two components IDENTICAL on channel 0 (K_2 == K_1 there), so M_0 is
    rank-deficient (min/max eig ~ 0 < eig_rel_floor) and the channel is culled. The
    code's documented contract for a culled-but-SUPPORTED channel is NOT the gain=0
    baseline but a rank-1 proper-prior evidence on the top eigenpair of M_f:

        emx   = lambda_max(M_f)        (top eigenvalue)
        bproj = b_f . v_top            (data projected on the top eigenvector)
        Ac    = emx + sigma^2 / s2
        lnZ_f = -0.5 (S_dd/sigma^2 - (bproj^2/Ac)/sigma^2)
                - 0.5 T ln(2 pi sigma^2) - 0.5 ln(1 + (s2/sigma^2) emx)

    We replicate that rank-1 value directly and assert the culled channel's
    contribution matches it -- and crucially that it is strictly BELOW the gain=0
    baseline (the merge is Occam-penalized, not rewarded), the invariant the branch
    exists to protect.
    """
    N, F, T = 2, 3, 9
    s2 = 1e6  # large fixed s2: the regime where culling-to-baseline would REWARD a merge
    # Well-conditioned channels 1,2; channel 0 collinear (component 2 == component 1).
    kernels, data, noise_std, valid = _make_well_conditioned_case(seed=99, N=N, F=F, T=T)
    kernels[1, 0, :] = kernels[0, 0, :]  # collapse component 2 onto 1 on channel 0

    model = _FakeModel(kernels, data, noise_std, valid)
    lnZ, diag = _gain_marginal_multi_band(model, *_dummy_params_keys(N), s2=s2)

    assert np.isfinite(lnZ)
    # Exactly one channel (channel 0) is not full-rank-N -> culled.
    assert diag["frac_culled"] == pytest.approx(1.0 / F)
    assert diag["n_supported"] == F - 1  # the two well-conditioned channels

    # --- Recompute the expected per-channel contributions independently. ---
    # Well-conditioned channels via the brute Woodbury (rows 1,2).
    K = kernels
    expected_total = 0.0
    f_cull = 0
    for f in range(F):
        if f == f_cull:
            continue
        Kf = K[:, f, :]  # (N, T)
        d = data[f]
        sig2 = float(noise_std[f]) ** 2
        Sigma = sig2 * np.eye(T) + s2 * (Kf.T @ Kf)
        _, logdet = np.linalg.slogdet(Sigma)
        quad = d @ np.linalg.solve(Sigma, d)
        expected_total += -0.5 * quad - 0.5 * (logdet + T * np.log(2.0 * np.pi))

    # Culled channel 0: rank-1 proper evidence on the top eigenpair of M_0.
    Kc = K[:, f_cull, :]  # (N, T)
    d0 = data[f_cull]
    sig2_0 = float(noise_std[f_cull]) ** 2
    M0 = Kc @ Kc.T  # (N, N) = sum_t K_i K_j
    b0 = Kc @ d0  # (N,)   = sum_t d K_i
    S_dd0 = float(d0 @ d0)
    evals0, evecs0 = np.linalg.eigh(M0)
    emx = evals0[-1]
    vtop = evecs0[:, -1]
    bproj = float(b0 @ vtop)
    Ac = emx + sig2_0 / s2
    quadc = (bproj * bproj) / Ac
    occ_c = np.log1p((s2 / sig2_0) * max(emx, 0.0))
    lnZ_cull_expected = (
        -0.5 * (S_dd0 / sig2_0 - quadc / sig2_0)
        - 0.5 * T * np.log(2.0 * np.pi * sig2_0)
        - 0.5 * occ_c
    )
    expected_total += lnZ_cull_expected

    np.testing.assert_allclose(lnZ, expected_total, rtol=1e-7, atol=1e-7)

    # Documented invariant: the rank-1 fallback sits strictly BELOW the gain=0
    # baseline, so a degenerate merge is penalized (not the +0.5 ln(s2/var) reward
    # that culling-to-baseline would hand back at large s2).
    baseline_cull = -0.5 * S_dd0 / sig2_0 - 0.5 * T * np.log(2.0 * np.pi * sig2_0)
    assert lnZ_cull_expected < baseline_cull


def test_no_valid_channels_returns_neg_inf():
    """All-invalid band -> -inf with frac_culled=1.0 (boundary contract)."""
    N, F, T = 2, 3, 6
    kernels, data, noise_std, _ = _make_well_conditioned_case(seed=5, N=N, F=F, T=T)
    valid = np.zeros(F, dtype=bool)
    model = _FakeModel(kernels, data, noise_std, valid)
    lnZ, diag = _gain_marginal_multi_band(model, *_dummy_params_keys(N), s2=1.0)
    assert lnZ == -np.inf
    assert diag["frac_culled"] == 1.0
    assert diag["n_supported"] == 0


def test_valid_mask_subsets_channels():
    """A partial `valid` mask must restrict the evidence to the unmasked channels.

    Confirms the function honours `valid` (passes it as freq_subset to the model
    and indexes data/noise by it) -- lnZ equals the brute force over valid only.
    """
    N, F, T = 2, 4, 8
    s2 = 1.1
    kernels, data, noise_std, _ = _make_well_conditioned_case(seed=321, N=N, F=F, T=T)
    valid = np.array([True, False, True, True])
    model = _FakeModel(kernels, data, noise_std, valid)
    lnZ, diag = _gain_marginal_multi_band(model, *_dummy_params_keys(N), s2=s2)
    brute = _brute_lnZ(kernels, data, noise_std, valid, s2)
    np.testing.assert_allclose(lnZ, brute, rtol=1e-9, atol=1e-9)
    assert diag["n_supported"] == int(np.count_nonzero(valid))
