"""Tests for emptycell/cell.py. Each one checks something the experiment depends on."""
import itertools


import numpy as np
import pytest

from emptycell import cell


def haar_random_pure(n, rng):
    v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    v /= np.linalg.norm(v)
    return np.outer(v, v.conj())


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))



def _haar_unitary(d, rng):
    z = (rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))) / np.sqrt(2)
    Q, R = np.linalg.qr(z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def _rho(g, d, rng):
    return (1 - g) * haar_random_pure(int(np.log2(d)), rng) + g * np.eye(d) / d


# The sampling shortcut gives the same numbers as real random unitaries

def test_dirichlet_shortcut_matches_explicit_haar_basis():
    """Chance that two shots match, once with a real random state and unitary, once with the shortcut.
    Both should have the same distribution."""
    from scipy.stats import ks_2samp
    rng = np.random.default_rng(1)
    d, g, T = 8, 0.3, 4000
    a = np.empty(T)
    b = np.empty(T)
    for t in range(T):
        rho = _rho(g, d, rng)
        U = _haar_unitary(d, rng)
        P = np.real(np.diag(U @ rho @ U.conj().T))
        a[t] = np.sum(P ** 2)
        q = cell.dirichlet_flat(d, rng)
        b[t] = np.sum(((1 - g) * q + g / d) ** 2)
    assert ks_2samp(a, b).pvalue > 1e-3
    # its average should be (1 + Tr ρ²)/(d + 1)
    assert abs(b.mean() - (1 + cell.purity_of(g, d)) / (d + 1)) < 4 * b.std() / np.sqrt(T)


def test_single_shot_per_fresh_instance_carries_no_information():
    """One copy of a sample gives a uniformly random outcome no matter what g is, because the average of ρ over
    ψ is I/d. That's why the classical agents need more than one copy per sample."""
    from scipy.stats import chisquare
    rng = np.random.default_rng(2)
    d, T = 16, 20000
    for g in (0.0, 0.5):
        s = np.array([cell._shots(cell.dirichlet_flat(d, rng), g, 1, rng)[0] for _ in range(T)])
        assert chisquare(np.bincount(s, minlength=d)).pvalue > 1e-3


# Every purity guess is right on average

@pytest.mark.parametrize("g", [0.0, 0.2, 0.5])
def test_collision_label_unbiased(g):
    rng = np.random.default_rng(3)
    d, k, T = 16, 12, 6000
    y = np.array([cell.collision_label(cell.dirichlet_flat(d, rng), g, k, rng) for _ in range(T)])
    assert abs(y.mean() - cell.purity_of(g, d)) < 4 * y.std() / np.sqrt(T)


@pytest.mark.parametrize("g", [0.0, 0.5])
def test_localrm_label_unbiased(g):
    rng = np.random.default_rng(4)
    n, k, T = 3, 10, 6000
    W = cell._hamming_weights(n)
    y = np.array([cell.localrm_label(n, g, k, W, rng) for _ in range(T)])
    assert abs(y.mean() - cell.purity_of(g, 2 ** n)) < 4 * y.std() / np.sqrt(T)


def test_compbasis_label_unbiased():
    rng = np.random.default_rng(5)
    d, g, k, T = 32, 0.35, 5, 40000
    y = np.array([cell.compbasis_label(d, g, k, rng) for _ in range(T)])
    assert abs(y.mean() - cell.purity_of(g, d)) < 4 * y.std() / np.sqrt(T)


def test_swap_label_unbiased():
    rng = np.random.default_rng(6)
    f = np.full(200000, cell.purity_of(0.4, 64))
    y = cell.swap_labels(f, rng)
    assert abs(y.mean() - f[0]) < 4 / np.sqrt(f.size)


def test_purity_formula_matches_density_matrix():
    rng = np.random.default_rng(7)
    for g in (0.0, 0.17, 0.5):
        assert abs(purity(_rho(g, 8, rng)) - cell.purity_of(g, 8)) < 1e-12


# The SWAP test is a real measurement

def _bell_even_singlet_projector(n):
    """Measure each qubit pair (A_i, B_i) in the Bell basis and accept if the number of singlets is even."""
    singlet = np.array([0, 1, -1, 0]) / np.sqrt(2)
    Ps = np.outer(singlet, singlet)
    Pt = np.eye(4) - Ps
    D = 4 ** n
    proj = np.zeros((D, D))
    for pattern in itertools.product([0, 1], repeat=n):
        if sum(pattern) % 2:
            continue
        op = np.array([[1.0]])
        for p in pattern:
            op = np.kron(op, Ps if p else Pt)
        proj += op
    # reorder from A0 B0 A1 B1 ... to A0 A1 ... B0 B1 ...
    perm = [2 * i for i in range(n)] + [2 * i + 1 for i in range(n)]
    t = proj.reshape([2] * (4 * n))
    t = np.transpose(t, perm + [2 * n + p for p in perm])
    return t.reshape(D, D)


def test_swap_test_acceptance_from_bell_measurement():
    rng = np.random.default_rng(8)
    n = 2
    d = 2 ** n
    P = _bell_even_singlet_projector(n)
    for g in (0.0, 0.3, 0.5):
        rho = _rho(g, d, rng)
        acc = np.real(np.trace(np.kron(rho, rho) @ P))
        assert abs(acc - (1 + purity(rho)) / 2) < 1e-12


# The inequality the lower bound rests on (Chen, Cotler, Huang, Li, eq. 120 and 121)

def _perm(G):
    k = G.shape[0]
    return sum(np.prod([G[i, p[i]] for i in range(k)]) for p in itertools.permutations(range(k)))


def _leaf_ratio_exact(phis, d):
    """E_ψ Π_t d|⟨φ_t|ψ⟩|² = d^k perm(G) / (d(d+1)...(d+k-1)), with G the Gram matrix of the φ's."""
    k = len(phis)
    G = np.array([[np.vdot(a, b) for b in phis] for a in phis])
    return float(np.real(_perm(G))) * d ** k / np.prod([d + i for i in range(k)])


def test_leaf_ratio_bound_and_tightness():
    rng = np.random.default_rng(9)
    d, k = 6, 4
    floor = np.prod([1 - t / d for t in range(k)])            # eq. 121
    for _ in range(50):
        phis = [v / np.linalg.norm(v) for v in rng.normal(size=(k, d)) + 1j * rng.normal(size=(k, d))]
        assert _leaf_ratio_exact(phis, d) >= floor - 1e-12
    ortho = list(np.eye(d)[:k].astype(complex))                 # orthonormal φ's hit the minimum exactly
    assert abs(_leaf_ratio_exact(ortho, d) - d ** k / np.prod([d + i for i in range(k)])) < 1e-12
    assert floor >= 1 - k * (k - 1) / (2 * d)


def test_leaf_ratio_exact_formula_matches_monte_carlo():
    rng = np.random.default_rng(10)
    d, k, T = 4, 3, 200000
    phis = [v / np.linalg.norm(v) for v in rng.normal(size=(k, d)) + 1j * rng.normal(size=(k, d))]
    psi = rng.normal(size=(T, d)) + 1j * rng.normal(size=(T, d))
    psi /= np.linalg.norm(psi, axis=1, keepdims=True)
    vals = np.prod([d * np.abs(psi @ p.conj()) ** 2 for p in phis], axis=0)
    assert abs(vals.mean() - _leaf_ratio_exact(phis, d)) < 5 * vals.std() / np.sqrt(T)


def test_mixture_leaf_ratio_bound():
    """The mixed state version from the proof: for ρ_t = (1 - g_t)ψ + g_t I/d with the same ψ,
    E_ψ Π_t d⟨φ_t|ρ_t|φ_t⟩ ≥ 1 - k(k - 1)/(2d). Computed exactly by expanding over subsets of copies."""
    rng = np.random.default_rng(11)
    d, k = 5, 4
    for _ in range(20):
        phis = [v / np.linalg.norm(v) for v in rng.normal(size=(k, d)) + 1j * rng.normal(size=(k, d))]
        gs = rng.random(k)
        total = 0.0
        for S in itertools.product([0, 1], repeat=k):
            w = np.prod([(1 - gs[t]) if S[t] else gs[t] for t in range(k)])
            sub = [phis[t] for t in range(k) if S[t]]
            total += w * (_leaf_ratio_exact(sub, d) if sub else 1.0)
        assert total >= 1 - k * (k - 1) / (2 * d) - 1e-12


# The learner

def test_histogram_learner_noiseless_converges():
    rng = np.random.default_rng(12)
    g = cell.random_g(rng)
    theta = rng.random(5000)
    y = cell.purity_of(g(theta), 1024)
    pred = cell.histogram_fit(theta, y, 16, 1024)
    assert cell._rmse(pred, cell.purity_of(g(cell.GRID), 1024)) < 0.05


def test_random_g_range():
    rng = np.random.default_rng(13)
    for _ in range(100):
        v = cell.random_g(rng)(np.linspace(0, 1, 1001))
        assert v.min() >= -1e-12 and v.max() <= 0.5 + 1e-12


def test_grouped_shared_label_unbiased():
    """Shared ψ0, G samples measured in one basis: the pooled guess is still right on average, because shots
    from different samples see the same ψ0."""
    rng = np.random.default_rng(14)
    d, g, G, k, T = 16, 0.3, 4, 3, 6000
    y = np.empty(T)
    for t in range(T):
        q = cell.dirichlet_flat(d, rng)
        s = np.concatenate([cell._shots(q, g, k, rng) for _ in range(G)])
        y[t] = (d + 1) * cell._pairs_colliding(s) / (s.size * (s.size - 1) / 2) - 1
    assert abs(y.mean() - cell.purity_of(g, d)) < 4 * y.std() / np.sqrt(T)
