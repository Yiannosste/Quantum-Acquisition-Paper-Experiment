"""The task, the agents and the learner.

A sample is a dial setting θ (uniform on [0, 1], the agent sees it) plus a random state ψ (the agent never sees
it). Every copy of the sample is

    ρ = (1 - g(θ)) |ψ⟩⟨ψ| + g(θ) I/d

where the noise curve g goes from [0, 1] to [0, 1/2] and is unknown. The agent has to predict the purity
Tr ρ² = (1 - g)² + (2g - g²)/d for samples it never measured. The purity only depends on θ.

How ψ is handed out:
    "fresh"   new random ψ for every sample (the main version)
    "shared"  one random ψ0 for all samples
    "known"   ψ = |0...0⟩ and the agent knows it (control, no gap expected)

Agents. Each one spends N copies in total:
    swap       holds 2 copies and runs a SWAP test. Gives +1 with probability (1 + Tr ρ²)/2
    collision  one copy at a time. Measures k copies in one random basis and counts how often two outcomes match.
               This is the best one copy at a time method I know (it needs about √d copies). In the shared version
               it can also pool G samples into one basis.
    localrm    one copy at a time, random single-qubit bases (Brydges et al. 2019)
    compbasis  one copy at a time, only for the control. Measures in the computational basis
All agents use the same learner: average the labels in B equal bins of θ.

Shortcut used for speed: if ψ is random, the chance of outcome s in any fixed basis is (1 - g)·q_s + g/d, where
q is a uniformly random probability vector (normalised exponentials). So I never need a d×d unitary. The tests
check this against real random unitaries.
"""
from __future__ import annotations

import numpy as np

BINS = (1, 2, 3, 4, 6, 8, 12, 16)
GRID = (np.arange(512) + 0.5) / 512          # evenly spaced θ values used to compute the error on unseen samples


# The task

def purity_of(g, d: int):
    g = np.asarray(g, dtype=float)
    return (1 - g) ** 2 + (2 * g - g ** 2) / d


def random_g(rng: np.random.Generator):
    """A random smooth noise curve between 0 and 1/2 (three cosines, rescaled)."""
    a = rng.normal(size=3) / np.arange(1, 4)
    h = lambda t: sum(a[j] * np.cos((j + 1) * np.pi * np.asarray(t)) for j in range(3))
    scale = np.max(np.abs(h(np.linspace(0, 1, 2001))))
    return lambda t: np.clip(0.25 + 0.25 * h(t) / scale, 0.0, 0.5)     # clip only fixes rounding errors


def dirichlet_flat(d: int, rng: np.random.Generator) -> np.ndarray:
    """Outcome probabilities of a random pure state in any fixed basis."""
    e = rng.exponential(size=d)
    return e / e.sum()


# The agents

def swap_labels(f_vals: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One SWAP test per sample. The average label is the purity."""
    return 2.0 * (rng.random(f_vals.shape) < (1 + f_vals) / 2) - 1.0


def _shots(q: np.ndarray, g: float, k: int, rng: np.random.Generator) -> np.ndarray:
    d = q.size
    m_unif = rng.binomial(k, g)
    s = np.empty(k, dtype=np.int64)
    s[:m_unif] = rng.integers(d, size=m_unif)
    if k - m_unif:
        s[m_unif:] = np.minimum(np.searchsorted(np.cumsum(q), rng.random(k - m_unif)), d - 1)
    return s


def _pairs_colliding(s: np.ndarray) -> float:
    _, c = np.unique(s, return_counts=True)
    return float(np.sum(c * (c - 1)) / 2)


def collision_label(q: np.ndarray, g: float, k: int, rng: np.random.Generator) -> float:
    """Purity guess from k shots in one random basis. Right on average."""
    d = q.size
    s = _shots(q, g, k, rng)
    return (d + 1) * _pairs_colliding(s) / (k * (k - 1) / 2) - 1


def _hamming_weights(n: int) -> np.ndarray:
    s = np.arange(2 ** n)
    bits = (s[:, None] >> np.arange(n)) & 1
    D = (bits[:, None, :] != bits[None, :, :]).sum(-1)
    return (-2.0) ** (-D.astype(float))


def _haar_2x2(rng):
    z = (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) / np.sqrt(2)
    Q, R = np.linalg.qr(z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def localrm_label(n: int, g: float, k: int, W: np.ndarray, rng: np.random.Generator) -> float:
    """Brydges et al. 2019 purity guess from k shots in one random single-qubit basis per qubit."""
    d = 2 ** n
    v = rng.normal(size=d) + 1j * rng.normal(size=d)
    v = (v / np.linalg.norm(v)).reshape([2] * n)
    for q in range(n):                                  # random rotation on each qubit
        v = np.moveaxis(np.tensordot(_haar_2x2(rng), v, axes=([1], [q])), 0, q)
    amp2 = np.abs(v.reshape(d)) ** 2
    s = _shots(amp2 / amp2.sum(), g, k, rng)
    c = np.bincount(s, minlength=d).astype(float)       # bit order doesn't matter, Hamming distance is symmetric
    return d * (c @ W @ c - c.sum()) / (k * (k - 1))


def compbasis_label(d: int, g: float, k: int, rng: np.random.Generator) -> float:
    """Control only (ψ = |0...0⟩ is known). x counts shots that didn't land on 0...0. The formula below is
    right on average because purity = 1 - 2(αg) + (αg)²/α with α = 1 - 1/d."""
    alpha = 1 - 1 / d
    x = rng.binomial(k, alpha * g)
    return 1 - 2 * x / k + x * (x - 1) / (k * (k - 1) * alpha)


# The learner

def histogram_fit(theta: np.ndarray, y: np.ndarray, B: int, d: int, grid: np.ndarray = GRID) -> np.ndarray:
    """Average the labels in B equal bins of θ. An empty bin copies the closest bin that has data.
    Predictions are clipped to [1/d, 1], the range a purity can take."""
    b = np.minimum((theta * B).astype(int), B - 1)
    sums = np.bincount(b, weights=y, minlength=B)
    cnts = np.bincount(b, minlength=B)
    return _fill_and_predict(sums, cnts, B, d, grid)


def _fill_and_predict(sums, cnts, B, d, grid):
    filled = np.flatnonzero(cnts)
    if filled.size == 0:
        means = np.full(B, 0.5)
    else:
        means = np.empty(B)
        means[filled] = sums[filled] / cnts[filled]
        empty = np.flatnonzero(cnts == 0)
        if empty.size:
            nearest = filled[np.abs(empty[:, None] - filled[None, :]).argmin(1)]
            means[empty] = means[nearest]
    gb = np.minimum((grid * B).astype(int), B - 1)
    return np.clip(means[gb], 1 / d, 1.0)


# One run

def run_trial(agent: str, nuisance: str, n: int, N: int, M: int, rng: np.random.Generator,
              bins=BINS, W: np.ndarray | None = None, group: int = 1) -> dict:
    """Spend N copies on M training samples (N // M copies each) and return the error on unseen samples for
    each number of bins.

    `group` is only for the shared version: that many neighbouring samples get measured in one basis and their
    shots are pooled. Each group gets its own random basis, so its outcome probabilities are again a fresh random
    vector and `dirichlet_flat` is still exact.

    The error is computed exactly on GRID, which is the same as the average error on a new sample, because the
    purity doesn't depend on ψ.
    """
    d = 2 ** n
    g = random_g(rng)
    k = N // M
    theta = rng.random(M)
    gt = g(theta)
    f_grid = purity_of(g(GRID), d)
    bins = [B for B in bins if B <= M]
    out = {}
    if agent == "swap":
        if k != 2:
            raise ValueError("swap agent uses k = 2 (one SWAP test per sample)")
        y = swap_labels(purity_of(gt, d), rng)
        for B in bins:
            out[B] = _rmse(histogram_fit(theta, y, B, d), f_grid)
        return out
    if k < (1 if (agent == "collision" and nuisance == "shared" and group > 1) else 2):
        return {B: np.inf for B in bins}
    if agent == "collision" and nuisance == "fresh":
        y = np.array([collision_label(dirichlet_flat(d, rng), gt[i], k, rng) for i in range(M)])
    elif agent == "collision" and nuisance == "shared":
        order = np.argsort(theta)
        th, ys = [], []
        for start in range(0, M, group):
            idx = order[start:start + group]
            q = dirichlet_flat(d, rng)
            s = np.concatenate([_shots(q, gt[i], k, rng) for i in idx])
            th.append(theta[idx].mean())
            ys.append((d + 1) * _pairs_colliding(s) / (s.size * (s.size - 1) / 2) - 1)
        theta, y = np.array(th), np.array(ys)
        bins = [B for B in bins if B <= theta.size]
    elif agent == "localrm" and nuisance == "fresh":
        y = np.array([localrm_label(n, gt[i], k, W, rng) for i in range(M)])
    elif agent == "compbasis" and nuisance == "known":
        y = np.array([compbasis_label(d, gt[i], k, rng) for i in range(M)])
    else:
        raise ValueError(f"unsupported agent/nuisance pair: {agent}/{nuisance}")
    for B in bins:
        out[B] = _rmse(histogram_fit(theta, y, B, d), f_grid)
    return out


def _rmse(pred, truth) -> float:
    return float(np.sqrt(np.mean((pred - truth) ** 2)))
