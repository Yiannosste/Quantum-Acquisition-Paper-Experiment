"""Break the code on purpose and check that the tests notice.

Passing tests don't mean much if the tests can't fail. Each change below plants one bug. The script prints
CAUGHT if the matching test fails, and MISSED if it doesn't.

Run: python scripts/mutation_check.py
"""
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from emptycell import cell  # noqa: E402

TESTS = str(ROOT / "tests" / "test_cell.py")


def wrong_collision(q, g, k, rng):
    """Uses d instead of d + 1 in the purity formula."""
    return q.size * cell._pairs_colliding(cell._shots(q, g, k, rng)) / (k * (k - 1) / 2) - 1


def wrong_swap(f, rng):
    """Accepts with probability Tr ρ² instead of (1 + Tr ρ²)/2."""
    return 2.0 * (rng.random(f.shape) < f) - 1.0


def leaky_shots(q, g, k, rng):
    """Every shot from ψ lands on outcome 0, so a single copy would give away information."""
    return np.where(rng.random(k) < g, rng.integers(q.size, size=k), 0)


BUGS = [
    ("shortcut uses a flat vector instead of a random one", "dirichlet_flat", lambda d, rng: np.full(d, 1 / d),
     "test_dirichlet_shortcut"),
    ("collision formula uses d instead of d + 1", "collision_label", wrong_collision, "test_collision_label_unbiased"),
    ("SWAP test accepts with the wrong probability", "swap_labels", wrong_swap, "test_swap_label_unbiased"),
    ("a single shot leaks information about psi", "_shots", leaky_shots, "test_single_shot"),
]

for name, attr, broken, test in BUGS:
    original = getattr(cell, attr)
    setattr(cell, attr, broken)
    failed = pytest.main(["-q", "-x", TESTS, "-k", test, "-p", "no:cacheprovider"]) != 0
    setattr(cell, attr, original)
    print(f"{'CAUGHT' if failed else 'MISSED'}: {name}")
