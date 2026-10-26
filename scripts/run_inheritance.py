"""Runs the whole experiment and writes the results.

Part 1. For each agent and each n, find N*: the smallest number of copies at which the error on unseen samples
is at most EPS in at least 90% of runs. For every budget, the agent's settings (how many samples, how many copies
each, how many bins) are picked on a few practice runs and then scored on fresh runs. Only the classical agents
have settings to pick, so this helps them.

Part 2. A small test of the proven limit. With N one-copy-at-a-time shots on one sample, tell g = 0 apart from
g = 1/2. No classical strategy is allowed to beat (1 + N(N - 1)/d)/2. The SWAP agent should beat it with 2
copies, which shows the test can actually fail.

Writes results/empty_cell_inheritance_<date>.csv, .md and .png.
Run: python scripts/run_inheritance.py [--quick]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from emptycell import cell  # noqa: E402

SEED = 20261026
EPS, DELTA = 0.10, 0.10
ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

# (label, agent, how ψ is handed out, n values, settings to try). A setting is M (number of samples) or (M, G)
# with G samples per basis in the shared version. None means the SWAP agent: every sample gets 2 copies.
JOBS = [
    ("Q swap (fresh ψ)",            "swap",      "fresh",  range(2, 17, 2), None),
    ("C collision (fresh ψ)",       "collision", "fresh",  range(2, 17, 2), (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)),
    ("C local-RM (fresh ψ)",        "localrm",   "fresh",  range(2, 9, 2),  (16, 64, 256, 1024)),
    ("Q swap (shared ψ0)",          "swap",      "shared", range(2, 15, 2), None),
    ("C pooled collision (shared ψ0)", "collision", "shared", range(2, 15, 2),
     tuple((M, G) for M in (16, 64, 256, 1024, 4096) for G in (1, 4, 16, 64) if G * 8 <= M)),
    ("Q swap (control: ψ=|0…0⟩)",   "swap",      "known",  range(2, 15, 2), None),
    ("C comp-basis (control)",      "compbasis", "known",  range(2, 15, 2), (16, 32, 64, 128, 256, 512)),
]
BUDGETS = [int(2 * round(2 ** (j / 2) / 2)) for j in range(10, 41)]   # 32 up to 2^20, each step ×√2


def _eval(agent, nuisance, n, N, M, trials, rng, W, G=1):
    """Error on unseen samples for every run and every bin count (inf if there are more bins than samples)."""
    rows = np.full((trials, len(cell.BINS)), np.inf)
    for t in range(trials):
        r = cell.run_trial(agent, nuisance, n, N, M, rng, W=W, group=G)
        for j, B in enumerate(cell.BINS):
            rows[t, j] = r.get(B, np.inf)
    return rows


def find_nstar(job, n, pilot, evals, seed):
    label, agent, nuisance, _, Ms = job
    rng = np.random.default_rng(seed)
    W = cell._hamming_weights(n) if agent == "localrm" else None
    passes, records = 0, []
    t0 = time.time()
    for N in BUDGETS:
        settings = [(N // 2, 1)] if Ms is None else [c if isinstance(c, tuple) else (c, 1) for c in Ms]
        # each sample needs 2 copies (1 is fine when samples are pooled), and skip settings that are too slow
        grid = [(M, G) for M, G in settings if N // M >= (1 if G > 1 else 2)
                and (Ms is None or (M // G) * 2 ** n <= 2 ** 24)]
        best = None
        for M, G in grid:                                       # practice runs: keep the best setting
            rows = _eval(agent, nuisance, n, N, M, pilot, rng, W, G)
            usable = np.isfinite(rows).all(0)                   # skip bin counts that gave inf
            if not usable.any():
                continue
            q90 = np.where(usable, np.quantile(np.where(usable, rows, 0.0), 1 - DELTA, axis=0), np.inf)
            j = int(np.argmin(q90))
            if best is None or q90[j] < best[0]:
                best = (q90[j], M, G, cell.BINS[j])
        if best is None:
            continue
        _, M, G, B = best
        rm = _eval(agent, nuisance, n, N, M, evals, rng, W, G)[:, cell.BINS.index(B)]   # scored on fresh runs
        ok = rm <= EPS
        records.append(dict(N=N, M=M, k=N // M, G=G, B=B, success=float(ok.mean()),
                            rmse_median=float(np.median(rm)), rmse_q90=float(np.quantile(rm, 1 - DELTA)), ok=ok))
        passes = passes + 1 if ok.mean() >= 1 - DELTA else 0
        if passes == 3 or time.time() - t0 > 5400:
            break
    return label, n, records


def nstar_from(records, rng=None):
    """N* is the first budget where it and the next budget both pass. With rng, resample the runs first."""
    succ = []
    for r in records:
        ok = r["ok"] if rng is None else rng.choice(r["ok"], size=r["ok"].size)
        succ.append(ok.mean() >= 1 - DELTA)
    for i in range(len(succ) - 1):
        if succ[i] and succ[i + 1]:
            return records[i]["N"]
    return np.nan


# Part 2

def lower_bound_game(n, ratios, T, seed):
    rng = np.random.default_rng(seed)
    d = 2 ** n
    out = []
    for r in ratios:
        N = max(2, int(round(r * np.sqrt(d))))
        stats = []
        for g in (0.0, 0.5):
            stats.append(np.array([cell._pairs_colliding(cell._shots(cell.dirichlet_flat(d, rng), g, N, rng))
                                   for _ in range(T)]))
        hi = int(max(s.max() for s in stats)) + 1
        p0, p1 = (np.bincount(s.astype(int), minlength=hi) / T for s in stats)
        tv = 0.5 * np.abs(p0 - p1).sum()
        out.append(dict(n=n, N=N, N_over_sqrt_d=N / np.sqrt(d), success=(1 + tv) / 2,
                        bound=(1 + min(1.0, N * (N - 1) / d)) / 2))
    # the SWAP agent in the same game, with 2 copies
    f0, f1 = cell.purity_of(0.0, d), cell.purity_of(0.5, d)
    out.append(dict(n=n, N=2, N_over_sqrt_d=2 / np.sqrt(d), success=(1 + abs(f0 - f1) / 2) / 2,
                    bound=(1 + 2 / d) / 2, control="swap"))
    return out


# Main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--replot", metavar="STEM", help="redo the figure and the .md from saved CSVs, e.g. 20261026")
    a = ap.parse_args()
    if a.replot:
        return replot(a.replot)
    pilot, evals, T_game = (10, 40, 2000) if a.quick else (30, 200, 20000)
    ss = np.random.SeedSequence(SEED)
    tasks = [(job, n) for job in JOBS for n in job[3]]
    seeds = ss.spawn(len(tasks) + 8)
    t0 = time.time()
    with ProcessPoolExecutor() as ex:
        futs = [ex.submit(find_nstar, job, n, pilot, evals, seeds[i]) for i, (job, n) in enumerate(tasks)]
        game_ns = (4, 6, 8, 10, 12)
        gfuts = [ex.submit(lower_bound_game, n, (0.25, 0.5, 1, 2, 4), T_game, seeds[len(tasks) + i])
                 for i, n in enumerate(game_ns)]
        results = [f.result() for f in futs]
        game = [row for f in gfuts for row in f.result()]
    wall = time.time() - t0

    brng = np.random.default_rng(SEED + 1)
    summary = []
    for label, n, recs in results:
        ns = nstar_from(recs)
        boots = np.array([nstar_from(recs, brng) for _ in range(400)])
        lo, hi = np.nanquantile(boots, [0.05, 0.95]) if np.isfinite(boots).any() else (np.nan, np.nan)
        rec = next((r for r in recs if r["N"] == ns), None)
        summary.append(dict(agent=label, n=n, d=2 ** n, N_star=ns, ci_lo=lo, ci_hi=hi,
                            M=rec and rec["M"], k=rec and rec["k"], G=rec and rec["G"], B=rec and rec["B"],
                            success=rec and rec["success"], budgets_tried=len(recs)))

    stamp = dt.date.today().strftime("%Y%m%d")
    tag = "_quick" if a.quick else ""
    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / f"empty_cell_inheritance_{stamp}{tag}"
    with open(base.with_suffix(".csv"), "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# script scripts/run_inheritance.py, seed {SEED}, eps {EPS}, delta {DELTA}, "
                 f"pilot {pilot}, evals {evals}\n")
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    curves = [dict(agent=l, n=n, **{k: v for k, v in r.items() if k != "ok"}) for l, n, rs in results for r in rs]
    with open(OUT / f"empty_cell_inheritance_curves_{stamp}{tag}.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# script scripts/run_inheritance.py, seed {SEED}; every budget visited\n")
        w = csv.DictWriter(fh, fieldnames=list(curves[0]))
        w.writeheader()
        w.writerows(curves)
    with open(OUT / f"empty_cell_lower_bound_game_{stamp}{tag}.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# script scripts/run_inheritance.py, seed {SEED}, trials/hypothesis {T_game}\n")
        w = csv.DictWriter(fh, fieldnames=["n", "N", "N_over_sqrt_d", "success", "bound", "control"])
        w.writeheader()
        w.writerows([{**r, "control": r.get("control", "")} for r in game])

    fits = slope_fits(summary, brng)
    write_md(base.with_suffix(".md"), summary, game, fits, pilot, evals, T_game, wall)
    plot(base.with_suffix(".png"), summary, game)
    print(f"wrote {base}.* in {wall:.0f} s")


def slope_fits(summary, rng, n_min=6):
    """Slope of log2 N* against n for n ≥ n_min. The range comes from resampling the points, which is rough."""
    out = {}
    for label in dict.fromkeys(s["agent"] for s in summary):
        pts = [(s["n"], np.log2(s["N_star"])) for s in summary
               if s["agent"] == label and s["n"] >= n_min and np.isfinite(s["N_star"])]
        if len(pts) < 3:
            continue
        x, y = map(np.array, zip(*pts))
        b = np.polyfit(x, y, 1)[0]
        bs = []
        for _ in range(1000):
            i = rng.integers(len(x), size=len(x))
            if np.unique(x[i]).size >= 2:
                bs.append(np.polyfit(x[i], y[i], 1)[0])
        out[label] = (b, *np.quantile(bs, [0.05, 0.95]), len(pts))
    return out


def write_md(path, summary, game, fits, pilot, evals, T_game, wall=None):
    timing = f" Took {wall:.0f} s." if wall else ""
    L = [f"# Results ({dt.date.today()}, seed {SEED})", "",
         f"Made by `scripts/run_inheritance.py`. Target: error on unseen samples (RMSE of the purity) at most {EPS} "
         f"in at least {1 - DELTA:.0%} of runs. {pilot} practice runs per setting, {evals} scored runs, "
         f"{T_game} runs per side in the test game.{timing}", "",
         "N* is the fewest copies that pass, at that budget and the next one. The range is 5% to 95% from "
         "resampling the scored runs. M = samples, k = copies per sample, G = samples sharing a basis (shared "
         "version only), B = bins.", "",
         "| agent | n | d | N* | range | M | k | G | B | pass rate |", "|---|---|---|---|---|---|---|---|---|---|"]
    for s in summary:
        rate = round(s["success"], 3) if s["success"] is not None and s["success"] == s["success"] else ""
        L.append(f"| {s['agent'].replace('|', chr(92) + '|')} | {s['n']:.0f} | {s['d']:.0f} | {s['N_star']:.0f} | {s['ci_lo']:.0f} to "
                 f"{s['ci_hi']:.0f} | {s['M']:.0f} | {s['k']:.0f} | {s['G']:.0f} | {s['B']:.0f} | {rate} |")
    L += ["", "## How fast N* grows", "",
          "Slope of log2 N* against n, using n ≥ 6. The proof says at least 0.5 for classical agents in the long "
          "run, and 0 for the SWAP agent.", "",
          "| agent | slope | 90% range | points |", "|---|---|---|---|"]
    for k, (b, lo, hi, m) in fits.items():
        L.append(f"| {k.replace('|', chr(92) + '|')} | {b:.3f} | {lo:.3f} to {hi:.3f} | {m} |")
    L += ["", "## Test game: g = 0 vs g = 1/2 on one sample", "",
          "Success is the best a test on the number of matching shots can do, estimated from the runs. That "
          "estimate comes out a bit high, which only makes it harder to stay under the limit. The limit is "
          "(1 + N(N - 1)/d)/2, capped at 1, and it is proven for every one-copy-at-a-time agent. The SWAP rows "
          "show the test can fail.", "",
          "| n | N | N/√d | success | limit | over the limit? |", "|---|---|---|---|---|---|"]
    for r in game:
        tag = " (SWAP)" if r.get("control") else ""
        L.append(f"| {r['n']:.0f} | {r['N']:.0f}{tag} | {r['N_over_sqrt_d']:.2f} | {r['success']:.3f} | "
                 f"{r['bound']:.3f} | {'yes' if r['success'] > r['bound'] + 1e-9 else 'no'} |")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _read_csv(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    num = lambda v: float(v) if v not in ("", "None", "nan") else np.nan
    return [{k: (v if k in ("agent", "control") else num(v)) for k, v in r.items()} for r in rows]


def replot(stem):
    summary = _read_csv(OUT / f"empty_cell_inheritance_{stem}.csv")
    game = [{**r, "control": r["control"] or None} for r in _read_csv(OUT / f"empty_cell_lower_bound_game_{stem}.csv")]
    plot(OUT / f"empty_cell_inheritance_{stem}.png", summary, game)
    head = open(OUT / f"empty_cell_inheritance_{stem}.csv", encoding="utf-8").readline()
    pilot, evals = (int(head.split(f"{key} ")[1].split(",")[0]) for key in ("pilot", "evals"))
    T_game = int(open(OUT / f"empty_cell_lower_bound_game_{stem}.csv", encoding="utf-8").readline().split()[-1])
    write_md(OUT / f"empty_cell_inheritance_{stem}.md", summary, game,
             slope_fits(summary, np.random.default_rng(SEED + 1)), pilot, evals, T_game)
    print("redrew", stem)


def plot(path, summary, game):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    blue, orange, aqua = "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.2), facecolor="#fcfcfb")
    axs[1].sharey(axs[0])                                       # same y axis on the first three panels
    axs[2].sharey(axs[0])
    panels = [
        ("Fresh ψ per instance (main family)",
         [("Q swap (fresh ψ)", blue), ("C collision (fresh ψ)", orange), ("C local-RM (fresh ψ)", aqua)]),
        ("Shared ψ0 (θ ↦ ρ_θ Lipschitz, §3.6)",
         [("Q swap (shared ψ0)", blue), ("C pooled collision (shared ψ0)", orange)]),
        ("Control: ψ = |0…0⟩ known",
         [("Q swap (control: ψ=|0…0⟩)", blue), ("C comp-basis (control)", orange)]),
    ]
    for i, (ax, (title, series)) in enumerate(zip(axs[:3], panels)):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=grid, lw=0.8)
        for label, col in series:
            pts = [(s["n"], s["N_star"], s["ci_lo"], s["ci_hi"]) for s in summary
                   if s["agent"] == label and np.isfinite(s["N_star"])]
            if not pts:
                continue
            x, y, lo, hi = map(np.array, zip(*pts))
            ax.fill_between(x, lo, hi, color=col, alpha=0.15, lw=0)
            ax.plot(x, y, "-o", color=col, lw=2, ms=6, mec="#fcfcfb", mew=2, label=label)
            ax.annotate(label.split(" (")[0], (x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                        color=ink, fontsize=8, va="center")
        ax.set_yscale("log", base=2)
        ax.set_xlabel("qubits n")
        ax.set_xticks(range(2, 17, 2))
        ax.set_title(title, color=ink, fontsize=10, loc="left")
        ax.legend(frameon=False, fontsize=8, loc="upper left" if i == 2 else "lower right")
    axs[0].set_ylabel(f"total copies N* to held-out RMSE ≤ {EPS}")
    ax = axs[3]
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=grid, lw=0.8)
    xs = np.linspace(0.05, 4.5, 200)
    ax.plot(xs, (1 + np.minimum(1, xs ** 2)) / 2, color=muted, lw=2, ls="--", label="proved ceiling")
    pts = [r for r in game if not r.get("control")]
    ax.plot([r["N_over_sqrt_d"] for r in pts], [r["success"] for r in pts], "o", color=orange, ms=6,
            mec="#fcfcfb", mew=2, label="single-copy, n = 4…12")
    ctl = [r for r in game if r.get("control")]
    ax.plot([r["N_over_sqrt_d"] for r in ctl], [r["success"] for r in ctl], "s", color=blue, ms=7,
            mec="#fcfcfb", mew=2, label="SWAP, 2 copies")
    ax.set_xlabel("total copies N / √d")
    ax.set_ylabel("best success distinguishing g≡0 vs g≡½")
    ax.set_title("Lower-bound check", color=ink, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 0.95))
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=fig.get_facecolor())


if __name__ == "__main__":
    main()
