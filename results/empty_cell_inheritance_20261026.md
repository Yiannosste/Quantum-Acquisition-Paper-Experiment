# Results (2026-10-26, seed 20261026)

Made by `scripts/run_inheritance.py`. Target: error on unseen samples (RMSE of the purity) at most 0.1 in at least 90% of runs. 30 practice runs per setting, 200 scored runs, 20000 runs per side in the test game.

N* is the fewest copies that pass, at that budget and the next one. The range is 5% to 95% from resampling the scored runs. M = samples, k = copies per sample, G = samples sharing a basis (shared version only), B = bins.

| agent | n | d | N* | range | M | k | G | B | pass rate |
|---|---|---|---|---|---|---|---|---|---|
| Q swap (fresh ψ) | 2 | 4 | 1448 | 1448 to 1448 | 724 | 2 | 1 | 6 | 0.93 |
| Q swap (fresh ψ) | 4 | 16 | 2048 | 2048 to 2048 | 1024 | 2 | 1 | 6 | 0.945 |
| Q swap (fresh ψ) | 6 | 64 | 2048 | 2048 to 2896 | 1024 | 2 | 1 | 8 | 0.92 |
| Q swap (fresh ψ) | 8 | 256 | 2048 | 2048 to 2896 | 1024 | 2 | 1 | 8 | 0.91 |
| Q swap (fresh ψ) | 10 | 1024 | 2896 | 2896 to 2896 | 1448 | 2 | 1 | 8 | 0.98 |
| Q swap (fresh ψ) | 12 | 4096 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.955 |
| Q swap (fresh ψ) | 14 | 16384 | 2896 | 2048 to 4096 | 1448 | 2 | 1 | 12 | 0.915 |
| Q swap (fresh ψ) | 16 | 65536 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.975 |
| C collision (fresh ψ) | 2 | 4 | 5792 | 5792 to 8192 | 512 | 11 | 1 | 6 | 0.925 |
| C collision (fresh ψ) | 4 | 16 | 16384 | 11586 to 16384 | 1024 | 16 | 1 | 6 | 0.965 |
| C collision (fresh ψ) | 6 | 64 | 23170 | 23170 to 23170 | 256 | 90 | 1 | 12 | 0.99 |
| C collision (fresh ψ) | 8 | 256 | 16384 | 16384 to 23170 | 64 | 256 | 1 | 8 | 0.92 |
| C collision (fresh ψ) | 10 | 1024 | 32768 | 32768 to 32768 | 64 | 512 | 1 | 8 | 0.975 |
| C collision (fresh ψ) | 12 | 4096 | 46340 | 46340 to 46340 | 32 | 1448 | 1 | 12 | 0.94 |
| C collision (fresh ψ) | 14 | 16384 | 92682 | 92682 to 92682 | 32 | 2896 | 1 | 8 | 0.95 |
| C collision (fresh ψ) | 16 | 65536 | 185364 | 185364 to 262144 | 16 | 11585 | 1 | 16 | 0.92 |
| C local-RM (fresh ψ) | 2 | 4 | 8192 | 8192 to 8192 | 1024 | 8 | 1 | 6 | 0.945 |
| C local-RM (fresh ψ) | 4 | 16 | 32768 | 32768 to 32768 | 1024 | 32 | 1 | 8 | 0.995 |
| C local-RM (fresh ψ) | 6 | 64 | 46340 | 46340 to 46340 | 1024 | 45 | 1 | 8 | 0.935 |
| C local-RM (fresh ψ) | 8 | 256 | 65536 | 65536 to 65536 | 256 | 256 | 1 | 8 | 0.94 |
| Q swap (shared ψ0) | 2 | 4 | 1448 | 1448 to 2048 | 724 | 2 | 1 | 6 | 0.92 |
| Q swap (shared ψ0) | 4 | 16 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.975 |
| Q swap (shared ψ0) | 6 | 64 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.97 |
| Q swap (shared ψ0) | 8 | 256 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 12 | 0.95 |
| Q swap (shared ψ0) | 10 | 1024 | 2896 | 2896 to 2896 | 1448 | 2 | 1 | 8 | 0.975 |
| Q swap (shared ψ0) | 12 | 4096 | 2048 | 2048 to 2896 | 1024 | 2 | 1 | 8 | 0.92 |
| Q swap (shared ψ0) | 14 | 16384 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 12 | 0.96 |
| C pooled collision (shared ψ0) | 2 | 4 | 8192 | 8192 to 8192 | 4096 | 2 | 4 | 8 | 0.975 |
| C pooled collision (shared ψ0) | 4 | 16 | 16384 | 11586 to 16384 | 4096 | 4 | 4 | 8 | 0.97 |
| C pooled collision (shared ψ0) | 6 | 64 | 23170 | 23170 to 23170 | 1024 | 22 | 4 | 8 | 0.97 |
| C pooled collision (shared ψ0) | 8 | 256 | 16384 | 16384 to 16384 | 1024 | 16 | 16 | 8 | 0.945 |
| C pooled collision (shared ψ0) | 10 | 1024 | 16384 | 16384 to 16384 | 1024 | 16 | 64 | 8 | 0.94 |
| C pooled collision (shared ψ0) | 12 | 4096 | 32768 | 23170 to 32768 | 1024 | 32 | 64 | 8 | 0.99 |
| C pooled collision (shared ψ0) | 14 | 16384 | 46340 | 46340 to 46340 | 1024 | 45 | 64 | 8 | 0.95 |
| Q swap (control: ψ=\|0…0⟩) | 2 | 4 | 1448 | 1448 to 2048 | 724 | 2 | 1 | 6 | 0.915 |
| Q swap (control: ψ=\|0…0⟩) | 4 | 16 | 2048 | 2048 to 2896 | 1024 | 2 | 1 | 8 | 0.92 |
| Q swap (control: ψ=\|0…0⟩) | 6 | 64 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.985 |
| Q swap (control: ψ=\|0…0⟩) | 8 | 256 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.975 |
| Q swap (control: ψ=\|0…0⟩) | 10 | 1024 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.97 |
| Q swap (control: ψ=\|0…0⟩) | 12 | 4096 | 2896 | 2048 to 2896 | 1448 | 2 | 1 | 8 | 0.975 |
| Q swap (control: ψ=\|0…0⟩) | 14 | 16384 | 2048 | 2048 to 2896 | 1024 | 2 | 1 | 8 | 0.9 |
| C comp-basis (control) | 2 | 4 | 512 | 512 to 512 | 64 | 8 | 1 | 6 | 0.95 |
| C comp-basis (control) | 4 | 16 | 1024 | 1024 to 1024 | 256 | 4 | 1 | 8 | 0.99 |
| C comp-basis (control) | 6 | 64 | 1024 | 724 to 1024 | 256 | 4 | 1 | 8 | 0.995 |
| C comp-basis (control) | 8 | 256 | 1024 | 724 to 1024 | 512 | 2 | 1 | 8 | 0.99 |
| C comp-basis (control) | 10 | 1024 | 1024 | 1024 to 1024 | 64 | 16 | 1 | 8 | 0.975 |
| C comp-basis (control) | 12 | 4096 | 1024 | 724 to 1024 | 128 | 8 | 1 | 8 | 0.98 |
| C comp-basis (control) | 14 | 16384 | 1024 | 724 to 1024 | 256 | 4 | 1 | 8 | 0.985 |

## How fast N* grows

Slope of log2 N* against n, using n ≥ 6. The proof says at least 0.5 for classical agents in the long run, and 0 for the SWAP agent.

| agent | slope | 90% range | points |
|---|---|---|---|
| Q swap (fresh ψ) | 0.057 | 0.000 to 0.098 | 6 |
| C collision (fresh ψ) | 0.329 | 0.184 to 0.436 | 6 |
| Q swap (shared ψ0) | -0.025 | -0.125 to 0.000 | 5 |
| C pooled collision (shared ψ0) | 0.150 | -0.109 to 0.359 | 5 |
| Q swap (control: ψ=\|0…0⟩) | -0.050 | -0.125 to 0.000 | 5 |
| C comp-basis (control) | 0.000 | -0.000 to 0.000 | 5 |

## Test game: g = 0 vs g = 1/2 on one sample

Success is the best a test on the number of matching shots can do, estimated from the runs. That estimate comes out a bit high, which only makes it harder to stay under the limit. The limit is (1 + N(N - 1)/d)/2, capped at 1, and it is proven for every one-copy-at-a-time agent. The SWAP rows show the test can fail.

| n | N | N/√d | success | limit | over the limit? |
|---|---|---|---|---|---|
| 4 | 2 | 0.50 | 0.522 | 0.562 | no |
| 4 | 2 | 0.50 | 0.521 | 0.562 | no |
| 4 | 4 | 1.00 | 0.568 | 0.875 | no |
| 4 | 8 | 2.00 | 0.626 | 1.000 | no |
| 4 | 16 | 4.00 | 0.705 | 1.000 | no |
| 4 | 2 (SWAP) | 0.50 | 0.676 | 0.562 | yes |
| 6 | 2 | 0.25 | 0.504 | 0.516 | no |
| 6 | 4 | 0.50 | 0.529 | 0.594 | no |
| 6 | 8 | 1.00 | 0.579 | 0.938 | no |
| 6 | 16 | 2.00 | 0.644 | 1.000 | no |
| 6 | 32 | 4.00 | 0.743 | 1.000 | no |
| 6 | 2 (SWAP) | 0.25 | 0.685 | 0.516 | yes |
| 8 | 4 | 0.25 | 0.510 | 0.523 | no |
| 8 | 8 | 0.50 | 0.534 | 0.609 | no |
| 8 | 16 | 1.00 | 0.579 | 0.969 | no |
| 8 | 32 | 2.00 | 0.650 | 1.000 | no |
| 8 | 64 | 4.00 | 0.765 | 1.000 | no |
| 8 | 2 (SWAP) | 0.12 | 0.687 | 0.504 | yes |
| 10 | 8 | 0.25 | 0.509 | 0.527 | no |
| 10 | 16 | 0.50 | 0.534 | 0.617 | no |
| 10 | 32 | 1.00 | 0.584 | 0.984 | no |
| 10 | 64 | 2.00 | 0.658 | 1.000 | no |
| 10 | 128 | 4.00 | 0.783 | 1.000 | no |
| 10 | 2 (SWAP) | 0.06 | 0.687 | 0.501 | yes |
| 12 | 16 | 0.25 | 0.512 | 0.529 | no |
| 12 | 32 | 0.50 | 0.538 | 0.621 | no |
| 12 | 64 | 1.00 | 0.585 | 0.992 | no |
| 12 | 128 | 2.00 | 0.659 | 1.000 | no |
| 12 | 256 | 4.00 | 0.791 | 1.000 | no |
| 12 | 2 (SWAP) | 0.03 | 0.687 | 0.500 | yes |
