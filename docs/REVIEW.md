# Review notes

Before writing this up I had the proofs reviewed as harshly as possible, like a referee who wants to reject the
paper. Every step of both proofs was re-derived and the cited lemmas were checked against the source. No step
broke. The review did find ten problems with how things were stated. Here they are, and what I did about each.

| # | Problem | What I changed |
|---|---|---|
| 1 | "Sample" in Definition 2 could be read as "instance", and then the SWAP agent counts as classical | THEORY section 2 now argues from the paper's own wording that it means copy |
| 2 | The two hard curves were constants, so in the shared version every sample was the same state | Added a pair that changes with θ (g₀ = cθ, g₁ = g₀ + ¼), which works for ε < 1/16 |
| 3 | Three things were assumed without saying so: worst case over g, worst case budget, and risk vs RMSE | All three are now stated in THEORY section 2, plus a remark for average budgets |
| 4 | The proof missed agents that measure a copy, keep it, and come back to it later | Added a paragraph to Step 1 showing the bound still holds |
| 5 | The conjecture was quoted without its last sentence | Now quoted in full, and weak spot 1 deals with that sentence |
| 6 | Missed Huang, Kueng and Preskill (arXiv:2101.02464), Theorem 1, which already covers linear targets | Added to the README and as weak spot 2. It turned out to be the best way to describe the result |
| 7 | "The ratio is about 2^(n/2)" oversold it, since the constant is tiny | THEORY section 6 now gives the real numbers and says where the proof kicks in |
| 8 | Wrong citation for the quantum side of the gap | Fixed |
| 9 | Small slips: a product that could go negative, B = 0 when L = 0, unset value for empty bins | Fixed |
| 10 | The list of weak spots was in the wrong order and missed 1, 2 and 6 | Reordered |

Still open after the review:

- A proper search for whether someone already did this. The quick search found nothing, but it was too small to
  count.
- Checking the Huang, Kueng and Preskill quote against the actual paper, since it came from a PDF text dump.
