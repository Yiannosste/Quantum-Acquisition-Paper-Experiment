# Theory

Draft. These proofs haven't been checked line by line yet, so don't quote them as settled.

## 1. What the paper claims

From section 2.4 of the paper:

> Conjecture. Let $\mathcal{T}$ be a task family admitting an exponential *acquisition* separation (coherent
> quantum memory with entangling joint measurement, against single-copy measurement with classical recording) for
> a functional class $\mathcal{F}$. Then the separation does not transfer to a downstream that is non-trivial in
> the sense of Definition 1: there is no cross-instance generalisation guarantee, for a learner acting on the
> acquisition record, that inherits the exponential advantage in the sense of Definition 2. **The
> acquisition-memory advantage is confined to single-instance estimation.**

Definition 1 asks for "a statistical learner over a distribution $\mathcal{D}$ of distinct problem instances,
evaluated by expected risk (generalisation error) on instances not measured during acquisition".

Definition 2 counts the "total number of physical interactions with instances drawn from $\mathcal{D}$
(instances measured, times copies or channel uses consumed per instance)" needed to get expected risk at most ε
with probability at least 1 − δ. $N_C$ is the best count for agents that measure "each sample singly and record
the outcome classically", $N_Q$ the best count for agents that may hold several samples in quantum memory. The
guarantee "inherits" the advantage if $N_C/N_Q$ grows exponentially in n at fixed ε and δ.

## 2. How I read the definitions

- **"Sample" means copy.** One could read it as "instance", and then the SWAP agent would count as classical and
  there'd be no gap. The paper rules that out: it describes the losing side as "single-copy measurement" (2.4),
  says the separations are bounds on "copies or channel uses consumed per instance" (2.4), and says about this
  exact construction that "the single-copy learner pays the acquisition cost at every sampled θ" (3.6).
- **Where the quantum part ends.** Section 3.1 says a protocol that measures first and learns from the records
  after holds no quantum memory. If the line were drawn after the SWAP test, the memory learner of Chen, Cotler,
  Huang and Li (arXiv:2111.05881, Definition 4.16) would also have none, and the paper's own row of quantum
  results would be empty. So the line has to sit before the joint measurement.
- **Worst case over the unknown.** If the unknown curve g were fixed and known, an agent could just write down
  the answer with zero copies, for any task. "PAC-style" in Definition 1 means the agent has to work for every g
  in the class.
- **Budget.** N is the most copies an agent ever uses. An agent that only promises an average budget can be cut
  off early (see the remark after Proposition 1).
- **Error units.** I use RMSE ≤ ε, which means risk ≤ ε².

## 3. The task

A sample is a pair (θ, ψ). θ is uniform on [0, 1] and the agent sees it. ψ is a Haar random state on d = 2ⁿ
dimensions, new for every sample, and the agent never sees it. Each copy of the sample is

    ρ = (1 − g(θ)) |ψ⟩⟨ψ| + g(θ) I/d

where g goes from [0, 1] to [0, ½] and is L-Lipschitz. The target is the purity, which only depends on θ:

    f(θ) = Tr ρ² = φ(g(θ)),   φ(g) = (1 − g)² + (2g − g²)/d

Loss is squared error, averaged over θ.

Useful facts: φ'(g) = −2(1 − g)(1 − 1/d), so |φ'| ≤ 2 and f is 2L-Lipschitz. For g ≤ ½ and d ≥ 2 we also have
|φ'| ≥ ½. And since g ≤ ½, the purity is at least ¼, so absolute and relative error are within a factor 4 of
each other.

**Shared version.** Use one ψ₀ for every sample. Then θ ↦ ρ_θ is a fixed map with
‖ρ_θ − ρ_θ'‖₁ ≤ 2L|θ − θ'|, which doesn't get worse with n. That's the construction from section 3.6 of the
paper. Everything below works for it too.

## 4. Proposition 1: classical agents need about 2^(n/2) copies

Any agent that measures one copy at a time and reaches RMSE ≤ ε with probability ≥ 1 − δ, for every g, uses N
copies with

    N(N − 1) ≥ (1 − 2δ) d

as long as δ < ½ and either ε < 1/16 (any L > 0, both versions) or ε < 3/16 (fresh ψ version, any L ≥ 0).

**Proof.**

*Step 1: write the agent as a tree.* Same as Chen, Cotler, Huang and Li (arXiv:2111.05881, Definition 4.15 and
section 5). Each node either draws a new sample (the agent sees θ, and could even pick it, since the proof never
uses how θ is drawn) or measures one copy of a sample it already has, with any measurement that depends on
everything so far. It can switch between samples and stop whenever it wants. Split every measurement outcome into
rank-one pieces E_t = w_t d |φ_t⟩⟨φ_t| with w_t = Tr E_t / d (their Remark 4.19). Finer outcomes only help the
agent, so a bound for these trees covers everything else. Extra helper qubits on one copy still give a
measurement on that copy. Random choices are a mix of trees, and the bound holds for each tree.

One more case: an agent could measure a copy gently, measure another copy, then come back to the first one. Each
step still touches one copy, but the whole thing isn't a tree like above. It still works: the probability of
each leaf is then a positive sum of products of one-copy pieces, so it is a positive sum of the rank-one terms
from Step 2, and the same bound holds.

*Step 2: compare one leaf to pure noise.* Fix a leaf. Say sample j got k_j copies along the path, with
Σ k_j ≤ N. Compare with the made-up case where every copy is I/d. In that case the leaf has probability
(θ terms) × Π w_t. In the real case, because the ψ's are independent,

    P(leaf) / P_noise(leaf) = Π_j E_ψ Π_{t in sample j} [ (1 − g_j) X_t + g_j ],   X_t = d |⟨φ_t|ψ⟩|²

Multiply out the product for one sample over subsets S of its copies. The weights (1 − g)^|S| g^(k−|S|) add up
to 1. Each E Π_{t∈S} X_t is exactly the pure state case of Chen et al. (their eq. 120), and their Lemma 5.12 says
it is at least Π_{i<s} 1/(1 + i/d) ≥ 1 − s(s − 1)/(2d), with s = |S|. So each sample gives at least
1 − k_j(k_j − 1)/(2d), and

    P(leaf) / P_noise(leaf) ≥ 1 − Σ_j k_j(k_j − 1)/(2d) ≥ 1 − N(N − 1)/(2d)

(For the shared version, do the same expansion over all N copies at once. The weights still add up to 1.)

*Step 3: from leaves to the whole tree.* If every leaf has P ≥ (1 − η) Q, then the total variation distance
TV(P, Q) is at most η (their Lemma 5.4). So for any two curves g and g', TV ≤ N(N − 1)/d.

*Step 4: learning means telling two curves apart.* Take two curves g₀ and g₁ whose purities are Δ apart.

- For L > 0: g₀(θ) = cθ with c = min(L, ¼), and g₁ = g₀ + ¼. Both are allowed, both change with θ (so in the
  shared version every sample really is a different state), and since |φ'| ≥ ½ their purities differ by at
  least ⅛ everywhere. So Δ ≥ ⅛.
- For the fresh ψ version: g₀ = 0 and g₁ = ½. Then Δ = ¾(1 − 1/d) ≥ 3/8.

Guess "g₀" if the agent's answer is within ε of f₀. If the truth is g₀, this is right with probability ≥ 1 − δ.
If the truth is g₁, the answer is within ε of f₁ with probability ≥ 1 − δ, so it is at least Δ − ε > ε away
from f₀ when ε < Δ/2, and the guess is right again. Being right both ways with probability 1 − δ needs
TV ≥ 1 − 2δ. Put that together with Step 3 and you get N(N − 1) ≥ (1 − 2δ)d. ∎

**Remark.** If an agent only promises that on average it uses at most N̄ copies, stop it after N̄/η copies. It
still works with probability ≥ 1 − δ − η, so N̄ ≥ η√((1 − 2δ − 2η)d).

**Why it works.** Averaged over ψ, a single copy is exactly I/d, so it says nothing. Information only shows up
when copies of the same sample collide, and N copies give about N²/d collisions.

## 5. Proposition 2: the quantum agent's cost doesn't depend on n

The agent draws M samples, runs one SWAP test on two copies of each, and writes down y = ±1, where
P(+1) = (1 + Tr ρ²)/2, so the average of y is the purity. Then it splits [0, 1] into B equal bins and predicts
the mean of y in each bin (any fixed value for an empty bin), clipped to [1/d, 1]. Its average risk is

    E[risk] ≤ 8L²/B² + 5B/M

With B = max(1, ⌈(L²M)^(1/3)⌉) that's at most 18(L/M)^(2/3) + 13/M. By Markov's inequality, risk ≤ ε_R with
probability ≥ 1 − δ once M ≥ max{ L(36/(δε_R))^(3/2), 26/(δε_R) }. So the number of copies is

    N_Q ≤ 2M = O( L/(δε_R)^(3/2) + 1/(δε_R) )

with no n anywhere.

**Proof.** In a bin with m ≥ 1 samples, the error at θ splits into noise (mean of y minus mean of f over the
bin) plus bias (mean of f over the bin minus f(θ)). The bias is at most 2L/B. The noise has variance at most
1/m, since |y| ≤ 1 and the y's are independent. So the bin's squared error is at most 2/m + 8L²/B². An empty bin
costs at most 1. With m ~ Binomial(M, 1/B) we have E[1/m when m ≥ 1] ≤ E[2/(m + 1)] ≤ 2B/(M + 1) and
P(m = 0) ≤ e^(−M/B) ≤ B/M. Averaging over bins gives the bound. For the choice of B: if L²M ≥ 1 then
(L²M)^(1/3) ≤ B ≤ 2(L²M)^(1/3), and the two terms are at most 8(L/M)^(2/3) and 10(L/M)^(2/3). If not, B = 1 and
the bound is 8L² + 5/M < 13/M. ∎

The SWAP test here is a Bell measurement on each qubit pair, accepting when the number of singlets is even. The
tests check this against the actual projector.

## 6. Putting them together

At fixed ε, δ and L, N_C / N_Q ≥ c · 2^(n/2). That grows exponentially in n, which is what Definition 2 calls
inheriting. The learner predicts for samples it never measured, which is Definition 1. So both parts of the
paper's falsification test (section 3.6) are met.

The constant c is tiny, though. At ε = δ = 0.1 and L = 1, Proposition 2 asks for M ≈ 6.8 million, so
c ≈ 6.5 × 10⁻⁸. The proven ratio only passes 1 at n = 48. Against the quantum agent's real cost from the
simulation (about 2,900 copies), the proven classical limit wins from n ≈ 24. I don't know how the best
classical cost grows as ε shrinks; the proven limit √d has no ε in it.

## 7. Weak spots

1. **The gain is really a one-sample gain.** The generalisation part is just a 1D regression. The quantum edge
   comes from measuring each sample better, and Definition 2's copy counting carries it through. That matches
   the conjecture's last sentence ("confined to single-instance estimation") but breaks its formal statement.
   The paper has to keep one of the two. Changing the definitions now to rule this out is what section 2.4 of
   the paper itself calls circular.
2. **Linear targets are already done.** Huang, Kueng and Preskill (arXiv:2101.02464, Theorem 1) show no
   exponential gap for targets tr(Oρ(x)). Purity is quadratic, so it isn't covered. If you double the state to
   ρ ⊗ ρ and take O = SWAP, their classical model is exactly the SWAP agent here, so the two results agree.
3. **ψ hides a lot of parameters.** All the hardness comes from ψ, which has about 2^(n+1) real parameters. The
   conjecture doesn't limit that, and the dial θ is one number.
4. **Haar states can't be made efficiently.** The paper already counts a result on unpreparable states as valid
   (section 3.3, about Gyurik, Molteni and Dunjko). A version with easy to make states is open. Random stabilizer
   states probably aren't enough, because the proof uses moments of ψ up to about √d copies.
5. **Lemma 1, condition (a).** The SWAP test is the same for every sample, and the result still works. Lemma 1
   only gives conditions that are enough, so nothing breaks, but the "Limitation" box in section 3.5 is worded
   too strongly.

## 8. What would break it

- A one-copy-at-a-time agent that beats N(N − 1) ≥ (1 − 2δ)d. The test game in the simulation checks for this.
- A reading of "sample" or of where the quantum part ends that excludes the SWAP agent without also emptying the
  paper's own list of quantum results.
- Deciding that predicting a property from a visible setting θ isn't generalising across samples. But the
  classical results the paper puts in its populated cell (Huang et al. 2022, Lewis et al. 2024) have the same
  shape.

## References

- Chen, Cotler, Huang, Li. Exponential separations between learning with and without quantum memory.
  arXiv:2111.05881. Definitions 4.15 and 4.16, Remark 4.19, Lemma 5.4, Theorem 5.11, eq. 120 and 121, Lemma 5.12.
- Huang, Kueng, Preskill. Information-theoretic bounds on quantum advantage in machine learning.
  arXiv:2101.02464, Theorem 1. (The quote I used came from a text extraction of the PDF; check it against the
  paper before citing.)
- Montanaro, de Wolf. A survey of quantum property testing. SWAP test for purity.
- Brydges et al. Probing Rényi entanglement entropy via randomized measurements. Science 364, 260 (2019).
