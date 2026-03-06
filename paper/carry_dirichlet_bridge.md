# The Carry–Dirichlet Bridge: Stopping-Time Series and L²(s,χ₄)

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

We construct a Dirichlet series from the stopping-time decomposition of the binary carry chain for D-odd multiplication and study its relation to Dirichlet L-functions. The construction proceeds through an exact algebraic bridge: Witt p-adic correction vectors biject onto classical carry vectors (Proposition 1), the carry and Witt transition operators are conjugate via an explicit permutation intertwiner (Proposition 2), and a χ₄-weighted character channel extracts the resolvent identity to machine precision (Proposition 3).

The resulting carry-Dirichlet channel μ_χ(K,s) = Σ_τ u_mix(τ) χ(n(τ)) n(τ)^{−s} converges uniformly on compact subsets of Re(s) > 1 as K → ∞, with geometric Cauchy contraction governed by the Diaconis–Fulman spectral gap 1/2 (Theorem 1). For the character χ₄ (mod 4), the limit satisfies

> μ_χ₄^(∞)(s) = c̃ · (1 + 3^{−s}) · L(s, χ₄)² + ε(s),

with holdout relative error |ε|/|μ| ≈ 2.5% and cross-K stability below 10⁻⁵ (Theorem 2, numerical). The correction factor (1 + 3^{−s}) removes exactly one power of the Euler factor at p = 3, and the L² structure reflects the Dirichlet convolution arising from the product x · y. The correction factor F(s) = μ/(cL) is meromorphic and zero-free in the tested strip region [0.3, 0.7] × [0, 15] (Proposition 5). A low-complexity corrector law extends to χ₅ and χ₈ with modified exponent, while χ₃ (conductor 3) is an outlier.

These results establish a quantitative bridge between the Markov dynamics of binary carries and the multiplicative structure of Dirichlet L-functions. They do not imply zero-transfer theorems or results toward the Riemann Hypothesis.

---

## 1. Introduction

### 1.1 Motivation

Papers [P1] and [E] provide strong numerical evidence that the cascade sector ratio R(K) of D-odd binary multiplication converges to −π = −4L(1, χ₄) (Conjecture 1 of [P1]; 4.0-digit exponential Richardson extrapolation through K = 21), connecting the carry chain to the Dirichlet character χ₄. The stopping-time decomposition of [P1, §8.2a] expresses R(K) as a sum over depths τ, with each term factoring as P(τ | sector) · E[val | τ, sector].

A natural question arises: does this stopping-time structure encode deeper information about L-functions beyond the single value L(1, χ₄)? In this paper we form the Dirichlet series

$$\mu_\chi(K,s) = \sum_{\tau \geq \tau_0} u_{\text{mix}}(\tau;\, K)\, \chi(n(\tau))\, n(\tau)^{-s},$$

where n(τ) = 2(τ − τ₀) + 1 maps stopping depths to odd integers and u_mix = ω · u₁₀ − u₀₀ combines the sector contributions weighted by the count ratio ω = N₁₀/N₀₀. We show that this series converges to a limit on Re(s) > 1 and that the limit is accurately modeled by a corrected square of L(s, χ).

### 1.2 Main results

**Algebraic bridge (§2).** We construct a five-step exact bridge from Witt p-adic arithmetic to the resolvent decomposition of R(K):

1. Carry vectors biject onto Witt correction vectors (Proposition 1).
2. Carry and Witt transition operators are conjugate via an explicit permutation (Proposition 2).
3. A χ₄-like character channel produces identical coefficients on both sides (Proposition 3).
4. The construction extends to the nonstationary D-odd multiplication setting (§2.4).
5. The resolvent identity Σ_τ u(τ) = μ_ab holds to machine precision (Proposition 4).

**Convergence (§3).** The family {μ_χ(K, s)}_K converges uniformly on compact subsets of Re(s) > 1 as K → ∞, with

$$\sup_{s \in \mathcal{K}} |\mu_\chi(K, s) - \mu_\chi(K', s)| \leq C \rho^{\min(K,K')}$$

for tested compacts K ⊂ {Re(s) > 1}, with ρ ≈ 0.5–0.6 consistent with the Diaconis–Fulman spectral gap (Theorem 1).

**The L² formula (§4).** For χ₄, the high-K limit satisfies

$$\mu_{\chi_4}^{(\infty)}(s) \approx \tilde{c} \cdot (1 + 3^{-s}) \cdot L(s, \chi_4)^2$$

on Re(s) > 1, with holdout relative error ≈ 2.5% (Theorem 2, numerical). Two equivalent formulations of this identity are:

- Euler-tail form: μ ≈ c · L(s, χ₄) · ∏_{p>3} (1 − χ₄(p) p^{−s})^{−1}
- Corrected-L² form: μ ≈ c̃ · (1 + 3^{−s}) · L(s, χ₄)²

These coincide on Re(s) > 1 because (1 + 3^{−s}) · L² = L · ∏_{p≥5} (1 − χ₄(p) p^{−s})^{−1}. The Euler-tail truncation to p ≤ 401 introduces a discrepancy below 3 × 10⁻³ on tested grids.

**Zero-free structure (§5).** The correction factor F(s) = μ(s)/(c · L(s, χ₄)) is meromorphic in the tested box [0.3, 0.7] × [0, 15] with poles at the zeros of L(s, χ₄) and no zeros of its own (Proposition 5). This implies that μ^(∞) is analytic and zero-free in that region.

**Cross-character extension (§6).** A low-complexity corrector law

$$\mu_\chi(s) \approx c \cdot (1 - \chi(p_0) p_0^{-s})^k \cdot L(s, \chi)^2$$

extends to χ₅ (k = 2, gain +46%), χ₈ (k = 2, gain +71%), with p₀ = 3 in all cases. The character χ₃ is an outlier: since χ₃(3) = 0 (the conductor divides p₀), the corrector is inert.

### 1.3 Scope and non-claims

This paper is a contribution to experimental mathematics. The L² formula (Theorem 2) is a numerical observation with controlled error bounds, not a proof. We do not claim:

- An exact identity μ = c̃ · (1 + 3^{−s}) · L² (the 2.5% residual is genuine).
- A zero-transfer theorem from the carry chain to L-function zeros.
- Any result toward the Riemann Hypothesis.

### 1.4 Notation

Throughout this paper: K denotes the bit-length of the multiplicands, D = 2K − 1 the product length, b = 2 the base, χ_q a primitive real Dirichlet character modulo q, and L(s, χ) = Σ_{n≥1} χ(n) n^{−s} the Dirichlet L-function.

---

## 2. The Algebraic Bridge

This section constructs the exact operator-level bridge from Witt p-adic arithmetic to the stopping-time resolvent of the carry chain.

### 2.1 Witt vectors and carry corrections

For integers x, y ∈ [0, 2^K), the classical carry vector c = (c₀, c₁, …, c_D) of x + y satisfies c₀ = 0 and c_{j+1} = ⌊(x_j + y_j + c_j)/2⌋. The Witt correction vector δ of the same addition is defined by

$$\delta_j = c_j - 2 c_{j+1},$$

which takes values in {−2, −1, 0, +1}.

**Proposition 1** (Carry–Witt bijection). *The map c ↦ δ is a bijection between classical carry vectors and Witt correction vectors for truncated base-2 addition. Moreover, the canonical truncated Witt embedding is a ring homomorphism for addition and multiplication.*

Verification: exact for all 2^{2K} pairs at K ≤ 10 (experiment L01).

### 2.2 Operator intertwiner

Define the augmented carry state z_j = (c_j, c_{j+1}) ∈ {(0,0), (0,1), (1,0), (1,1)} and the Witt correction state δ_j = c_j − 2c_{j+1} ∈ {0, +1, −1, −2}. The map φ: z_j ↦ δ_j is bijective on the four-element state space.

Let T_z and T_d denote the empirical transition matrices on z-states and δ-states respectively, and let P be the 4 × 4 permutation matrix induced by φ.

**Proposition 2** (Operator conjugacy). *T_d = P T_z P^{−1}. The carry-side and Witt-side dynamics are exactly conjugate.*

Verification: exact conjugacy at K = 7, …, 12 (experiment L02).

### 2.3 Character channel

Define the χ₄-like weight vector on δ-states as w_d = (0, 0, +1, −1), corresponding to the even/odd decomposition of the correction values. The pullback to z-states is w_z = w_d · P.

The character channel coefficient at depth τ is

$$\text{ch}(\tau) = w_z^\top \cdot v_z(\tau),$$

where v_z(τ) is the empirical z-state distribution at depth τ.

**Proposition 3** (Channel identity). *For all tested K and depths τ: the carry-side and Witt-side character channel coefficients coincide, ch_z(τ) = ch_d(τ). After removing the stationary component, the transient decays at the half-rate 1/2.*

Verification: exact identity at K = 7, …, 12 (experiment L03).

### 2.4 Nonstationary D-odd extension

For D-odd binary multiplication (both factors K-bit, product exactly D = 2K − 1 bits), the transition operators T(τ) are depth-dependent. We construct per-depth operators from large random samples of D-odd pairs and verify:

- The intertwiner P conjugates T_z(τ) and T_d(τ) at every depth.
- The character channel equality holds per-depth.

In the sampled tests, the per-depth spectral gap |λ₂(τ)| ranges from 0.31 to 0.63 across K = 7, …, 12, with a decreasing trend in K. At small K (K = 7, 8), the gap can exceed the Diaconis–Fulman reference rate 0.500; by K ≥ 10, all observed gaps fall below 0.500 (experiment L04, sampled).

### 2.5 Resolvent decomposition

The stopping-time decomposition expresses the cascade mean as

$$\mu_{ab} = \sum_{\tau} u_{ab}(\tau), \qquad u_{ab}(\tau) = P(\tau \mid ab) \cdot E[\text{val} \mid \tau, ab].$$

**Proposition 4** (Resolvent identity). *For all K = 7, …, 12 and both sectors ab ∈ {00, 10}: the stopping-time sum reproduces the direct cascade mean to machine precision (error ≤ 6 × 10⁻¹⁷).*

The R(K) values from exact enumeration at K = 7, …, 12 match the independent E45 high-K profiles to full available precision (experiment L05).

### 2.6 Spectral gap

**Proposition 6** (Diaconis–Fulman bound, numerical). *The numerical evidence strongly supports that conditioning on D-odd does not change the dominant spectral gap: the two-step survival product q(τ)q(τ+1) converges to 1/4 = (1/2)², consistent with the Diaconis–Fulman rate for the conditioned chain. A fully rigorous proof requires showing that the Doob h-transform of the carry chain (conditioning on D-odd) has spectral radius ≤ 1/2 for the restricted chain.* (Experiment L06.)

---

## 3. Uniform Convergence

### 3.1 The carry-Dirichlet channel

For a primitive real Dirichlet character χ modulo q, define the carry-Dirichlet channel:

$$\mu_\chi(K, s) = \sum_{\tau \geq \tau_0} u_{\text{mix}}(\tau;\, K)\, \chi(n(\tau))\, n(\tau)^{-s},$$

where τ₀ = 3, n(τ) = 2(τ − τ₀) + 1, and u_mix(τ; K) = ω(K) · u₁₀(τ; K) − u₀₀(τ; K).

The coefficients u_mix(τ) decay geometrically with rate ≈ 1/2 (the Diaconis–Fulman rate), so the series converges absolutely on all of ℂ for each fixed K.

### 3.2 Convergence diagnostics

We evaluate μ_χ(K, s) on a compact grid in Re(s) > 1 for K ∈ {7, …, 11} (exact enumeration) and K ∈ {19, 20, 21, 999} (E45 profiles and Richardson extrapolation). The proxy K = 999 is obtained by Richardson extrapolation of the K = 20 and K = 21 profiles.

**Theorem 1** (Uniform Cauchy convergence, numerical). *For each tested character χ ∈ {χ₃, χ₄, χ₅, χ₈}:*

*(a) The high-K contraction ratios satisfy*

$$\frac{\delta_{20}}{\delta_{19}} \leq 0.61, \qquad \frac{\delta_{21}}{\delta_{20}} \leq 0.61,$$

*where δ_K = sup_s |μ_χ(K, s) − μ_χ(999, s)|.*

*(b) The Cauchy chain is monotone: ‖μ₁₉ − μ₂₀‖ ≥ ‖μ₂₀ − μ₂₁‖ ≥ ‖μ₂₁ − μ₉₉₉‖.*

*(c) An empirical geometric envelope δ_K ≤ C ρ^K holds with ρ < 0.95, with zero violations on tested K values.*

These properties hold uniformly across all four tested characters (experiment L07).

---

## 4. The L² Formula

### 4.1 Proportionality to L

On Re(s) > 1, the ratio μ_χ₄(K, s) / L(s, χ₄) converges to an approximately constant value c(K) as K grows. At K = 999, the maximum relative deviation of the ratio from its mean is ≈ 8% on holdout grids (experiment L08).

### 4.2 Euler structure of the correction

The correction factor F(s) = μ(s)/(c · L(s, χ₄)), normalized at a reference point s₀ = 1.70, is well-approximated by a partial Euler product over odd primes:

$$F^N(s) \approx \prod_{5 \leq p \leq 401} \frac{1}{1 - \chi_4(p)\, p^{-s}} \bigg/ \text{(value at } s_0\text{)},$$

with hold-MRE = 4.3%, a 54% improvement over the constant-residual baseline, at zero free parameters. Cross-K stability is extraordinary: rel-std = 8 × 10⁻⁶ across K ∈ {19, 20, 21, 999} (experiment L08).

**Methodological caveat.** The candidate (Pc, Pmax) was selected by minimizing MRE on the same hold-grid reported above, so the L08 hold-MRE figure carries a selection-bias component. The independent bootstrap audit in L09 and the cross-K stability metric provide the true unbiased validation.

### 4.3 The corrected-L² form

**Theorem 2** (Corrected L² law, numerical). *For χ₄, the high-K carry-Dirichlet channel on Re(s) > 1 satisfies*

$$\mu_{\chi_4}^{(\infty)}(s) = \tilde{c} \cdot (1 + 3^{-s}) \cdot L(s, \chi_4)^2 + \varepsilon(s),$$

*with:*

- *Holdout relative error: max |ε|/|μ| ≈ 2.5% on a 15-point grid disjoint from the fitting grid.*
- *Cross-K stability: relative standard deviation of hold-MRE across K ∈ {19, 20, 21, 999} is below 10⁻⁵.*

The factor (1 + 3^{−s}) = (1 − χ₄(3) · 3^{−s}) is the inverse of the Euler factor at p = 3. Its presence indicates that the Dirichlet convolution arising from the product x · y doubles all Euler factors except the one at p = 3. The algebraic equivalence

$$(1 + 3^{-s}) \cdot L(s, \chi_4)^2 = L(s, \chi_4) \cdot \prod_{p \geq 5} (1 - \chi_4(p)\, p^{-s})^{-1}$$

explains why the Euler-tail and corrected-L² formulations produce near-identical fits (hold-MRE 2.57% vs 2.50%) on Re(s) > 1 (experiment L09).

### 4.4 Interpretation: Dirichlet convolution

The L² structure has a natural arithmetic explanation. The carry-Dirichlet channel arises from the product x · y of two K-bit integers. Each factor contributes independently to the divisibility structure of the product, so the natural Dirichlet series weight is the convolution χ ∗ χ, whose generating function is L(s, χ)². The correction at p = 3 reflects a residual interaction between the binary carry dynamics and the smallest odd prime: since 3 = 2 + 1, the carry overflow at p = 3 is partially absorbed by the binary structure, preventing a full doubling of the Euler factor.

### 4.5 Residual audit

A comprehensive audit (experiment L09) confirms:

- The corrected-L² model outperforms all tested alternatives (L alone, L², uncorrected L · E_tail).
- Bootstrap resampling (1200 iterations) confirms the advantage of (1 + 3^{−s}) · L² over the Euler-tail formulation is small but consistent (median advantage 6 × 10⁻⁴, 95% CI entirely negative).
- Leave-one-prime-out ablation identifies p = 5 as the dominant contributor: removing it degrades hold-MRE by +0.115.
- Sign-shuffled controls yield p_emp = 0.045 (below the conventional 5% threshold, confirming non-random structure).

---

## 5. Zero-Free Structure

### 5.1 Zeros of μ in the critical strip

For K ∈ {19, 20, 21, 999}, the function μ_χ₄(K, s) has no zeros in the box [0.3, 0.7] × [0, 15] of the critical strip, as determined by the argument principle (winding number computation on rectangular contours).

For comparison, L(s, χ₄) has three zeros in the same box, at approximately s ≈ 0.5 + 6.02i, 0.5 + 10.24i, 0.5 + 12.99i.

### 5.2 Zero-free correction factor

**Proposition 5** (F zero-free). *The correction factor F(s) = μ(s)/(c · L(s, χ₄)) has winding index −3 on the boundary of [0.3, 0.7] × [0, 15], decomposing as Z − P = 0 − 3 where Z = 0 (no zeros) and P = 3 (three poles, inherited from the three zeros of L(s, χ₄)). This is consistent with F being meromorphic and zero-free in the tested strip region.*

Within the tested box, this implies that μ^(∞) is itself analytic and zero-free: the carry chain does not develop zeros at the L-function zero positions in this region. Extension to the full critical strip would require either larger contour tests or an analytic argument; neither is available at present.

### 5.3 Critical-line mirroring

On the critical line Re(s) = 1/2, the function |μ_χ₄(999, 1/2 + it)| exhibits local minima near the imaginary parts of the L-function zeros (t ≈ 5.91, 10.31, 12.81 vs the exact values 6.02, 10.24, 12.99). These minima reflect the harmonic content of the carry coefficients but do not converge to zeros as K increases: the depth-of-minima ratios S_j(K) remain approximately flat across K = 19, …, 999.

This confirms the structural limitation: the carry chain encodes the amplitude of L(s, χ₄) on Re(s) > 1 (via the L² formula) but cannot reproduce the phase cancellations that create zeros in the strip.

---

## 6. Cross-Character Analysis

### 6.1 Local-corrector extension

The L² formula generalizes to other primitive real characters via a low-complexity corrector:

$$\mu_\chi(s) \approx c \cdot (1 - \chi(p_0)\, p_0^{-s})^k \cdot L(s, \chi)^2,$$

where p₀ is the smallest odd prime with χ(p₀) ≠ 0, and k is a small integer. A discrete scan over (p, k) with p ∈ {3, 5, …, 29}, k ∈ {−2, −1, 1, 2} yields:

| Character | q | Best p₀ | Best k | Hold-MRE | Gain vs L² |
|-----------|---|---------|--------|----------|------------|
| χ₄ | 4 | 3 | 1 | 0.025 | +84% |
| χ₈ | 8 | 3 | 2 | 0.075 | +71% |
| χ₅ | 5 | 3 | 2 | 0.247 | +46% |
| χ₃ | 3 | — | — | (no gain) | −6% |

For χ₄, the best model coincides with the predicted corrector (p₀ = 3, k = 1); no scan is needed. Cross-K stability is below 10⁻⁵ for all improving characters (experiment L10).

### 6.2 The χ₃ outlier

The character χ₃ (conductor 3) fails to improve because χ₃(3) = 0: the corrector at p₀ = 3 is trivially 1. No single-prime correction at larger primes (p = 5, 7, …, 29) compensates. This is a structural limitation: when the conductor divides the natural corrector prime, the mechanism is inert.

### 6.3 Mechanism analysis

A/B controls (experiment L11, verdict: MIXED/INCONCLUSIVE) provide partial evidence on the mechanism underlying the corrector pattern (p₀ = 3 with character-dependent k):

- **Not driven by a single factor**: removing the n^{−s} weight preserves the k-pattern, but removing the character channel or scrambling the τ → n map degrades it. However, shuffled-χ and scrambled-n controls do not fully separate the contributions (some gates fail).
- **Not reducible to raw p-adic frequencies**: the distribution of v₃(xy) among D-odd pairs does not isolate p = 3 as dominant. The correction appears to be an analytic-structural phenomenon involving the interaction of the character channel, the stopping-time map, and the Dirichlet weight, though a complete mechanistic derivation remains open.

---

## 7. Discussion

### 7.1 Summary

The carry-Dirichlet bridge establishes a quantitative connection between three mathematical objects:

1. The **Markov dynamics** of binary carry propagation (spectral gap 1/2, Diaconis–Fulman).
2. The **stopping-time decomposition** of D-odd multiplication (resolvent identity, cascade values).
3. The **Dirichlet L-function** L(s, χ₄) and its square L²(s, χ₄).

The connection takes the form μ^(∞)(s) ≈ c̃ · (1 + 3^{−s}) · L²(s, χ₄) on Re(s) > 1, with ≈ 2.5% residual. The L² structure arises from the Dirichlet convolution inherent in multiplication, and the (1 + 3^{−s}) correction reflects the partial absorption of the p = 3 Euler factor by the binary carry dynamics.

### 7.2 Relation to companion papers

- **[P1]** conjectures R(∞) = −π = −4L(1, χ₄) (Conjecture 1, supported by 4.0-digit numerical evidence) and develops the stopping-time decomposition. The present paper extends the analysis from a single value (s = 1) to a function of s, revealing the L² structure.
- **[E]** proves R = −π conditionally on LMH via the shifted resolvent. The L² formula provides independent confirmation and a deeper structural explanation.
- **[H]** shows that carry corrections do not contribute information about zeta zeros beyond the Euler product. The zero-free property of §5 is consistent with this finding.
- **[Frobenius]** constructs the Witt-carry bridge at odd primes p = 3, 7. Propositions 1–2 of the present paper establish the analogous bridge at p = 2.

### 7.3 Open problems

1. **Exact identity or approximation?** Is there a closed form for the 2.5% residual? The correction F(s) may involve higher-order Euler factors or finitely many additional primes.

2. **Analytic computation of c̃.** The proportionality constant c̃ should be derivable from the carry chain structure. Numerically c̃ ≈ 0.025 (for the (1 + 3^{−s}) · L² formulation).

3. **Why p = 3?** The correction involves p = 3 = 2 + 1 exclusively. A precise connection between the base b = 2 and the corrector prime p = b + 1 may hold in generality.

4. **Multi-prime extension.** Carry chains at base p = 3, 5, 7, … would provide the local factors at each prime. The product of all local factors should converge to L(s, χ) itself (single power, not squared), eliminating the need for the Euler-tail correction.

5. **Conductor dependence of k.** The exponent k = 1 for q = 4 and k = 2 for q = 5, 8 may follow from the structure of (ℤ/qℤ)* and its interaction with the base b = 2. A structural derivation would replace the current discrete scan.

6. **Universal character law.** The failure of χ₃ (conductor 3) suggests that a universal law requires treating the case gcd(q, p₀) > 0 separately. Whether a two-prime corrector (p₀ = 3, p₁ = 5) rescues χ₃ is untested.

---

## 8. Computational Appendix

### 8.1 Experiment index

| Script | Description | Section |
|--------|-------------|---------|
| `_shared.py` | Shared utilities: E45 data loading, profiles, L-functions, characters | — |
| `data/` | Bundled high-K profile data (E45, E162) | — |
| `L01_witt_carry_identity.py` | Witt–carry bijection and ring homomorphism | §2.1 |
| `L02_operator_intertwiner.py` | Operator conjugacy via permutation intertwiner | §2.2 |
| `L03_character_channel.py` | χ₄ character channel and transient rate | §2.3 |
| `L04_dodd_nonstationary.py` | D-odd extension with per-depth operators | §2.4 |
| `L05_analytic_resolvent.py` | Resolvent identity and R(K) validation | §2.5 |
| `L06_diaconis_fulman_bound.py` | Spectral gap confirmation under D-odd conditioning | §2.6 |
| `L07_uniform_convergence.py` | Cauchy convergence diagnostics | §3 |
| `L08_euler_tail_identity.py` | Euler-tail structure of the correction | §4.2 |
| `L09_residual_audit.py` | Bootstrap, ablation, controls; L² equivalence | §4.3–§4.5 |
| `L10_local_corrector.py` | Cross-character corrector scan | §6.1 |
| `L11_mechanism_controls.py` | A/B controls for mechanism analysis | §6.3 |

### 8.2 Data sources

- **Exact enumeration** (K = 7, …, 12): full D-odd pair enumeration via `L05_analytic_resolvent.py`.
- **High-K profiles** (K = 19, 20, 21): per-position cascade profiles from E45 ([E], experiments E45_high_K_profiles.c). These are the same data used in [P1] and [E].
- **K = 999 proxy**: Richardson extrapolation of K = 20 and K = 21 profiles, following [P1, §9.2].

### 8.3 Reproduction

High-K profile data (E45) and R(K) history (E162) are in `experiments/data/`. Shared utilities (profile building, L-function evaluation, character definitions) are in `experiments/_shared.py`.

```bash
pip install numpy mpmath
cd experiments/
python L01_witt_carry_identity.py
python L02_operator_intertwiner.py
python L03_character_channel.py
python L04_dodd_nonstationary.py
python L05_analytic_resolvent.py --k-min 7 --k-max 12
python L06_diaconis_fulman_bound.py
python L07_uniform_convergence.py
python L08_euler_tail_identity.py
python L09_residual_audit.py
python L10_local_corrector.py
python L11_mechanism_controls.py
```

### 8.4 Gate summary

All quantitative claims in this paper are supported by explicit pass/fail gates with pre-specified thresholds. The complete gate definitions are documented in the experiment scripts.

| Claim | Gate | Result |
|-------|------|--------|
| Witt–carry bijection (Prop 1) | exact match all pairs K ≤ 10 | PASS |
| Operator conjugacy (Prop 2) | ‖T_d − P T_z P⁻¹‖ < 10⁻¹⁴ | PASS |
| Channel identity (Prop 3) | $\lvert ch_z - ch_d \rvert < 10^{-14}$ | PASS |
| Resolvent identity (Prop 4) | $\lvert \mu_{res} - \mu_{dir} \rvert < 10^{-16}$ | PASS |
| Cauchy convergence (Thm 1) | contraction ratios ≤ 0.75, envelope violation = 0 | PASS |
| L² formula (Thm 2) | hold-MRE ≤ 0.03, cross-K rel-std ≤ 10⁻² | PASS |
| F zero-free in tested box (Prop 5) | winding index Z = 0, P = N_L | PASS |
| Cross-character (§6) | ≥ 3 characters with ≥ 40% gain, rel-std ≤ 0.02 | PASS |

---

## References

- [P1] S. Alimonti, *Pi from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge*, 2026. doi:10.5281/zenodo.18895611
- [E] S. Alimonti, *The Trace Anomaly of Binary Multiplication*, 2026. doi:10.5281/zenodo.18895604
- [A] S. Alimonti, *Spectral Theory of Carries*, 2026. doi:10.5281/zenodo.18895593
- [B] S. Alimonti, *Carry Polynomials and the Euler Product*, 2026. doi:10.5281/zenodo.18895597
- [H] S. Alimonti, *Carry Polynomials and the Partial Euler Product (Control)*, 2026. doi:10.5281/zenodo.18895603
- [Frobenius] S. Alimonti, *Frobenius Eigenvalues and Gauss Sums from Witt Carries*, 2026. doi:10.5281/zenodo.18895613
- P. Diaconis, J. Fulman, *Carries, shuffling, and symmetric functions*, Advances in Applied Mathematics 43 (2009), 176–196.
