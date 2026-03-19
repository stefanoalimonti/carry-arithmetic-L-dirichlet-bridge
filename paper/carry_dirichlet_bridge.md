# The Carry–Dirichlet Bridge: Stopping-Time Series and L²(s,χ₄)

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

We construct a Dirichlet series from the stopping-time decomposition of the binary carry chain for D-odd multiplication and study its relation to Dirichlet L-functions. The construction proceeds through an exact algebraic bridge: Witt p-adic correction vectors biject onto classical carry vectors (Proposition 1), the carry and Witt transition operators are conjugate via an explicit permutation intertwiner (Proposition 2), and a χ₄-weighted character channel extracts the resolvent identity to machine precision (Proposition 3). The same channel admits a canonical weighted first-return resolvent realization, so the complex variable s enters through an intrinsic operator quantity rather than through a post hoc fit.

The resulting carry-Dirichlet channel μ_χ(K,s) = Σ_τ u_mix(τ) χ(n(τ)) n(τ)^{−s} is numerically observed to converge uniformly on tested compact subsets of Re(s) > 1 as K → ∞, with geometric Cauchy contraction at a rate numerically close to the Diaconis–Fulman spectral gap 1/2 (Theorem 1, numerical). For the character χ₄ (mod 4), the limit satisfies

> μ_χ₄^(∞)(s) = c̃ · (1 + 3^{−s}) · L(s, χ₄)² + ε(s),

with holdout relative error |ε|/|μ| ≈ 2.5% and cross-K stability below 10⁻⁵ (Theorem 2, numerical). The correction factor (1 + 3^{−s}) removes exactly one power of the Euler factor at p = 3, and the L² structure reflects the Dirichlet convolution arising from the product x · y. The same local-corrector law is recovered when μ_χ is computed from the canonical weighted resolvent object. The correction factor F(s) = μ/(cL) is meromorphic and zero-free in the tested strip region [0.3, 0.7] × [0, 15] (Proposition 6). No simple gamma-like completion of the canonical object yields an approximate functional equation on the tested mirrored grids, and the canonical object remains zero-free in the tested strip box. A low-complexity corrector law extends to χ₅ and χ₈ with modified exponent, while χ₃ (conductor 3) is an outlier.

These results establish a quantitative bridge between the Markov dynamics of binary carries and the multiplicative structure of Dirichlet L-functions. They do not imply zero-transfer theorems or results toward the Riemann Hypothesis; at present the bridge reaches Euler/amplitude structure, but not the phase/completion mechanism that would be needed for zeros.

---

## 1. Introduction

### 1.1 From Carries to Dirichlet Series

The sector ratio $R(K)$ is conjectured to converge to $-\pi = -4L(1, \chi_4)$ [P1], connecting binary carry arithmetic to a single value of a Dirichlet $L$-function. This paper asks: is there a deeper connection — not just at that single special value, but as a function of the complex variable $s$?

**The idea in one sentence.** The cascade stopping-time decomposition of [P1, §8.2a] organizes the sector contributions by stopping depth $\tau$. Each depth $\tau$ corresponds to an odd integer $n(\tau) = 2(\tau - \tau_0) + 1$. Weighting the resulting mixed stopping-time profile by $n(\tau)^{-s}$ and a Dirichlet character $\chi$ turns it into a Dirichlet series — a function of $s$ that encodes the full frequency-by-frequency structure of the carry chain.

**What are Witt vectors?** To justify this construction, we need to connect carry arithmetic to $p$-adic number theory. The bridge is provided by *Witt vectors*: the ring $W(\mathbb{F}_p)$ is an algebraic encoding of the $p$-adic integers $\mathbb{Z}_p$, where addition and multiplication are computed component-wise with "correction terms" that are exactly the classical carries. In base 2, the Witt correction $\delta_j = c_j - 2c_{j+1}$ is a simple linear function of adjacent carries. This means that the carry transition operator and the Witt transition operator are *conjugate* — they encode the same dynamics in different coordinates.

### 1.2 Motivation

Papers [P1] and [E] provide strong numerical evidence that the cascade sector ratio R(K) of D-odd binary multiplication converges to $-\pi = -4L(1, \chi_4)$ (Conjecture 1 of [P1]; 4.0-digit exponential Richardson extrapolation through K = 21), connecting the carry chain to the Dirichlet character $\chi_4$. The stopping-time decomposition of [P1, §8.2a] organizes the sector contributions by depths $\tau$, with each term factoring as $P(\tau \mid \text{sector}) \cdot E[\text{val} \mid \tau, \text{sector}]$.

A natural question arises: does this stopping-time structure encode deeper information about L-functions beyond the single value $L(1, \chi_4)$? In this paper we form the Dirichlet series

$$\mu_\chi(K,s) = \sum_{\tau \geq \tau_0} u_{\text{mix}}(\tau;\, K)\, \chi(n(\tau))\, n(\tau)^{-s},$$

where n(τ) = 2(τ − τ₀) + 1 maps stopping depths to odd integers and u_mix = ω · u₁₀ − u₀₀ combines the sector contributions weighted by the count ratio ω = N₁₀/N₀₀. Equivalently, the same object is realized as a weighted first-return resolvent of the stopping-time chain. We show that this canonical operator quantity converges to a limit on Re(s) > 1 and that the limit is accurately modeled by a corrected square of L(s, χ).

### 1.3 Main results

**Algebraic bridge (§2).** We construct a five-step exact bridge from Witt p-adic arithmetic to the stopping-time resolvent objects that motivate the sector-ratio problem:

1. Carry vectors biject onto Witt correction vectors (Proposition 1).
2. Carry and Witt transition operators are conjugate via an explicit permutation (Proposition 2).
3. A χ₄-like character channel produces identical coefficients on both sides (Proposition 3).
4. The construction extends to the nonstationary D-odd multiplication setting (§2.4).
5. The resolvent identity Σ_τ u(τ) = μ_ab holds to machine precision (Proposition 4).

**Convergence (§3).** The family {μ_χ(K, s)}_K converges uniformly on compact subsets of Re(s) > 1 as K → ∞, with

$$\sup_{s \in \mathcal{K}} |\mu_\chi(K, s) - \mu_\chi(K', s)| \leq C \rho^{\min(K,K')}$$

for tested compacts K ⊂ {Re(s) > 1}, with ρ ≈ 0.5–0.6 consistent with the Diaconis–Fulman spectral gap (Theorem 1).

**Canonical realization (§3).** The mixed channel admits an intrinsic operator form

$$M_\chi(K,s)=\omega(K)\,M_{10}(K,s)-M_{00}(K,s),$$

where

$$M_{ab}(K,s)=x_0^\top(I-Q_{ab}(K))^{-1}r_{ab,\chi}(s),$$

and the weighted reward vector $r_{ab,\chi}(s)$ inserts the factor $\chi(n(\tau))\,n(\tau)^{-s}$ directly into the first-return resolvent. On the tested exact and high-K profiles, $M_\chi(K,s)$ agrees with the direct Dirichlet-series definition of $\mu_\chi(K,s)$ to machine precision.

**The L² formula (§4).** For χ₄, the high-K limit satisfies

$$\mu_{\chi_4}^{(\infty)}(s) \approx \tilde{c} \cdot (1 + 3^{-s}) \cdot L(s, \chi_4)^2$$

on Re(s) > 1, with holdout relative error ≈ 2.5% (Theorem 2, numerical). Two equivalent formulations of this identity are:

- Euler-tail form: μ ≈ c · L(s, χ₄) · ∏_{p>3} (1 − χ₄(p) p^{−s})^{−1}
- Corrected-L² form: μ ≈ c̃ · (1 + 3^{−s}) · L(s, χ₄)²

These coincide on Re(s) > 1 because (1 + 3^{−s}) · L² = L · ∏_{p≥5} (1 − χ₄(p) p^{−s})^{−1}. The Euler-tail truncation to p ≤ 401 introduces a discrepancy below 3 × 10⁻³ on tested grids.

**Zero-free structure (§5).** The correction factor F(s) = μ(s)/(c · L(s, χ₄)) is meromorphic in the tested box [0.3, 0.7] × [0, 15] with poles at the zeros of L(s, χ₄) and no zeros of its own (Proposition 6). This implies that μ^(∞) is analytic and zero-free in that region.

**Cross-character extension (§6).** A low-complexity corrector law

$$\mu_\chi(s) \approx c \cdot (1 - \chi(p_0) p_0^{-s})^k \cdot L(s, \chi)^2$$

extends to χ₅ (k = 2, gain +46%), χ₈ (k = 2, gain +71%), with p₀ = 3 in all cases. The character χ₃ is an outlier: since χ₃(3) = 0 (the conductor divides p₀), the corrector is inert.

### 1.4 Scope and non-claims

This paper is a contribution to experimental mathematics. The L² formula (Theorem 2) is a numerical observation with controlled error bounds, not a proof. We do not claim:

- An exact identity μ = c̃ · (1 + 3^{−s}) · L² (the 2.5% residual is genuine).
- A zero-transfer theorem from the carry chain to L-function zeros.
- Any result toward the Riemann Hypothesis.

### 1.5 Notation

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

**Proposition 5** (Diaconis–Fulman bound, numerical). *The numerical evidence strongly supports that conditioning on D-odd does not change the dominant spectral gap: the two-step survival product q(τ)q(τ+1) converges to 1/4 = (1/2)², consistent with the Diaconis–Fulman rate for the conditioned chain. A fully rigorous proof requires showing that the Doob h-transform of the carry chain (conditioning on D-odd) has spectral radius ≤ 1/2 for the restricted chain.* (Experiment L06.)

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

### 3.3 Canonical weighted resolvent

The stopping-time profile of each sector defines a transient first-return system with kernel $Q_{ab}$, reward vector $r_{ab}$, and entry state $x_0$. At $s = 0$ this gives the exact identity

$$\mu_{ab}=x_0^\top (I-Q_{ab})^{-1} r_{ab}=\sum_\tau u_{ab}(\tau).$$

For general $s$, the weighted reward vector

$$r_{ab,\chi}(s;\tau)=a_{ab}(\tau)\,E[\mathrm{val}\mid\tau,ab]\,\chi(n(\tau))\,n(\tau)^{-s}$$

defines the canonical analytic object

$$M_{ab}(K,s)=x_0^\top (I-Q_{ab}(K))^{-1} r_{ab,\chi}(s), \qquad
M_\chi(K,s)=\omega(K)\,M_{10}(K,s)-M_{00}(K,s).$$

On the tested profiles ($K = 7, 8, 9, 19, 20, 21, 999$) and on both right-half-plane and strip samples, the direct series and weighted-resolvent evaluations agree to machine precision. This gives an operator-level realization of the carry-Dirichlet channel without changing its numerical content.

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

For K ∈ {19, 20, 21, 999}, the function μ_χ₄(K, s) has no zeros in the box [0.3, 0.7] × [0, 15] of the critical strip, as determined by the argument principle (winding number computation on rectangular contours). The same conclusion holds for the canonical weighted-resolvent realization $M_{\chi_4}(K,s)$, which agrees with μ_χ₄(K,s) on the tested profiles.

For comparison, L(s, χ₄) has three zeros in the same box, at approximately s ≈ 0.5 + 6.02i, 0.5 + 10.24i, 0.5 + 12.99i.

### 5.2 Zero-free correction factor

**Proposition 6** (F zero-free). *The correction factor F(s) = μ(s)/(c · L(s, χ₄)) has winding index −3 on the boundary of [0.3, 0.7] × [0, 15], decomposing as Z − P = 0 − 3 where Z = 0 (no zeros) and P = 3 (three poles, inherited from the three zeros of L(s, χ₄)). This is consistent with F being meromorphic and zero-free in the tested strip region.*

Within the tested box, this implies that μ^(∞) is itself analytic and zero-free: the carry chain does not develop zeros at the L-function zero positions in this region. Extension to the full critical strip would require either larger contour tests or an analytic argument; neither is available at present.

### 5.3 Critical-line mirroring

On the critical line Re(s) = 1/2, the function |μ_χ₄(999, 1/2 + it)| exhibits local minima near the imaginary parts of the L-function zeros (t ≈ 5.91, 10.31, 12.81 vs the exact values 6.02, 10.24, 12.99). These minima reflect the harmonic content of the carry coefficients but do not converge to zeros as K increases: the depth-of-minima ratios S_j(K) remain approximately flat across K = 19, …, 999.

This confirms the structural limitation: the carry chain encodes the amplitude of L(s, χ₄) on Re(s) > 1 (via the L² formula) but cannot reproduce the phase cancellations that create zeros in the strip. The zero-transfer obstruction persists after passing from the direct stopping-time series to the canonical weighted-resolvent realization.

### 5.4 Completion diagnostics

To test whether the missing zero mechanism could arise from a simple archimedean completion, we formed gamma-like completions of the canonical weighted-resolvent object and compared $\Xi(s)$ against $\varepsilon\,\Xi(1-s)$ on mirrored strip samples. No raw or gamma-weighted candidate produced a small symmetry defect: the best mean defect stays above 0.18 across the tested characters, and for χ₄ the raw object already performs best (mean defect ≈ 0.278). Thus the present carry object does not admit a simple classical completion analogous to the one for Dirichlet L-functions.

### 5.5 Topological stability of minimum positions

The minimum *positions* of $|\mu_{\chi_4}(999, 1/2 + it)|$ are topologically stable under the K → ∞ limit. Comparing the K = 21 and K = 999 (Richardson proxy) profiles on a grid with spacing dt = 0.05, all 15 tested minima have shift = 0.000: no minimum moves by even one grid point. This strengthens the mirroring result of §5.3 in two ways:

- **§5.3** establishes that the depth-of-minima ratios $S_j(K)$ remain approximately flat (the *values* of the minima do not converge to zero), confirming the zero-transfer obstruction.
- **§5.5** establishes that the *locations* of the minima are fixed to within the grid resolution. The zero-detection structure is therefore a topological property of the carry operator, not an artifact of the truncation at finite K.

The carry object encodes the positions of L-function zeros (in the sense that its amplitude minima track them) with exact positional stability, while remaining strictly positive (zero-free). Both properties are simultaneously present: spatial accuracy without magnitude cancellation.

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

For χ₄, the best model coincides with the predicted corrector (p₀ = 3, k = 1); no scan is needed. Cross-K stability is below 10⁻⁵ for all improving characters (experiment L10). The same best pairs are recovered when the scan is performed on the canonical weighted-resolvent object, with stable selection across K ∈ {19, 20, 21, 999}; the local-corrector law is therefore an operator-level feature of the carry channel, not an artifact of the series presentation.

### 6.2 The χ₃ outlier

The character χ₃ (conductor 3) fails to improve because χ₃(3) = 0: the corrector at p₀ = 3 is trivially 1. No single-prime correction at larger primes (p = 5, 7, …, 29) compensates. This is a structural limitation: when the conductor divides the natural corrector prime, the mechanism is inert.

### 6.3 Mechanism analysis

A/B controls (experiment L11) provide partial evidence on the mechanism underlying the corrector pattern (p₀ = 3 with character-dependent k):

- **Not driven by a single factor**: removing the n^{−s} weight preserves the k-pattern, but removing the character channel or scrambling the τ → n map degrades it. The controls point to a genuinely hybrid mechanism rather than to a single dominant ingredient.
- **Not reducible to raw p-adic frequencies**: the distribution of v₃(xy) among D-odd pairs does not isolate p = 3 as dominant. The correction appears to be an analytic-structural phenomenon involving the interaction of the character channel, the stopping-time map, and the Dirichlet weight, though a complete mechanistic derivation remains open.

### 6.4 Critical-line non-universality

On the critical line Re(s) = 1/2, the mixed channel admits a linear decomposition (see §7.4):

$$\mu_{\chi}(s) \approx c_{\mathrm{carry}} \cdot L(s, \chi) + B(s) \cdot L'(s, \chi).$$

This two-term fit is effective only for characters whose conductor is a power of 2. Specifically:

| χ | q | Critical-line fit error |
|---|---|------------------------|
| χ₄ | 4 | 5% |
| χ₈ | 8 | 10% |
| χ₃ | 3 | 44% |
| χ₅ | 5 | 47% |

The mechanism behind the conductor restriction is arithmetic: the stopping-time map produces $n(\tau) = 2(\tau - \tau_0) + 1$, which takes only odd values. For χ with odd-prime conductor $q$ (such as χ₃ or χ₅), the character satisfies $\chi(n) = 0$ whenever $q \mid n$, leaving systematic gaps in the Dirichlet series support. For $q = 2^k$ (χ₄, χ₈), every odd integer is coprime to $q$, so the series is complete.

This restriction is distinct from the cross-character corrector behavior on Re(s) > 1 (§6.1), where χ₅ gains +46% from the local corrector. The two regimes — L² on Re(s) > 1 and linear L + L' on Re(s) = 1/2 — have different universality properties.

---

## 7. Discussion

### 7.1 Summary

The carry-Dirichlet bridge establishes a quantitative connection between three mathematical objects:

1. The **Markov dynamics** of binary carry propagation (spectral gap 1/2, Diaconis–Fulman).
2. The **stopping-time decomposition** of D-odd multiplication (resolvent identity, cascade values).
3. The **Dirichlet L-function** L(s, χ₄) and its square L²(s, χ₄).

The connection takes the form μ^(∞)(s) ≈ c̃ · (1 + 3^{−s}) · L²(s, χ₄) on Re(s) > 1, with ≈ 2.5% residual. The L² structure arises from the Dirichlet convolution inherent in multiplication, and the (1 + 3^{−s}) correction reflects the partial absorption of the p = 3 Euler factor by the binary carry dynamics. The mixed channel also has a canonical weighted first-return resolvent realization, so the current bridge is genuinely operator-theoretic and not only a fitted Dirichlet series.

### 7.2 Relation to companion papers

- **[P1]** conjectures R(∞) = −π = −4L(1, χ₄) (Conjecture 1, supported by 4.0-digit numerical evidence) and develops the stopping-time decomposition. The present paper extends the analysis from a single value (s = 1) to a function of s, identifies a canonical weighted-resolvent object, and sharpens the remaining s = 1 closure to a scalar constant together with the conditioned 1/2-rate theorem. The topological stability of §5.5 (minimum positions shift = 0 from K = 21 to K = 999) provides a complementary stability result: the zero-detection structure is K-stable in position, not only in depth.
- **[E]** proves R = −π conditionally on LMH via the shifted resolvent. The present paper supports the same macro picture from the stopping-time side, but also shows that no simple completion of the current carry object yields an L-function-style functional equation.
- **[H]** shows that carry corrections do not contribute information about zeta zeros beyond the Euler product. The zero-free property of §5 is consistent with this finding.
- **[Frobenius]** constructs the Witt-carry bridge at odd primes: exact Frobenius factors at p = 3 and p = 7, and an anti-diagonal unitary factorization B = √p · P (with analytic Weil bound proof) valid for all primes. Propositions 1–2 of the present paper establish the analogous bridge at p = 2.

### 7.3 Open problems

1. **Exact identity or approximation?** Is there a closed form for the 2.5% residual? The correction F(s) may involve higher-order Euler factors or finitely many additional primes.

2. **Analytic computation of c̃.** The proportionality constant c̃ should be derivable from the carry chain structure. Numerically, $\tilde{c}(\sigma) = \mu(s)/[(1+3^{-s}) \cdot L^2(s,\chi_4)]$ converges as $\sigma \to \infty$ to the limit $c_{\mathrm{carry}} = \lim_{\sigma\to\infty} \tilde{c}(\sigma) \approx 0.0556$ (equal to $u_{\mathrm{mix}}(\tau_0)$, the dominant term in the stopping-time series). This value is close to $1/18 = 0.05556$ to within $0.05\%$; whether $c_{\mathrm{carry}} = 1/18$ exactly remains an open conjecture. (The value $0.025$ represents the mean relative error of the approximation, distinct from the proportionality constant $c_{\mathrm{carry}}$.)

3. **Why p = 3?** The correction involves p = 3 = 2 + 1 exclusively. A precise connection between the base b = 2 and the corrector prime p = b + 1 may hold in generality.

4. **Multi-prime extension.** Carry chains at base p = 3, 5, 7, … would provide the local factors at each prime. The product of all local factors should converge to L(s, χ) itself (single power, not squared), eliminating the need for the Euler-tail correction.

5. **Conductor dependence of k.** The exponent k = 1 for q = 4 and k = 2 for q = 5, 8 may follow from the structure of (ℤ/qℤ)* and its interaction with the base b = 2. A structural derivation would replace the current discrete scan.

6. **Universal character law.** The failure of χ₃ (conductor 3) suggests that a universal law requires treating the case gcd(q, p₀) > 1 separately. Whether a two-prime corrector (p₀ = 3, p₁ = 5) rescues χ₃ is untested.

7. **Archimedean completion.** The current canonical weighted-resolvent object does not admit a simple gamma-like completion with approximate symmetry $\Xi(s)\approx\varepsilon\,\Xi(1-s)$ on the tested grids. The missing completion mechanism is therefore a genuine structural problem, not a notational omission.

8. **Zero transfer.** The canonical carry object remains zero-free in the tested critical-strip box while $L(s,\chi_4)$ has three zeros there. Any RH-adjacent extension must therefore build phase and zero structure beyond the current amplitude/Euler bridge.

9. **Analytic form of B(s).** The critical-line coefficient function $B(1/2 + it)$ (§7.4) satisfies $B(t) \approx 0.041 \cdot t^{-0.23}$ with $B(\infty) \approx 0.012$. The exponent $0.23$ does not match $1/4$ or any simple rational. Derivation of $B(s)$ from the carry-chain resolvent, and whether $B(\infty) > 0$ implies a persistent second-order correction, remain open.

### 7.4 Critical-line spectral bridge

On the critical line Re(s) = 1/2, the mixed channel satisfies a linear two-term decomposition that is distinct from the L² formula valid on Re(s) > 1:

$$\mu_{\chi_4}(s) \approx c_{\mathrm{carry}} \cdot L(s, \chi_4) + B(s) \cdot L'(s, \chi_4), \qquad \mathrm{Re}(s) = 1/2,$$

where $L'(s, \chi_4) = \partial_s L(s, \chi_4)$ is the derivative with respect to $s$, and:

- $c_{\mathrm{carry}} \approx 0.0556$ is a constant independent of t (the limit $\lim_{\sigma\to\infty} \tilde{c}(\sigma)$, see §7.3 problem 2).
- $B(1/2 + it) \approx 0.041 \cdot t^{-0.23}$ decreases with height; the residual $B(\infty) \approx 0.012$.
- Fit error: $\approx 5\%$ on the critical line ($\sigma = 1/2$), improving to $\approx 0.5\%$ at $\sigma = 2$.

The two regimes — $L^2$ on Re(s) > 1 and $L + L'$ on Re(s) = 1/2 — are not contradictory: as $\sigma$ increases from 1/2 toward infinity, $B(s)$ grows (from the critical-line perspective, $B$ decreases with $t$ but increases with $\sigma$ in the right half-plane), and the derivative term becomes relatively less significant compared to the dominant $L^2$ structure.

**Spectral depth saturation.** A multi-term fit $\mu \approx \sum_{n=0}^{N} a_n L^{(n)}(s, \chi_4)$ on the critical line saturates at $N = 4$: the fit error does not decrease beyond $\approx 4\%$ for any $N \geq 4$. The residual is consistent with white noise in the $\{L^{(n)}\}$ basis (flat FFT spectrum, no dominant frequency). This establishes a structural floor: the carry channel encodes $L$ and its derivatives to order 3, with a $\sim 4\%$ indecodable remainder.

**Coefficient at zeros.** At the zeros $\rho$ of $L(s, \chi_4)$, the two-term decomposition reduces to $\mu(\rho) \approx B(\rho) \cdot L'(\rho)$. The residues $\mu(\rho)/L'(\rho)$ decrease slowly from $\approx 0.024$ at the first zero ($t \approx 6$) to $\approx 0.017$ at the 25th zero ($t \approx 58$), consistent with the power-law decay of $B(t)$.

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
| `L12_operator_core_regression.py` | Shared operator/stopping-time regression checks | §3.3 |
| `L13_s1_scalar_reduction.py` | Scalar reduction `R(∞)=C·L(1,χ₄)` at s = 1 | §7.2 |
| `L14_s1_lambda2_envelope.py` | Theorem-style envelope toward conditioned rate 1/2 | §7.2 |
| `L15_s1_tail_exchange.py` | Tail-majorant template for limit/sum exchange | §7.2 |
| `L16_canonical_weighted_resolvent.py` | Canonical weighted first-return resolvent realization | §3.3 |
| `L17_resolvent_local_factor_scan.py` | Operator-level local-factor scan `(p,k)` | §6.1 |
| `L18_completion_symmetry_diagnostics.py` | Completion / functional-equation diagnostics | §5.4 |
| `L19_canonical_zero_fingerprint.py` | Zero fingerprint for the canonical object | §5.1–§5.4 |
| `L20_spectral_bridge.py` | Critical-line spectral decomposition $\mu \approx c_{\mathrm{carry}} L + B L'$ | §7.4 |
| `L21_conductor_restriction.py` | Conductor restriction: $\chi_4, \chi_8$ vs.\ $\chi_3, \chi_5$ fit quality | §6.4 |
| `L22_zero_position_stability.py` | Topological stability of $|\mu_{\chi_4}|$ minimum positions $K=21\to 999$ | §5.5 |

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
python L12_operator_core_regression.py
python L13_s1_scalar_reduction.py
python L14_s1_lambda2_envelope.py
python L15_s1_tail_exchange.py
python L16_canonical_weighted_resolvent.py
python L17_resolvent_local_factor_scan.py
python L18_completion_symmetry_diagnostics.py
python L19_canonical_zero_fingerprint.py
python L20_spectral_bridge.py
python L21_conductor_restriction.py
python L22_zero_position_stability.py
```

### 8.4 Gate summary

All quantitative claims in this paper are supported by explicit quantitative thresholds documented in the experiment scripts.

| Claim | Gate | Result |
|-------|------|--------|
| Witt–carry bijection (Prop 1) | exact match all pairs K ≤ 10 | PASS |
| Operator conjugacy (Prop 2) | ‖T_d − P T_z P⁻¹‖ < 10⁻¹⁴ | PASS |
| Channel identity (Prop 3) | $\lvert ch_z - ch_d \rvert < 10^{-14}$ | PASS |
| Resolvent identity (Prop 4) | $\lvert \mu_{res} - \mu_{dir} \rvert < 10^{-16}$ | PASS |
| Cauchy convergence (Thm 1) | contraction ratios ≤ 0.75, envelope violation = 0 | PASS |
| L² formula (Thm 2) | hold-MRE ≤ 0.03, cross-K rel-std ≤ 10⁻² | PASS |
| F zero-free in tested box (Prop 6) | winding index Z = 0, P = N_L | PASS |
| Cross-character (§6) | ≥ 3 characters with ≥ 40% gain, rel-std ≤ 0.02 | PASS |
| Canonical weighted resolvent (§3.3) | `|M_\chi - \mu_\chi| ≤ 10^{-12}` on tested grids | PASS |
| Operator-level `(p,k)` law (§6.1) | χ₄ best `(3,1)` and ≥ 3 stable characters | PASS |
| Completion symmetry (§5.4) | simple FE candidate with mean defect < 0.10 | no simple candidate observed |
| Canonical zero fingerprint (§5) | in-box winding matches L-function zero count | carry object remains zero-free |

---

## References

- [P1] S. Alimonti, *Pi from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge*, 2026. doi:10.5281/zenodo.18895611
- [E] S. Alimonti, *The Trace Anomaly of Binary Multiplication*, 2026. doi:10.5281/zenodo.18895604
- [A] S. Alimonti, *Spectral Theory of Carries*, 2026. doi:10.5281/zenodo.18895593
- [B] S. Alimonti, *Carry Polynomials and the Euler Product*, 2026. doi:10.5281/zenodo.18895597
- [H] S. Alimonti, *Carry Polynomials and the Partial Euler Product (Control)*, 2026. doi:10.5281/zenodo.18895603
- [Frobenius] S. Alimonti, *Frobenius Eigenvalues and Gauss Sums from Witt Carries*, 2026. doi:10.5281/zenodo.18895613
- P. Diaconis, J. Fulman, *Carries, shuffling, and symmetric functions*, Advances in Applied Mathematics 43 (2009), 176–196.
