# Reader's Guide: The Carry–Dirichlet Bridge

**Stefano Alimonti** — March 2026

*This guide explains how binary carry arithmetic gives rise to a Dirichlet series that approximates $L(s, \chi_4)^2$, how the same object admits a canonical weighted first-return resolvent realization, and why Witt vectors provide the algebraic justification. It assumes familiarity with Dirichlet $L$-functions and basic $p$-adic concepts.*

---

## 1. The Setup in One Paragraph

When two $K$-bit integers are multiplied in binary and the product has exactly $D = 2K - 1$ bits (the "D-odd" condition), the carry chain forms a bridge: $c_0 = c_D = 0$. The cascade valuation scans this bridge from top to bottom, stopping at each nonzero carry position $\tau$. The sector ratio $R(K) = \sigma_{10}/\sigma_{00}$ is conjectured to converge to $-\pi = -4L(1, \chi_4)$ [P1]. This paper promotes the mixed stopping-time profile from its unweighted sum ($s = 0$, trivial character) to a full Dirichlet series in $s$, and finds that the limit is approximately $L(s, \chi_4)^2$.

---

## 2. The Stopping-Time Dirichlet Series

The cascade stopping-time decomposition [P1, §8.2a] expresses the sector contributions as sums over stopping depths $\tau$. Defining

$$u_{\mathrm{mix}}(\tau) = \omega \cdot u_{10}(\tau) - u_{00}(\tau)$$

where $\omega = N_{10}/N_{00}$ is the count ratio, the stopping-time weights form the mixed linear combination naturally associated with the sector normalization. Each stopping depth $\tau$ corresponds to an odd integer $n(\tau) = 2(\tau - \tau_0) + 1$.

**The Dirichlet series.** Introducing the complex variable $s$ and a Dirichlet character $\chi$:

$$\mu_\chi(K, s) = \sum_{\tau \geq \tau_0} u_{\mathrm{mix}}(\tau;\, K) \cdot \chi(n(\tau)) \cdot n(\tau)^{-s}$$

At $s = 0$ and $\chi = \chi_0$ (trivial character), this reduces to the stopping-time sum $\sum u_{\mathrm{mix}}(\tau;\, K)$. The series $\mu_\chi(K, s)$ is a "carry-Dirichlet series" that encodes the full frequency structure of the stopping-time distribution.

The same channel can be written as a weighted first-return resolvent

$$M_\chi(K,s)=\omega(K)\,x_0^\top(I-Q_{10})^{-1}r_{10,\chi}(s)-x_0^\top(I-Q_{00})^{-1}r_{00,\chi}(s),$$

where the reward vectors $r_{ab,\chi}(s)$ insert the factor $\chi(n(\tau))\,n(\tau)^{-s}$ directly into the stopping-time operator. On the tested exact and high-K profiles, $M_\chi(K,s)$ and $\mu_\chi(K,s)$ agree to machine precision.

---

## 3. The Witt Vector Bridge

**Why Witt vectors?** Carries in base-$p$ addition are not an ad hoc construction — they are the *defining structure* of the $p$-adic integers when expressed in Witt coordinates. The ring of truncated Witt vectors $W_n(\mathbb{F}_p)$ is isomorphic to $\mathbb{Z}/p^n\mathbb{Z}$, and the Witt addition law reproduces classical carry arithmetic exactly.

**The five-step bridge:**

1. **Carry-Witt bijection (Proposition 1).** The Witt correction vector $\delta_j = c_j - 2c_{j+1}$ bijects onto the classical carry vector. Values: $\delta_j \in \lbrace{}-2, -1, 0, +1\rbrace{}$.

2. **Operator conjugacy (Proposition 2).** The transition matrix on carry states $(c_j, c_{j+1})$ and the transition matrix on Witt states $\delta_j$ are conjugate: $T_\delta = P \cdot T_c \cdot P^{-1}$ via an explicit $4 \times 4$ permutation matrix.

3. **Character channel (Proposition 3).** A $\chi_4$-weighted projection $w^\top \cdot v(\tau)$ produces identical coefficients on both carry and Witt sides. After removing the stationary component, the transient decays at rate $1/2$ — the Diaconis-Fulman spectral gap.

4. **D-odd extension (§2.4).** The conjugacy and channel identity extend to the nonstationary D-odd setting (depth-dependent operators).

5. **Resolvent identity (Proposition 4).** The sum $\sum_\tau u(\tau)$ equals the carry-side sector contribution to machine precision.

The bridge justifies the Dirichlet series construction: the stopping-time weights $u(\tau)$ have a $p$-adic origin in the Witt ghost map, and the character channel naturally selects $\chi_4$.

---

## 4. The Main Result: μ ≈ L(s, χ₄)²

**Theorem 1 (Convergence, numerical).** The family $\lbrace{}\mu_{\chi}(K, s)\rbrace{}_K$ is numerically observed to converge uniformly on tested compact subsets of $\mathrm{Re}(s) > 1$ with geometric Cauchy contraction:

$$\sup_{s \in \mathcal{K}} |\mu(K, s) - \mu(K', s)| \leq C \rho^{\min(K, K')}$$

with $\rho \approx 0.5$-$0.6$, consistent with the Diaconis-Fulman spectral gap.

**Theorem 2 (L² formula, numerical).** For $\chi_4$:

$$\mu_{\chi_4}^{(\infty)}(s) \approx \tilde{c} \cdot (1 + 3^{-s}) \cdot L(s, \chi_4)^2$$

with holdout relative error $\approx 2.5\%$ and cross-$K$ stability below $10^{-5}$. The same corrected-L² law is recovered when the channel is computed from the canonical weighted-resolvent object, so the law is intrinsic to the operator formulation.

The correction factor $(1 + 3^{-s})$ removes exactly one power of the Euler factor at $p = 3$. The $L^2$ structure reflects the Dirichlet convolution arising from the product $x \cdot y$: in the character ring, the multiplicative structure of the two factors produces a convolution of $L$-functions.

---

## 5. Why $L^2$ and Not $L$?

The $L^2$ structure is natural from two perspectives:

**Algebraic:** The Dirichlet series associated to a product $N = p \cdot q$ of two random integers is a convolution of two Dirichlet series (one for each factor). Under the Dirichlet character $\chi$, the convolution becomes $L(s, \chi) * L(s, \chi) = L(s, \chi)^2$ (up to finitely many Euler factors).

**Analytic:** The correction factor $F(s) = \mu(s)/(\tilde{c} \cdot L(s, \chi_4))$ is meromorphic and zero-free in the tested strip $[0.3, 0.7] \times [0, 15]$ (Proposition 5). Its poles coincide with the zeros of $L(s, \chi_4)$, confirming that the second $L$-factor comes from the correction, not from the base series. The direct and canonical carry objects themselves remain zero-free in that box.

---

## 6. Cross-Character Extension

The construction extends beyond $\chi_4$:

| Character $\chi$ | Conductor | Corrector | Improvement |
|---|---|---|---|
| $\chi_4$ (mod 4) | 4 | $(1 + 3^{-s})$ | baseline |
| $\chi_5$ (mod 5) | 5 | $(1 - \chi_5(3) \cdot 3^{-s})^2$ | +46% |
| $\chi_8$ (mod 8) | 8 | $(1 - \chi_8(3) \cdot 3^{-s})^2$ | +71% |
| $\chi_3$ (mod 3) | 3 | inert ($\chi_3(3) = 0$) | outlier |

The corrector always involves $p_0 = 3$ (the first odd prime). The conductor-3 character is an outlier because $3 | 3$ makes the corrector trivial. The same $(p_0,k)$ pattern is recovered from the canonical weighted-resolvent scan, so it is not an artifact of the series presentation.

---

## 7. What This Does Not Claim

- The $L^2$ formula is a **numerical observation** with controlled error bounds, not a proof.
- No **zero-transfer theorem** from the carry chain to $L$-function zeros.
- No simple **gamma-like completion** of the canonical object produces an approximate functional equation on the tested mirrored grids.
- No result toward the **Riemann Hypothesis** or Generalized Riemann Hypothesis.
- The 2.5% residual is **genuine** — the formula is approximate, not exact.

---

## References

1. [P1] S. Alimonti, "π from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge," this series. doi:[10.5281/zenodo.18895611](https://doi.org/10.5281/zenodo.18895611) — [GitHub](https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral)
2. [E] S. Alimonti, "The Trace Anomaly of Binary Multiplication," this series. doi:[10.5281/zenodo.18895604](https://doi.org/10.5281/zenodo.18895604) — [GitHub](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly)
3. [A] S. Alimonti, "Spectral Theory of Carries in Positional Multiplication," this series. doi:[10.5281/zenodo.18895593](https://doi.org/10.5281/zenodo.18895593) — [GitHub](https://github.com/stefanoalimonti/carry-arithmetic-A-spectral-theory)

---

*CC BY 4.0*
