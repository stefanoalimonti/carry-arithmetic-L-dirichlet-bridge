# Experiments — Paper L

| Script | Description | Referenced in |
|--------|-------------|---------------|
| `L01_witt_carry_identity.py` | Witt p-adic correction ↔ classical carry bijection; ring homomorphism verification for truncated base-2 addition and multiplication. K ≤ 10. | §2.1 (Proposition 1) |
| `L02_operator_intertwiner.py` | Augmented carry-state and Witt-state transition operators; explicit permutation intertwiner P; exact conjugacy T_d = P T_z P⁻¹. K = 7…12. | §2.2 (Proposition 2) |
| `L03_character_channel.py` | χ₄-weighted character channel on z-states; carry/Witt channel identity; transient half-rate 1/2. K = 7…12. | §2.3 (Proposition 3) |
| `L04_dodd_nonstationary.py` | D-odd multiplication carry chains; per-depth transition operators; intertwiner and channel checks at every depth. K = 7…12. | §2.4 |
| `L05_analytic_resolvent.py` | Full D-odd enumeration; stopping-time decomposition; resolvent identity Σu(τ) = μ to machine precision; spectral gap; R(K) computation. K = 7…12. | §2.5 (Proposition 4) |
| `L06_diaconis_fulman_bound.py` | Survival probabilities, two-step products, K-convergence and τ-convergence diagnostics. Confirms dominant 1/2 rate under D-odd conditioning. | §2.6 (Proposition 5) |
| `L07_uniform_convergence.py` | sup-norm Cauchy diagnostics for μ_χ(K,s) on compact grids in Re(s) > 1. Small-K exact + high-K E45 profiles. Geometric contraction and envelope. | §3 (Theorem 1) |
| `L08_euler_tail_identity.py` | Normalized F_obs vs Euler-product templates on odd primes. Zero-parameter shape match. Cross-K stability. | §4.2 |
| `L09_residual_audit.py` | Full residual model audit: M0–M3 comparison, bootstrap M1 vs M2, leave-one-prime-out ablation, sign-shuffle control, cross-character survey, direct μ modeling (cL, cL², c(1+3⁻ˢ)L², cL·E_tail). | §4.3–§4.5 |
| `L10_local_corrector.py` | Low-complexity (p, k) scan for μ ~ c·C_{p,k}·L² across χ₃, χ₄, χ₅, χ₈. Predicted-corrector checkpoint for χ₄. Cross-K stability. | §6.1 |
| `L11_mechanism_controls.py` | A/B ablation controls: full, no_n_weight, no_channel, scrambled_n, shuffled_chi. Identifies hybrid mechanism (channel + mapping + analytic weight). | §6.3 |
| `L12_operator_core_regression.py` | Regression checks for the shared operator/stopping-time core: exact and parsed first-return identities, channel-helper equivalence, unified bank coverage. | §3.3 |
| `L13_s1_scalar_reduction.py` | Formal `s=1` reduction `R(∞)=C·L(1,χ₄)`: exact χ₄ selector, rigorous alternating-series bracket for `L(1,χ₄)`, scalar closure interval for `C`. | Supporting evidence for [P1], [E] |
| `L14_s1_lambda2_envelope.py` | Conditioned-chain dominant-rate envelope toward `1/2` using `q2(τ)` and high-K decrement ratios. | Supporting evidence for [P1], [E] |
| `L15_s1_tail_exchange.py` | Tail-majorant template for exchanging `K→∞` with the stopping-time sum. Produces explicit `epsilon(T)` targets. | Supporting evidence for [P1], [E] |
| `L16_canonical_weighted_resolvent.py` | Defines the canonical `s`-dependent carry object as a weighted first-return resolvent and verifies exact agreement with the direct mixed Dirichlet channel on tested grids. | §3.3 |
| `L17_resolvent_local_factor_scan.py` | Re-runs the low-complexity `(p,k)` local-factor scan using the canonical weighted resolvent object, with cross-K coherence checks across the character family. | §6.1 |
| `L18_completion_symmetry_diagnostics.py` | Tests raw and gamma-like completions of the canonical weighted resolvent object against `Ξ(s) ≈ εΞ(1-s)` symmetry on mirrored strip samples. | §5.4 |
| `L19_canonical_zero_fingerprint.py` | Critical-strip winding test for the canonical weighted resolvent object, compared directly against `L(s,χ₄)` and the correction factor `F(s)=M/(cL)`. | §5.1–§5.4 |

## Shared Utilities

| File | Description |
|------|-------------|
| `_shared.py` | Shared operator/stopping-time core: E45 data loading, exact/high-K profile building, first-return resolvent system, unified low-K/high-K bank, character definitions, L-function evaluation, fitting utilities, D-odd enumeration. Imported by L06–L19. |
| `data/E45_K19_K20.txt` | Per-position cascade profiles for K = 19 and K = 20 (from Paper [E]). |
| `data/E45_K21.txt` | Per-position cascade profile for K = 21 (from Paper [E]). |
| `data/E162_K19_analysis_output.txt` | R(K) = σ₁₀/σ₀₀ historical values (from Paper [E]). |

No external repository access is required.

## Requirements

Python >= 3.8, NumPy, mpmath.
