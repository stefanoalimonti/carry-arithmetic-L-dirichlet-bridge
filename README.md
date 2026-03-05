# carry-arithmetic-L-dirichlet-bridge

**The Carry–Dirichlet Bridge: Stopping-Time Series and L²(s,χ₄)**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

The stopping-time decomposition of D-odd binary multiplication, projected onto the Dirichlet character χ₄, yields a Dirichlet series μ_χ(K,s) that converges uniformly on compact subsets of Re(s) > 1. For χ₄ the limit satisfies

> μ_χ₄^(∞)(s) ≈ c̃ · (1 + 3⁻ˢ) · L(s, χ₄)²

with holdout relative error ≈ 2.5% and cross-K stability below 10⁻⁵.

Key results:
- **Propositions 1–4** (unconditional): exact algebraic bridge from Witt p-adic arithmetic through carry operators, character channel, and resolvent identity.
- **Theorem 1** (numerical): uniform Cauchy convergence with geometric contraction ρ ≈ 0.5–0.6.
- **Theorem 2** (numerical): corrected-L² law with Euler-tail equivalence.
- **Proposition 5** (numerical): the correction factor F = μ/(cL) is zero-free in the tested strip region [0.3, 0.7] × [0, 15].

## Status

- **Proved:** Propositions 1–4 (algebraic bridge, resolvent identity — exact for all tested K).
- **Numerical:** Theorems 1–2, Propositions 5–6 (uniform convergence, L² law, zero-free structure, spectral gap — strong numerical evidence, rigorous proof pending).
- **Open target:** exact closed form for the 2.5% residual; analytic derivation of the corrector exponent k; universal character law.

## Repository Structure

```
paper/carry_dirichlet_bridge.md               The paper
experiments/
  _shared.py                                  Shared utilities (profiles, L-functions, characters)
  data/                                       Bundled high-K profile data (E45, E162)
  L01_witt_carry_identity.py                  Witt–carry bijection (§2.1)
  L02_operator_intertwiner.py                 Operator conjugacy (§2.2)
  L03_character_channel.py                    χ₄ character channel (§2.3)
  L04_dodd_nonstationary.py                   D-odd nonstationary extension (§2.4)
  L05_analytic_resolvent.py                   Resolvent identity and R(K) (§2.5)
  L06_diaconis_fulman_bound.py                Spectral gap under D-odd conditioning (§2.6)
  L07_uniform_convergence.py                  Uniform Cauchy convergence (§3)
  L08_euler_tail_identity.py                  Euler-tail structure of F (§4.2)
  L09_residual_audit.py                       Bootstrap, ablation, L² equivalence (§4.3–4.5)
  L10_local_corrector.py                      Cross-character corrector scan (§6)
  L11_mechanism_controls.py                   A/B mechanism controls (§6.3)
```

## Reproduction

High-K profile data and shared utilities are in `experiments/data/` and `experiments/_shared.py`.

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

## Dependencies

- Python >= 3.8, NumPy, mpmath

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [P1] | Pi from Pure Arithmetic | [`carry-arithmetic-P1-pi-spectral`](https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral) |
| [E] | The Trace Anomaly of Binary Multiplication | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [H] | Carry Polynomials and the Partial Euler Product (Control) | [`carry-arithmetic-H-euler-control`](https://github.com/stefanoalimonti/carry-arithmetic-H-euler-control) |
| [Frobenius] | Frobenius Eigenvalues and Gauss Sums from Witt Carries | [`carry-frobenius`](https://github.com/stefanoalimonti/carry-frobenius) |

### Citation

```bibtex
@article{alimonti2026carry_dirichlet,
  author  = {Alimonti, Stefano},
  title   = {The Carry--Dirichlet Bridge: Stopping-Time Series and $L^2(s,\chi_4)$},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-L-dirichlet-bridge}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
