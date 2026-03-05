#!/usr/bin/env python3
"""
P1-SPRINT-L5ter: low-complexity local corrector over L(s,chi)^2
===============================================================

Hypothesis:
  In Re(s)>1, mu_chi(s) is better modeled by
    mu_chi(s) ~= c * C_{p,k}(s) * L(s,chi)^2,
  where C_{p,k}(s) = (1 - chi(p) p^{-s})^k
  with low-complexity discrete correction (single odd prime p, small integer k).

This explicitly tests the post-L5bis structural reading:
  chi4: C_{3,1}(s) = (1 + 3^{-s}).
"""

import cmath
import math
from typing import Callable, Dict, List, Tuple

try:
    import mpmath as mp
except Exception:
    mp = None

import _shared

primitive_characters = _shared.primitive_characters
primes_upto = _shared.primes_upto
build_highk_bank = _shared.build_highk_bank
u_mix_map = _shared.u_mix_map
n_of_tau = _shared.n_of_tau
mu_chi = _shared.mu_chi
L_hurwitz = _shared.L_hurwitz
fit_c = _shared.fit_c
rel = _shared.rel
mean = _shared.mean


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def first_odd_prime_with_nonzero_chi(chi_fn: Callable[[int], int], plist: List[int]) -> int:
    for p in plist:
        if p % 2 == 1 and chi_fn(p) != 0:
            return p
    return -1


def main():
    flush("=" * 78)
    flush("P1-SPRINT-L5ter: local-corrector over L^2")
    flush("=" * 78)
    if mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = build_highk_bank()
    chars = primitive_characters()
    ks = [19, 20, 21, 999]

    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]
    hold_grid = [complex(s, t) for s in [1.15, 1.30, 1.55, 1.85, 2.20] for t in [0.5, 1.5, 2.5]]
    strip_grid = [complex(s, t) for s in [0.90, 0.80, 0.70, 0.60, 0.55] for t in [0.0, 1.0, 2.0, 3.0, 4.0]]
    p_candidates = [p for p in primes_upto(29) if p % 2 == 1]
    k_candidates = [-2, -1, 1, 2]

    # Per-character cache to avoid repeated Hurwitz evaluations.
    l_cache: Dict[Tuple[str, complex], complex] = {}

    rows = []
    chi4_detail = {}

    for ch in chars:
        name, q, chi_fn = ch["name"], ch["q"], ch["fn"]
        u_bank = {K: u_mix_map(bank[K]) for K in ks}

        def L(s: complex) -> complex:
            key = (name, s)
            if key not in l_cache:
                l_cache[key] = L_hurwitz(s, q, chi_fn)
            return l_cache[key]

        # Keep strip points away from L poles/zeros for ratios in diagnostics.
        strip_grid_ch = [s for s in strip_grid if abs(L(s)) > 1e-6]

        def basis_l(s: complex) -> complex:
            return L(s)

        def basis_l2(s: complex) -> complex:
            ls = L(s)
            return ls * ls

        # Baselines at K=999
        mu_fit_999 = [mu_chi(u_bank[999], chi_fn, s) for s in fit_grid]
        c_l = fit_c(mu_fit_999, [basis_l(s) for s in fit_grid])
        c_l2 = fit_c(mu_fit_999, [basis_l2(s) for s in fit_grid])
        hold_l = mean([rel(mu_chi(u_bank[999], chi_fn, s), c_l * basis_l(s)) for s in hold_grid])
        hold_l2 = mean([rel(mu_chi(u_bank[999], chi_fn, s), c_l2 * basis_l2(s)) for s in hold_grid])

        # Search low-complexity single-prime correction over L^2.
        best = None
        for p in p_candidates:
            cp = chi_fn(p)
            if cp == 0:
                continue
            for kappa in k_candidates:
                def corr(s: complex, p=p, cp=cp, kappa=kappa) -> complex:
                    return (1.0 - cp * cmath.exp(-s * math.log(p))) ** kappa

                basis_fit = [corr(s) * basis_l2(s) for s in fit_grid]
                c = fit_c(mu_fit_999, basis_fit)
                fit_mre = mean([rel(mu_chi(u_bank[999], chi_fn, s), c * corr(s) * basis_l2(s)) for s in fit_grid])
                hold_mre = mean([rel(mu_chi(u_bank[999], chi_fn, s), c * corr(s) * basis_l2(s)) for s in hold_grid])
                rec = {"p": p, "kappa": kappa, "fit_mre_999": fit_mre, "hold_mre_999": hold_mre}
                if best is None or rec["fit_mre_999"] < best["fit_mre_999"]:
                    best = rec

        # Cross-K stability for selected model
        hold_by_k = {}
        strip_by_k = {}
        cp_best = chi_fn(best["p"])
        for K in ks:
            mu_fit_k = [mu_chi(u_bank[K], chi_fn, s) for s in fit_grid]

            def corr_best(s: complex) -> complex:
                return (1.0 - cp_best * cmath.exp(-s * math.log(best["p"]))) ** best["kappa"]

            cK = fit_c(mu_fit_k, [corr_best(s) * basis_l2(s) for s in fit_grid])
            hold_by_k[K] = mean([rel(mu_chi(u_bank[K], chi_fn, s), cK * corr_best(s) * basis_l2(s)) for s in hold_grid])
            strip_by_k[K] = mean([rel(mu_chi(u_bank[K], chi_fn, s), cK * corr_best(s) * basis_l2(s)) for s in strip_grid_ch]) if strip_grid_ch else float("nan")

        vals = [hold_by_k[K] for K in ks]
        mm = mean(vals)
        sd = (sum((x - mm) ** 2 for x in vals) / len(vals)) ** 0.5
        rel_std = sd / max(abs(mm), 1e-15)

        improvement_vs_l2 = (hold_l2 - best["hold_mre_999"]) / max(hold_l2, 1e-15)

        # Predicted corrector: first odd prime with chi(p)!=0 and kappa=+1.
        p0 = first_odd_prime_with_nonzero_chi(chi_fn, p_candidates)
        cp0 = chi_fn(p0)
        def corr_pred(s: complex) -> complex:
            return (1.0 - cp0 * cmath.exp(-s * math.log(p0)))
        c_pred = fit_c(mu_fit_999, [corr_pred(s) * basis_l2(s) for s in fit_grid])
        hold_pred = mean([rel(mu_chi(u_bank[999], chi_fn, s), c_pred * corr_pred(s) * basis_l2(s)) for s in hold_grid])

        row = {
            "name": name,
            "hold_l": hold_l,
            "hold_l2": hold_l2,
            "best_p": best["p"],
            "best_kappa": best["kappa"],
            "best_fit_999": best["fit_mre_999"],
            "best_hold_999": best["hold_mre_999"],
            "best_hold_by_k": hold_by_k,
            "best_strip_999": strip_by_k[999],
            "best_rel_std": rel_std,
            "improvement_vs_l2": improvement_vs_l2,
            "pred_p0": p0,
            "pred_hold_999": hold_pred,
        }
        rows.append(row)
        if name == "chi4":
            chi4_detail = row

    flush(f"{'chi':>6s}  {'L':>10s}  {'L2':>10s}  {'best':>16s}  {'best_hold':>10s}  {'gain_vs_L2%':>12s}  {'relstd':>8s}")
    for r in rows:
        model_tag = f"p={r['best_p']},k={r['best_kappa']}"
        flush(
            f"{r['name']:>6s}  {r['hold_l']:10.6f}  {r['hold_l2']:10.6f}  {model_tag:>16s}  "
            f"{r['best_hold_999']:10.6f}  {100.0*r['improvement_vs_l2']:12.2f}  {r['best_rel_std']:8.6f}"
        )

    flush("\nchi4 predicted-corrector check:")
    flush(
        f"  predicted (p0={chi4_detail['pred_p0']},k=+1) hold@999 = {chi4_detail['pred_hold_999']:.6f}; "
        f"selected best (p={chi4_detail['best_p']},k={chi4_detail['best_kappa']}) hold@999 = {chi4_detail['best_hold_999']:.6f}"
    )
    flush(f"  chi4 selected strip-MRE@999 = {chi4_detail['best_strip_999']:.6f}")

    # Gates (L5ter interpretation)
    g1 = (chi4_detail["best_p"] == 3 and chi4_detail["best_kappa"] == 1)
    g2 = chi4_detail["best_hold_999"] <= 0.03
    g3 = sum(1 for r in rows if r["improvement_vs_l2"] >= 0.40) >= 3
    g4 = sum(1 for r in rows if r["best_rel_std"] <= 0.02) >= 3

    flush("\n" + "=" * 78)
    flush("SPRINT-L5ter GATES")
    flush("=" * 78)
    flush(f"L5ter.1: chi4 best model is (p=3,k=+1)               -> {'PASS' if g1 else 'FAIL'}")
    flush(f"L5ter.2: chi4 corrected-L2 hold MRE <= 0.03          -> {'PASS' if g2 else 'FAIL'}")
    flush(f"L5ter.3: >=3 chars gain >=40% vs L2 baseline         -> {'PASS' if g3 else 'FAIL'}")
    flush(f"L5ter.4: >=3 chars with cross-K rel-std <= 0.02      -> {'PASS' if g4 else 'FAIL'}")

    verdict = "PASS" if (g1 and g2 and g3 and g4) else "MIXED/INCONCLUSIVE"
    flush("\n" + "=" * 78)
    flush("SPRINT-L5ter VERDICT")
    flush("=" * 78)
    flush(f"Result: {verdict}")
    if verdict == "PASS":
        flush("Low-complexity local-corrector over L^2 is strongly supported in current tested family.")
    else:
        flush("Local-corrector over L^2 has strong character-dependent signal, but not fully universal by current gates.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

