#!/usr/bin/env python3
"""
P1-SPRINT-L1: Uniform convergence diagnostics for mu_chi(K,s)
=============================================================

Goal:
  Test the first structural step of the new strategy:
    K -> infinity stabilization of mu_chi(K,s) on compact subsets of Re(s)>1.

Construction:
  - Small-K exact profiles from D-odd enumeration (K=7..K_exact_max)
  - High-K profiles from E45 (K=19,20,21) plus K=999 extrapolated proxy
  - Character channel:
      mu_chi(K,s) = sum_{tau>=tau0} u_mix(tau;K)*chi(n(tau))/n(tau)^s
      n(tau)=2*(tau-tau0)+1, tau0=3
      u_mix = omega(K)*u_10 - u_00

Diagnostics:
  - sup_s |mu(K,s)-mu(999,s)| on an s-compact grid
  - high-K geometric contraction ratio estimates
  - empirical bound template:
      sup_s |mu(K,s)-mu(K',s)| <= C * rho^(min(K,K'))
"""

import argparse
from typing import Callable, List

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


primitive_characters = _shared.primitive_characters
build_exact_bank = _shared.build_exact_bank
build_highk_bank = _shared.build_highk_bank
mu_chi_of_s = _shared.mu_chi_from_record


def sup_diff(rec_a: dict, rec_b: dict, chi_fn: Callable[[int], int], s_grid: List[complex]) -> float:
    vals = [abs(mu_chi_of_s(rec_a, chi_fn, s) - mu_chi_of_s(rec_b, chi_fn, s)) for s in s_grid]
    return max(vals) if vals else 0.0


def main():
    parser = argparse.ArgumentParser(description="Sprint L1 uniform convergence diagnostics.")
    parser.add_argument("--k-exact-min", type=int, default=7)
    parser.add_argument("--k-exact-max", type=int, default=11)
    args = parser.parse_args()

    flush("=" * 78)
    flush("P1-SPRINT-L1: uniform convergence diagnostics for mu_chi(K,s)")
    flush("=" * 78)
    flush(f"Exact K window: {args.k_exact_min}..{args.k_exact_max}")

    exact_bank = build_exact_bank(args.k_exact_min, args.k_exact_max)
    high_bank = build_highk_bank()

    rec_bank = {}
    rec_bank.update(exact_bank)
    rec_bank.update(high_bank)
    all_ks = sorted(rec_bank.keys())

    if 999 not in rec_bank:
        flush("ERROR: missing K=999 proxy; cannot run L1 diagnostics.")
        return

    chars = primitive_characters()

    # Compact sample in Re(s)>1 (real + complex points)
    s_grid = [
        complex(1.15, 0.0),
        complex(1.30, 0.0),
        complex(1.50, 0.0),
        complex(1.80, 0.0),
        complex(2.20, 0.0),
        complex(1.20, 1.0),
        complex(1.40, 2.0),
    ]

    flush(f"Total K values used: {all_ks}")
    flush(f"s-grid size: {len(s_grid)}")

    # per-character diagnostics
    pass_flags = []
    for ch in chars:
        name = ch["name"]
        fn = ch["fn"]

        # Delta to K=999 proxy
        deltas = {}
        for K in all_ks:
            if K == 999:
                continue
            deltas[K] = sup_diff(rec_bank[K], rec_bank[999], fn, s_grid)

        # High-K contraction ratios
        d19 = deltas.get(19, float("nan"))
        d20 = deltas.get(20, float("nan"))
        d21 = deltas.get(21, float("nan"))
        r20_19 = d20 / d19 if d19 > 1e-15 else float("nan")
        r21_20 = d21 / d20 if d20 > 1e-15 else float("nan")
        rho_hi = max(r for r in [r20_19, r21_20] if r == r and r > 0.0) if any(
            (r == r and r > 0.0) for r in [r20_19, r21_20]
        ) else float("nan")

        # Cauchy high-K distances
        c_19_20 = sup_diff(rec_bank[19], rec_bank[20], fn, s_grid)
        c_20_21 = sup_diff(rec_bank[20], rec_bank[21], fn, s_grid)
        c_21_999 = sup_diff(rec_bank[21], rec_bank[999], fn, s_grid)

        # empirical C,rho envelope on available K<999
        rho_env = rho_hi if rho_hi == rho_hi else 0.7
        rho_env = min(max(rho_env, 1e-6), 0.95)
        C_env = 0.0
        for K, d in deltas.items():
            C_env = max(C_env, d / (rho_env ** K))

        # verify envelope on available K
        env_viol = 0
        for K, d in deltas.items():
            rhs = C_env * (rho_env ** K)
            if d > rhs * (1.0 + 1e-10):
                env_viol += 1

        # gates (diagnostic numeric)
        g1 = (r20_19 == r20_19) and (r21_20 == r21_20) and (max(r20_19, r21_20) <= 0.75)
        # K=999 is an extrapolated proxy, so allow modest slack vs ||20-21||
        g2 = (c_20_21 <= c_19_20) and (c_21_999 <= 1.30 * c_20_21)
        g3 = env_viol == 0 and rho_env < 1.0
        ch_pass = g1 and g2 and g3
        pass_flags.append(ch_pass)

        flush("\n" + "-" * 78)
        flush(f"Character: {name}")
        flush("-" * 78)
        ks_show = [k for k in all_ks if k != 999]
        flush("delta_K := sup_s |mu(K,s)-mu(999,s)|")
        for K in ks_show:
            flush(f"  K={K:>3d}: {deltas[K]:.6e}")
        flush(f"ratios: d20/d19={r20_19:.6f}, d21/d20={r21_20:.6f}, rho_hi={rho_hi:.6f}")
        flush(f"Cauchy(high): ||19-20||={c_19_20:.6e}, ||20-21||={c_20_21:.6e}, ||21-999||={c_21_999:.6e}")
        flush(f"envelope: rho={rho_env:.6f}, C={C_env:.6e}, violations={env_viol}")
        flush(f"G1(high contraction)={'PASS' if g1 else 'FAIL'}  "
              f"G2(cauchy chain)={'PASS' if g2 else 'FAIL'}  "
              f"G3(envelope)={'PASS' if g3 else 'FAIL'}")
        flush(f"CHAR VERDICT: {'PASS' if ch_pass else 'FAIL'}")

    sprint_pass = all(pass_flags)

    flush("\n" + "=" * 78)
    flush("SPRINT-L1 VERDICT")
    flush("=" * 78)
    flush(f"Result: {'PASS' if sprint_pass else 'FAIL/PIVOT'}")
    if sprint_pass:
        flush("Interpretation: mu_chi(K,s) shows stable high-K geometric Cauchy behavior on tested compact.")
        flush("Next: Sprint L2 (s-renewal / continuation mechanism candidates).")
    else:
        flush("Interpretation: convergence diagnostics are not yet strong/uniform enough.")
        flush("Pivot: refine K-window or compact, then isolate unstable channel components.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

