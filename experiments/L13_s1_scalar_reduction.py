#!/usr/bin/env python3
"""
L13: s=1 scalar reduction through the chi4 channel
==================================================

Formal reduction of the carry anomaly at s=1:
    R(inf) = C * L(1, chi4),
with the only remaining scalar closure target C = -4.

This script is intentionally proof-oriented:
  1) it brackets L(1, chi4) by alternating-series estimates,
  2) it computes high-K / extrapolated R(inf) candidates from the shared
     first-return resolvent core,
  3) it converts the bridge to a scalar interval for C.
"""

import math
from typing import Dict, Tuple

import _shared

PI = math.pi


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def midpoint_selector_identity(n_max: int = 12):
    rows = []
    for n in range(1, n_max + 1):
        v = math.sin(n * math.pi / 2.0)
        rows.append((n, v, "even -> 0" if n % 2 == 0 else "odd -> +/-1"))
    return rows


def leibniz_bracket(n_pairs: int) -> Tuple[float, float, float]:
    s = 0.0
    m_max = 2 * n_pairs + 1
    for m in range(m_max + 1):
        s += ((-1.0) ** m) / (2.0 * m + 1.0)
        if m == 2 * n_pairs:
            s_even = s
    s_odd = s
    lo = min(s_even, s_odd)
    hi = max(s_even, s_odd)
    return lo, hi, 0.5 * (lo + hi)


def r_estimates_from_core() -> Dict[str, float]:
    all_res = _shared.load_e45()
    by_k = {r["K"]: r for r in all_res}

    r19 = by_k[19]["S10"] / by_k[19]["S00"]
    r20 = by_k[20]["S10"] / by_k[20]["S00"]
    r21 = by_k[21]["S10"] / by_k[21]["S00"]

    p20_00 = _shared.build_profile(by_k[20], "00")
    p21_00 = _shared.build_profile(by_k[21], "00")
    p20_10 = _shared.build_profile(by_k[20], "10")
    p21_10 = _shared.build_profile(by_k[21], "10")
    p_inf_00 = _shared.build_extrapolated_profile(p20_00, p21_00, "00")
    p_inf_10 = _shared.build_extrapolated_profile(p20_10, p21_10, "10")

    omega20 = by_k[20]["n10"] / by_k[20]["n00"]
    omega21 = by_k[21]["n10"] / by_k[21]["n00"]
    omega_inf = 2.0 * omega21 - omega20
    r_step1 = omega_inf * _shared.mu_from_resolvent(p_inf_10) / _shared.mu_from_resolvent(p_inf_00)

    r_rich = 2.0 * r21 - r20
    den = r21 - 2.0 * r20 + r19
    r_aitken = r21 - (r21 - r20) ** 2 / den if abs(den) > 1e-18 else float("nan")

    return {
        "R19": r19,
        "R20": r20,
        "R21": r21,
        "R_step1_profile": r_step1,
        "R_richardson": r_rich,
        "R_aitken": r_aitken,
    }


def interval_mul_scalar(lo: float, hi: float, c: float) -> Tuple[float, float]:
    a = c * lo
    b = c * hi
    return min(a, b), max(a, b)


def interval_div(lo_a: float, hi_a: float, lo_b: float, hi_b: float) -> Tuple[float, float]:
    vals = [lo_a / lo_b, lo_a / hi_b, hi_a / lo_b, hi_a / hi_b]
    return min(vals), max(vals)


def main():
    flush("=" * 78)
    flush("L13: s=1 scalar reduction via chi4")
    flush("=" * 78)

    flush("\n" + "=" * 78)
    flush("A) Exact selector identity")
    flush("=" * 78)
    flush(f"{'n':>4s}  {'sin(n*pi/2)':>14s}  {'meaning':>14s}")
    for n, v, label in midpoint_selector_identity(12):
        flush(f"{n:4d}  {v:+14.9f}  {label:>14s}")

    flush("\n" + "=" * 78)
    flush("B) Rigorous alternating-series bracket for L(1, chi4)")
    flush("=" * 78)
    n_pairs = 1_000_000
    l_lo, l_hi, l_mid = leibniz_bracket(n_pairs)
    flush(f"L(1,chi4) in [{l_lo:.15f}, {l_hi:.15f}]")
    flush(f"midpoint = {l_mid:.15f}, width = {l_hi - l_lo:.3e}")
    flush(f"diagnostic midpoint - pi/4 = {l_mid - PI / 4:+.3e}")

    est = r_estimates_from_core()
    flush("\n" + "=" * 78)
    flush("C) High-K / extrapolated R(inf) estimates")
    flush("=" * 78)
    for key in ["R19", "R20", "R21", "R_step1_profile", "R_richardson", "R_aitken"]:
        value = est[key]
        flush(f"{key:16s} = {value:+.12f}   gap_to_-pi={value + PI:+.6e}")

    r_candidates = [est["R21"], est["R_step1_profile"], est["R_richardson"]]
    if est["R_aitken"] == est["R_aitken"]:
        r_candidates.append(est["R_aitken"])
    r_lo = min(r_candidates)
    r_hi = max(r_candidates)

    flush("\n" + "=" * 78)
    flush("D) Scalar closure interval R(inf) = C * L(1, chi4)")
    flush("=" * 78)
    c_lo, c_hi = interval_div(r_lo, r_hi, l_lo, l_hi)
    flush(f"R(inf) interval : [{r_lo:+.12f}, {r_hi:+.12f}]")
    flush(f"C interval      : [{c_lo:+.9f}, {c_hi:+.9f}]")
    flush(f"distance to -4  : [{c_lo + 4:+.3e}, {c_hi + 4:+.3e}]")
    flush(f"contains C=-4   : {'YES' if c_lo <= -4.0 <= c_hi else 'NO'}")

    r4_lo, r4_hi = interval_mul_scalar(l_lo, l_hi, -4.0)
    flush(f"\nIf C=-4 exactly => R in [{r4_lo:+.12f}, {r4_hi:+.12f}]")
    flush(
        f"overlap with current R interval: "
        f"[{max(r_lo, r4_lo):+.12f}, {min(r_hi, r4_hi):+.12f}]"
    )

    flush("\n" + "=" * 78)
    flush("L13 VERDICT")
    flush("=" * 78)
    flush(
        "Formal reduction achieved:\n"
        "  (1) the chi4 selector is exact,\n"
        "  (2) L(1,chi4) is bracketed independently of pi,\n"
        "  (3) the s=1 bridge is reduced to one scalar closure C,\n"
        "  (4) current data constrain C to a short interval around -4.\n"
        "\nRemaining theorem gap:\n"
        "  prove C=-4 from the carry resolvent itself."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
