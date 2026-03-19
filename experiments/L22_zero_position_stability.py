#!/usr/bin/env python3
"""
L22: Topological stability of zero positions
=============================================

Verifies that the positions of local minima of |mu_chi4(1/2 + it)| are
stable as K -> infinity.

Method: compare the K = 21 and K = 999 profiles of |mu_chi4(1/2 + it)|
on a fine grid (dt = 0.05). For each of the 15 tested minima near the
known zeros of L(s, chi4), compute the positional shift between K = 21
and K = 999.

Expected result: shift = 0.000 on all minima (no minimum moves by even
one grid point), confirming that zero-location is a structural property
of the carry operator rather than an artifact of finite truncation.

References: Paper L §5.5.
"""

import sys
import os
import math
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _shared import build_highk_bank, mu_chi_from_record, primitive_characters, L_hurwitz


# Known zeros of L(s, chi4) on the critical line (imaginary parts)
ZEROS_L_CHI4 = [
    6.0209, 10.2438, 12.9881, 16.3426, 18.2920,
    21.4506, 23.2784, 25.7288, 28.3596, 29.6564,
    32.5922, 34.2000, 36.1429, 38.5119, 40.3227,
]


def find_local_min_near(mu_arr, ts, t_center, window=1.5):
    """Return (t_min, |mu|_min) for the local minimum of |mu| nearest to t_center."""
    mask = (ts >= t_center - window) & (ts <= t_center + window)
    if not np.any(mask):
        return float('nan'), float('nan')
    sub_t = ts[mask]
    sub_v = np.abs(mu_arr[mask])
    idx = np.argmin(sub_v)
    return sub_t[idx], sub_v[idx]


def main():
    bank = build_highk_bank()
    chi4_fn = {c["name"]: c for c in primitive_characters()}["chi4"]["fn"]

    sigma = 0.5
    dt = 0.05
    ts = np.arange(1.0, 55.0, dt)

    print("=" * 70)
    print("L22: TOPOLOGICAL STABILITY OF ZERO POSITIONS")
    print("=" * 70)
    print(f"\n  Grid: t in [1, 55), dt = {dt}")
    print(f"  Comparing K = 21 and K = 999")

    print("\n  Computing |mu| profiles for K = 21 and K = 999 ...")
    mu_21  = np.array([mu_chi_from_record(bank[21],  chi4_fn, complex(sigma, t)) for t in ts])
    mu_999 = np.array([mu_chi_from_record(bank[999], chi4_fn, complex(sigma, t)) for t in ts])

    print(f"\n  {'#':>3}  {'t_zero_L':>10}  {'t_min K=21':>12}  {'t_min K=999':>12}"
          f"  {'shift':>8}  {'|mu|_21':>10}  {'|mu|_999':>10}")
    print("  " + "-" * 78)

    shifts = []
    for i, tz in enumerate(ZEROS_L_CHI4):
        t21,  v21  = find_local_min_near(mu_21,  ts, tz)
        t999, v999 = find_local_min_near(mu_999, ts, tz)
        shift = abs(t21 - t999)
        shifts.append(shift)
        print(f"  {i+1:>3d}  {tz:>10.4f}  {t21:>12.4f}  {t999:>12.4f}"
              f"  {shift:>8.4f}  {v21:>10.6f}  {v999:>10.6f}")

    mean_shift = np.mean(shifts)
    max_shift  = np.max(shifts)
    print(f"\n  Mean positional shift (K=21 -> K=999): {mean_shift:.4f}")
    print(f"  Max  positional shift:                 {max_shift:.4f}")
    print(f"  Grid spacing dt:                        {dt}")
    print(f"  Shift in grid units: {mean_shift/dt:.1f} (mean),  {max_shift/dt:.1f} (max)")

    print("\n--- Summary ---")
    if max_shift < dt:
        print("  STABLE: no minimum moved by even one grid point.")
        print("  Zero positions are a structural property of the carry operator,")
        print("  not an artifact of the truncation at finite K.")
    else:
        print(f"  WARNING: max shift = {max_shift:.4f} > dt = {dt}")


if __name__ == "__main__":
    main()
