#!/usr/bin/env python3
"""
L21: Conductor restriction of the critical-line fit
=====================================================

Tests whether the two-term linear approximation

    mu_chi(s) ~= alpha * L(s, chi) + beta * L'(s, chi)

holds for Dirichlet characters of different conductors.

Characters tested:
  chi4 (conductor 4 = 2^2) — baseline
  chi8 (conductor 8 = 2^3)
  chi3 (conductor 3, odd)
  chi5 (conductor 5, odd)

The carry stopping-time series n(tau) = 2*(tau - tau0) + 1 always yields
odd integers; characters with odd prime conductors (3, 5) assign chi(n) = 0
whenever n is divisible by the conductor, making the Dirichlet series sparse.
This predicts poor fit for odd conductors and good fit for conductors 2^k.

References: Paper L §6.4.
"""

import sys
import os
import cmath
import math
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _shared import build_highk_bank, mu_chi_from_record, primitive_characters, L_hurwitz

K_USE = 999


def L_prime(s, q, chi_fn, delta=1e-6):
    return (L_hurwitz(s + delta, q, chi_fn) - L_hurwitz(s - delta, q, chi_fn)) / (2 * delta)


def find_zeros(q, chi_fn, t_lo=0.5, t_hi=30.0, dt=0.05):
    """Find zeros of L(1/2+it, chi) by sign changes of Re or Im."""
    zeros = []
    t = t_lo
    prev_L = L_hurwitz(complex(0.5, t), q, chi_fn)
    while t < t_hi:
        t += dt
        cur_L = L_hurwitz(complex(0.5, t), q, chi_fn)
        for comp in ['real', 'imag']:
            pv = getattr(prev_L, comp)
            cv = getattr(cur_L, comp)
            if pv * cv < 0 and abs(cur_L) < 0.1:
                ta, tb, pv2 = t - dt, t, pv
                for _ in range(60):
                    tm = (ta + tb) / 2
                    vm = getattr(L_hurwitz(complex(0.5, tm), q, chi_fn), comp)
                    if pv2 * vm < 0:
                        tb = tm
                    else:
                        ta = tm; pv2 = vm
                t_zero = (ta + tb) / 2
                if abs(L_hurwitz(complex(0.5, t_zero), q, chi_fn)) < 0.01:
                    if not any(abs(t_zero - tz) < 0.5 for tz in zeros):
                        zeros.append(t_zero)
                break
        prev_L = cur_L
    return sorted(zeros)


def fit_two_term(bank, chi_fn, q, ts_grid):
    """Global two-term fit mu ~= A*L + B*L' on the critical line. Returns A, B, err%."""
    mu_v = np.array([mu_chi_from_record(bank[K_USE], chi_fn, complex(0.5, t))
                     for t in ts_grid])
    L_v  = np.array([L_hurwitz(complex(0.5, t), q, chi_fn) for t in ts_grid])
    Lp_v = np.array([L_prime(complex(0.5, t), q, chi_fn) for t in ts_grid])

    A_mat = np.column_stack([
        np.concatenate([L_v.real, L_v.imag]),
        np.concatenate([Lp_v.real, Lp_v.imag]),
    ])
    y_vec = np.concatenate([mu_v.real, mu_v.imag])
    cc, _, _, _ = np.linalg.lstsq(A_mat, y_vec, rcond=None)
    A_fit, B_fit = cc
    mu_approx = A_fit * L_v + B_fit * Lp_v
    err = np.mean(np.abs(mu_v - mu_approx)) / np.mean(np.abs(mu_v)) * 100
    return A_fit, B_fit, err


def main():
    bank = build_highk_bank()
    chars = {c["name"]: c for c in primitive_characters()}

    print("=" * 70)
    print("L21: CONDUCTOR RESTRICTION OF THE CRITICAL-LINE FIT")
    print("=" * 70)

    ts_grid = np.arange(1.0, 28.0, 0.15)

    char_list = [
        ("chi4", 4,  "2^2 (even)"),
        ("chi8", 8,  "2^3 (even)"),
        ("chi3", 3,  "3   (odd)"),
        ("chi5", 5,  "5   (odd)"),
    ]

    print(f"\n  {'character':>8}  {'conductor':>12}  {'A (=c_carry)':>14}  {'B':>12}  {'fit error':>10}")
    print("  " + "-" * 62)

    results = {}
    for name, q, desc in char_list:
        if name not in chars:
            print(f"  {name:>8}  {desc:>12}  (not available)")
            continue
        chi_fn = chars[name]["fn"]
        try:
            A, B, err = fit_two_term(bank, chi_fn, q, ts_grid)
            results[name] = (A, B, err)
            print(f"  {name:>8}  {desc:>12}  {A:>+14.6f}  {B:>+12.6f}  {err:>9.1f}%")
        except Exception as e:
            print(f"  {name:>8}  {desc:>12}  ERROR: {e}")

    print("\n--- Zeros of L(1/2+it, chi) ---")
    for name, q, desc in char_list:
        if name not in chars:
            continue
        chi_fn = chars[name]["fn"]
        zeros = find_zeros(q, chi_fn, t_hi=30.0)
        print(f"  {name} ({desc}): {len(zeros)} zeros up to t=30")
        for tz in zeros[:5]:
            print(f"    t = {tz:.5f}")

    print("\n--- Summary ---")
    print("  Character  Conductor  Fit error  Interpretation")
    for name, q, desc in char_list:
        if name not in results:
            continue
        A, B, err = results[name]
        quality = "good (2^k conductor)" if q in (4, 8) else "poor (odd conductor)"
        print(f"  {name:>8}  {q:>9}  {err:>8.1f}%  {quality}")

    print("\n  Mechanism: n(tau) = 2*(tau - tau0) + 1 is always odd.")
    print("  chi with odd prime conductor q: chi(n) = 0 for n divisible by q,")
    print("  making the Dirichlet series sparse and the linear fit unreliable.")


if __name__ == "__main__":
    main()
