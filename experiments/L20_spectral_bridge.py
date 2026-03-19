#!/usr/bin/env python3
"""
L20: Critical-line spectral bridge
====================================

Tests the two-term decomposition

    mu_chi4(s) ~= c_carry * L(s, chi4) + B(s) * L'(s, chi4)

on the critical line Re(s) = 1/2.

Measurements:
  1. c_carry = lim_{sigma->inf} c_tilde(sigma), where
     c_tilde(sigma) = mu_chi4(sigma) / ((1 + 3^{-sigma}) * L(sigma, chi4)^2)
  2. B(rho) = mu(rho) / L'(rho) at each zero rho of L(s, chi4)
  3. Power-law fit B(t) ~ a * t^{-p} and affine fit B(t) ~ B_inf + C/t
  4. Global fit error of the two-term formula on Re(s) = 1/2

References: Paper L §7.4 (spectral bridge), §7.3 problem 9.
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
PI = math.pi


def L_prime(s, q=4, chi_fn=None, delta=1e-6):
    return (L_hurwitz(s + delta, q, chi_fn) - L_hurwitz(s - delta, q, chi_fn)) / (2 * delta)


def find_zeros_chi4(t_lo=1.0, t_hi=60.0, dt=0.02):
    zeros = []
    from _shared import L_hurwitz as Lh
    chi4_fn = {c["name"]: c for c in primitive_characters()}["chi4"]["fn"]
    t = t_lo
    prev_L = Lh(complex(0.5, t), 4, chi4_fn)
    while t < t_hi:
        t += dt
        cur_L = Lh(complex(0.5, t), 4, chi4_fn)
        for comp in ['real', 'imag']:
            pv = getattr(prev_L, comp)
            cv = getattr(cur_L, comp)
            if pv * cv < 0 and abs(cur_L) < 0.5:
                ta, tb = t - dt, t
                pv2 = pv
                for _ in range(80):
                    tm = (ta + tb) / 2
                    Lm = Lh(complex(0.5, tm), 4, chi4_fn)
                    vm = getattr(Lm, comp)
                    if pv2 * vm < 0:
                        tb = tm
                    else:
                        ta = tm
                        pv2 = vm
                t_zero = (ta + tb) / 2
                L_check = Lh(complex(0.5, t_zero), 4, chi4_fn)
                if abs(L_check) < 1e-8:
                    if not any(abs(t_zero - tz) < 0.3 for tz in zeros):
                        zeros.append(t_zero)
                break
        prev_L = cur_L
    return sorted(zeros)


def main():
    bank = build_highk_bank()
    chars = {c["name"]: c for c in primitive_characters()}
    chi4_fn = chars["chi4"]["fn"]

    print("=" * 70)
    print("L20: CRITICAL-LINE SPECTRAL BRIDGE")
    print("=" * 70)

    # ── Part 1: c_carry = lim c_tilde(sigma) ──────────────────────────────
    print("\n--- Part 1: c_carry = lim_{sigma->inf} c_tilde(sigma) ---")
    sigmas = [3, 4, 5, 6, 8, 10, 15, 20]
    c_tilde_vals = {}
    print(f"\n  {'sigma':>6}  {'c_tilde':>14}  {'diff from 1/18':>16}  {'ratio*18':>10}")
    for sigma in sigmas:
        s = complex(float(sigma), 0.0)
        mu_s = mu_chi_from_record(bank[K_USE], chi4_fn, s)
        L_s = L_hurwitz(s, 4, chi4_fn)
        f_s = 1 + 3.0 ** (-sigma)
        c_t = (mu_s / (f_s * L_s ** 2)).real
        c_tilde_vals[sigma] = c_t
        print(f"  {sigma:>6d}  {c_t:>14.10f}  {c_t - 1/18:>+16.4e}  {c_t * 18:>10.8f}")

    c_carry = c_tilde_vals[20]
    print(f"\n  c_carry = c_tilde(sigma=20) = {c_carry:.10f}")
    print(f"  1/18                          = {1/18:.10f}")
    print(f"  gap c_carry - 1/18            = {c_carry - 1/18:+.4e}  ({abs(c_carry - 1/18)/(1/18)*100:.4f}%)")

    # ── Part 2: B(rho) at zeros of L ──────────────────────────────────────
    print("\n--- Part 2: B(rho) = mu(rho) / L'(rho) at zeros of L(s, chi4) ---")
    zeros = find_zeros_chi4(1.0, 60.0)
    print(f"\n  Found {len(zeros)} zeros up to t=60")
    B_at_zeros = []
    print(f"\n  {'#':>3}  {'t':>10}  {'B(rho)':>14}  {'|B|':>10}  {'B*t':>8}  {'B*t^0.23':>12}")
    for i, tz in enumerate(zeros):
        rho = complex(0.5, tz)
        mu_rho = mu_chi_from_record(bank[K_USE], chi4_fn, rho)
        Lp_rho = L_prime(rho, chi_fn=chi4_fn)
        B_rho = (mu_rho / Lp_rho).real
        B_at_zeros.append((tz, B_rho))
        print(f"  {i+1:>3d}  {tz:>10.4f}  {B_rho:>+14.8f}  {abs(B_rho):>10.6f}"
              f"  {B_rho * tz:>8.4f}  {B_rho * tz**0.23:>12.8f}")

    # ── Part 3: Functional form of B(t) ───────────────────────────────────
    print("\n--- Part 3: Functional form of B(t) ---")
    t_vals = np.array([tz for tz, _ in B_at_zeros])
    B_vals = np.array([B for _, B in B_at_zeros])

    if len(t_vals) >= 4:
        # Power law: B(t) = a / t^p
        log_t = np.log(t_vals)
        log_B = np.log(np.abs(B_vals))
        A_mat = np.column_stack([np.ones_like(log_t), log_t])
        coeffs = np.linalg.lstsq(A_mat, log_B, rcond=None)[0]
        a_pw = math.exp(coeffs[0])
        p_pw = -coeffs[1]
        err_pw = np.mean(np.abs(B_vals - a_pw * t_vals**(-p_pw))) / np.mean(np.abs(B_vals))
        print(f"\n  Power law:  B(t) = {a_pw:.6f} * t^{{-{p_pw:.4f}}}")
        print(f"    Relative error: {err_pw*100:.2f}%")

        # Affine: B(t) = B_inf + C/t
        A_inv = np.column_stack([np.ones_like(t_vals), 1.0 / t_vals])
        cc_inv = np.linalg.lstsq(A_inv, B_vals, rcond=None)[0]
        err_inv = np.mean(np.abs(B_vals - (cc_inv[0] + cc_inv[1] / t_vals))) / np.mean(np.abs(B_vals))
        print(f"\n  Affine:     B(t) = {cc_inv[0]:.6f} + {cc_inv[1]:.4f}/t")
        print(f"    B(inf) = {cc_inv[0]:.6f}")
        print(f"    Relative error: {err_inv*100:.2f}%")

    # ── Part 4: Global fit error on critical line ──────────────────────────
    print("\n--- Part 4: Two-term fit error on Re(s) = 1/2 ---")
    ts_grid = np.arange(1.0, 55.0, 0.1)
    mu_vals = np.array([mu_chi_from_record(bank[K_USE], chi4_fn, complex(0.5, t))
                        for t in ts_grid])
    L_vals = np.array([L_hurwitz(complex(0.5, t), 4, chi4_fn) for t in ts_grid])
    Lp_vals = np.array([L_prime(complex(0.5, t), chi_fn=chi4_fn) for t in ts_grid])

    # Fit: mu ~ A*L + B_coeff*Lp  (global constant coefficients)
    A_mat2 = np.column_stack([
        np.concatenate([L_vals.real, L_vals.imag]),
        np.concatenate([Lp_vals.real, Lp_vals.imag]),
    ])
    y_vec = np.concatenate([mu_vals.real, mu_vals.imag])
    cc2, _, _, _ = np.linalg.lstsq(A_mat2, y_vec, rcond=None)
    A_fit, B_fit = cc2
    mu_approx = A_fit * L_vals + B_fit * Lp_vals
    err_rel = np.mean(np.abs(mu_vals - mu_approx)) / np.mean(np.abs(mu_vals))

    print(f"\n  Global fit:  mu ~= A*L + B*L'")
    print(f"    A = {A_fit:.8f}  (= c_carry, expected ~{c_carry:.6f})")
    print(f"    B = {B_fit:.8f}")
    print(f"    Mean relative error: {err_rel*100:.2f}%")

    print("\n--- Summary ---")
    print(f"  c_carry = {c_carry:.8f}  (1/18 = {1/18:.8f}, gap = {c_carry - 1/18:+.2e})")
    if len(B_at_zeros) >= 4:
        print(f"  B(t) ~ {a_pw:.4f} * t^(-{p_pw:.4f})  (power law, error {err_pw*100:.1f}%)")
        print(f"  B(inf) = {cc_inv[0]:.4f}  (affine extrapolation)")
    print(f"  Two-term fit error on critical line: {err_rel*100:.2f}%")


if __name__ == "__main__":
    main()
