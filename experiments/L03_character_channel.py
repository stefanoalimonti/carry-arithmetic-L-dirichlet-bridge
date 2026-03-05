#!/usr/bin/env python3
"""
P1-STEP3C: Character-projected trace channel on the carry/Witt bridge
======================================================================

Step-3C goal:
  Define a chi_4-like character channel on the Step-3B intertwined operators
  and verify that carry-side and Witt-side resolvent channels coincide.

Channel definition (delta side):
  States delta in {0, -2, +1, -1} map to residues mod 4: {0,2,1,3}.
  Define chi_4-like weights:
    w(delta) = 0 for even residues (0,2),
    w(delta) = +1 for residue 1,
    w(delta) = -1 for residue 3.

  Hence in D-state order [0, -2, +1, -1]:
    w_d = [0, 0, +1, -1].

We track:
  c_n = w^T T^n e0
and its centered version:
  d_n = c_n - c_inf.

The centered channel isolates the transient spectral content (lambda_2 = 1/2).
"""

import argparse
import numpy as np

import L02_operator_intertwiner as L02


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def stationary_distribution_column_stochastic(T):
    """
    Return normalized right eigenvector for eigenvalue 1.
    T is column-stochastic, so stationary pi satisfies T pi = pi.
    """
    vals, vecs = np.linalg.eig(T)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    v = vecs[:, idx].real
    # enforce positivity orientation
    if v.sum() < 0:
        v = -v
    v = np.maximum(v, 0.0)
    s = v.sum()
    if s <= 1e-15:
        raise RuntimeError("Failed to extract stationary distribution.")
    return v / s


def channel_coeffs(T, w, e0, n_max):
    coeff = []
    v = e0.copy()
    for _ in range(n_max + 1):
        coeff.append(float(w @ v))
        v = T @ v
    return np.array(coeff, dtype=float)


def main():
    parser = argparse.ArgumentParser(description="Step3C character-projected trace channel.")
    parser.add_argument("--n-max", type=int, default=20, help="Max power n for channel coefficients.")
    parser.add_argument("--z-grid", type=str, default="0.2,0.4,0.6,0.8",
                        help="Comma-separated z values for resolvent checks.")
    args = parser.parse_args()

    Tz = L02.build_exact_z_operator()
    P = L02.permutation_z_to_d()
    Pinv = np.linalg.inv(P)
    Td = P @ Tz @ Pinv

    # chi_4-like channel on delta states
    w_d = np.array([0.0, 0.0, +1.0, -1.0], dtype=float)
    # map to z-side channel
    w_z = w_d @ P

    # initial state: no incoming carries at top (z=(0,0), delta=0)
    e0_z = np.zeros(4, dtype=float)
    e0_z[L02.Z_INDEX[(0, 0)]] = 1.0
    e0_d = P @ e0_z

    flush("=" * 78)
    flush("P1-STEP3C: character-projected trace channel")
    flush("=" * 78)
    flush(f"n_max={args.n_max}")
    flush(f"w_d={w_d.tolist()}  (chi_4-like on residues mod 4)")
    flush(f"w_z={w_z.tolist()}  (mapped via intertwiner)")

    # A) coefficient-level equality
    c_d = channel_coeffs(Td, w_d, e0_d, args.n_max)
    c_z = channel_coeffs(Tz, w_z, e0_z, args.n_max)
    coeff_err = float(np.max(np.abs(c_d - c_z)))

    flush("\n" + "=" * 78)
    flush("A) Coefficient channel equality (carry vs Witt)")
    flush("=" * 78)
    flush(f"max_n |c_d(n)-c_z(n)| = {coeff_err:.3e}")
    flush(f"{'n':>3s}  {'c_d(n)':>14s}  {'c_z(n)':>14s}  {'diff':>12s}")
    for n in range(min(args.n_max + 1, 12)):
        flush(f"{n:3d}  {c_d[n]:+14.9f}  {c_z[n]:+14.9f}  {c_d[n]-c_z[n]:+12.3e}")

    # B) stationary + centered decomposition
    pi_d = stationary_distribution_column_stochastic(Td)
    c_inf = float(w_d @ pi_d)
    d = c_d - c_inf

    flush("\n" + "=" * 78)
    flush("B) Stationary split and transient rate")
    flush("=" * 78)
    flush(f"stationary pi_d = {pi_d}")
    flush(f"c_inf = w_d^T pi_d = {c_inf:+.9f}")

    flush(f"\n{'n':>3s}  {'c_n':>14s}  {'d_n=c_n-c_inf':>16s}  {'ratio d_{n+1}/d_n':>18s}")
    for n in range(min(args.n_max, 12)):
        ratio = d[n + 1] / d[n] if abs(d[n]) > 1e-15 else float("nan")
        flush(f"{n:3d}  {c_d[n]:+14.9f}  {d[n]:+16.9f}  {ratio:+18.9f}")

    # expected transient mode 1/2
    ratios = []
    for n in range(2, min(args.n_max, 16)):
        if abs(d[n]) > 1e-14:
            ratios.append(d[n + 1] / d[n])
    ratio_tail = float(np.mean(ratios[-5:])) if len(ratios) >= 5 else float("nan")
    flush(f"\nMean tail ratio d_(n+1)/d_n: {ratio_tail:.9f}  (expected 0.5)")

    # C) resolvent channel at selected z
    z_values = [float(x.strip()) for x in args.z_grid.split(",") if x.strip()]
    flush("\n" + "=" * 78)
    flush("C) Resolvent channel equality and closed form (centered)")
    flush("=" * 78)
    flush(f"{'z':>6s}  {'C_d(z)':>14s}  {'C_z(z)':>14s}  {'diff':>12s}  {'Centered model':>16s}")
    for z in z_values:
        I = np.eye(4)
        Cd = float(w_d @ np.linalg.solve(I - z * Td, e0_d))
        Cz = float(w_z @ np.linalg.solve(I - z * Tz, e0_z))
        # Centered generating function: G(z)=sum d_n z^n = (c(z) - c_inf/(1-z))
        Gd = Cd - c_inf / (1.0 - z)
        # Exact centered sequence here has a short boundary transient:
        # d_0 = 1/4, and for n>=1: d_n = 2^{-n}.
        # Therefore:
        #   G(z)=sum_{n>=0} d_n z^n = 1/4 + sum_{n>=1} 2^{-n} z^n
        #       = 1/4 + z/(4-2z) = (2+z)/(4*(2-z)).
        G_model = (2.0 + z) / (4.0 * (2.0 - z))
        flush(f"{z:6.3f}  {Cd:+14.9f}  {Cz:+14.9f}  {Cd-Cz:+12.3e}  {G_model:+16.9f}")
        flush(f"       centered G_d(z)={Gd:+.9f}  model err={Gd-G_model:+.3e}")

    # D) verdict
    flush("\n" + "=" * 78)
    flush("STEP3C VERDICT")
    flush("=" * 78)
    pass_coeff = coeff_err < 1e-12
    pass_rate = abs(ratio_tail - 0.5) < 1e-9 if ratio_tail == ratio_tail else False

    flush(f"Carry/Witt channel equality: {'PASS' if pass_coeff else 'FAIL'}")
    flush(f"Centered transient rate (1/2): {'PASS' if pass_rate else 'FAIL'}")
    flush("Interpretation:")
    flush("  - Character-projected resolvent channel is identical on the two sides.")
    flush("  - After removing stationary contribution, the channel is purely 1/2-mode.")
    flush("  - This is an operator-level character trace channel in the Step-3 bridge.")
    flush("  - It does NOT yet prove the full D-odd multiplication chi_4->L(1,chi_4) step.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

