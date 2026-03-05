#!/usr/bin/env python3
"""
P1-STEP3B: Carry-side vs Witt-side operator intertwiner (base 2)
=================================================================

Step-3B objective:
  Build an explicit operator-level bridge between
  - carry-side dynamics, and
  - Witt-correction-side dynamics,
  and quantify mismatch.

Model:
  c_{j+1} = floor((x_j + y_j + c_j)/2), with x_j,y_j iid Bernoulli(1/2).
  Augmented carry state: z_j = (c_j, c_{j+1}) in {0,1}^2.
  Witt local correction state:
      delta_j = c_j - 2*c_{j+1} in {0, -2, +1, -1}.

The map phi : z -> delta is bijective, hence defines a permutation intertwiner.
"""

import argparse
import math
import numpy as np


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


Z_STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]
D_STATES = [0, -2, +1, -1]  # image of phi in Z_STATES order
Z_INDEX = {s: i for i, s in enumerate(Z_STATES)}
D_INDEX = {d: i for i, d in enumerate(D_STATES)}


def phi(z):
    c0, c1 = z
    return c0 - 2 * c1


def carry_next_prob(c_in):
    # c_out = 1 probability conditioned on carry-in.
    # c_in=0: P(sum>=2)=1/4. c_in=1: P(sum+1>=2)=3/4.
    return 0.25 if c_in == 0 else 0.75


def build_exact_z_operator():
    """
    Build exact transition on z_j=(c_j,c_{j+1}) -> z_{j+1}=(c_{j+1},c_{j+2}).
    Column-stochastic convention.
    """
    T = np.zeros((4, 4), dtype=float)
    for z in Z_STATES:
        src = Z_INDEX[z]
        c_j, c_j1 = z
        p1 = carry_next_prob(c_j1)
        for c_j2, prob in [(0, 1.0 - p1), (1, p1)]:
            z2 = (c_j1, c_j2)
            dst = Z_INDEX[z2]
            T[dst, src] += prob
    return T


def permutation_z_to_d():
    """
    Build permutation matrix P such that v_d = P v_z.
    Then T_d = P T_z P^{-1}.
    """
    P = np.zeros((4, 4), dtype=float)
    for z in Z_STATES:
        i_z = Z_INDEX[z]
        i_d = D_INDEX[phi(z)]
        P[i_d, i_z] = 1.0
    return P


def carries_for_pair(a, b, n_bits):
    c = [0] * (n_bits + 1)
    for j in range(n_bits):
        aj = (a >> j) & 1
        bj = (b >> j) & 1
        c[j + 1] = (aj + bj + c[j]) >> 1
    return c


def build_empirical_operators(n_bits):
    """
    Enumerate all pairs (a,b) mod 2^n and estimate:
      - Tz_emp on z_j=(c_j,c_{j+1}),
      - Td_emp on delta_j=c_j-2c_{j+1}.
    """
    mod = 1 << n_bits
    C_z = np.zeros((4, 4), dtype=float)
    C_d = np.zeros((4, 4), dtype=float)
    out_z = np.zeros((4,), dtype=float)
    out_d = np.zeros((4,), dtype=float)

    for a in range(mod):
        for b in range(mod):
            c = carries_for_pair(a, b, n_bits)
            # transitions for j=0..n_bits-2
            for j in range(n_bits - 1):
                z1 = (c[j], c[j + 1])
                z2 = (c[j + 1], c[j + 2])
                i1 = Z_INDEX[z1]
                i2 = Z_INDEX[z2]
                C_z[i2, i1] += 1.0
                out_z[i1] += 1.0

                d1 = phi(z1)
                d2 = phi(z2)
                k1 = D_INDEX[d1]
                k2 = D_INDEX[d2]
                C_d[k2, k1] += 1.0
                out_d[k1] += 1.0

    Tz = np.zeros((4, 4), dtype=float)
    Td = np.zeros((4, 4), dtype=float)
    for i in range(4):
        if out_z[i] > 0:
            Tz[:, i] = C_z[:, i] / out_z[i]
        if out_d[i] > 0:
            Td[:, i] = C_d[:, i] / out_d[i]
    return Tz, Td


def fro_norm(x):
    return float(np.linalg.norm(x, ord="fro"))


def spectral_info(T):
    eig = np.linalg.eigvals(T)
    eig_sorted = sorted(eig, key=lambda z: -abs(z))
    return eig_sorted


def fmt_eigs(eigs, k=4):
    out = []
    for e in eigs[:k]:
        out.append(f"{e.real:+.8f}{e.imag:+.8f}i |{abs(e):.8f}|")
    return out


def main():
    parser = argparse.ArgumentParser(description="Step3B operator intertwiner checks.")
    parser.add_argument("--n-bits", type=int, default=10, help="Enumeration modulus: 2^n")
    args = parser.parse_args()

    flush("=" * 78)
    flush("P1-STEP3B: carry-side vs Witt-side operator intertwiner")
    flush("=" * 78)
    flush(f"Enumeration level: n_bits={args.n_bits} (pairs={1 << (2*args.n_bits):,})")

    # Exact operators
    Tz_exact = build_exact_z_operator()
    P = permutation_z_to_d()
    Pinv = np.linalg.inv(P)
    Td_exact = P @ Tz_exact @ Pinv

    flush("\n" + "=" * 78)
    flush("A) Exact operator bridge")
    flush("=" * 78)
    flush("Tz_exact (carry augmented state z=(c_j,c_{j+1})):")
    flush(str(Tz_exact))
    flush("\nTd_exact (Witt correction state delta):")
    flush(str(Td_exact))

    intertwiner_exact = fro_norm(Td_exact @ P - P @ Tz_exact)
    flush(f"\nExact intertwiner residual ||Td P - P Tz||_F = {intertwiner_exact:.3e}")

    eig_z = spectral_info(Tz_exact)
    eig_d = spectral_info(Td_exact)
    flush("\nExact spectra:")
    flush("  eig(Tz): " + " ; ".join(fmt_eigs(eig_z)))
    flush("  eig(Td): " + " ; ".join(fmt_eigs(eig_d)))

    # Empirical operators
    flush("\n" + "=" * 78)
    flush("B) Empirical operator bridge (full enumeration)")
    flush("=" * 78)
    Tz_emp, Td_emp = build_empirical_operators(args.n_bits)

    residual_z = fro_norm(Tz_emp - Tz_exact)
    residual_d = fro_norm(Td_emp - Td_exact)
    intertwiner_emp = fro_norm(Td_emp @ P - P @ Tz_emp)
    conjugacy_emp = fro_norm(Td_emp - (P @ Tz_emp @ Pinv))

    flush(f"||Tz_emp - Tz_exact||_F = {residual_z:.3e}")
    flush(f"||Td_emp - Td_exact||_F = {residual_d:.3e}")
    flush(f"Empirical intertwiner residual ||Td_emp P - P Tz_emp||_F = {intertwiner_emp:.3e}")
    flush(f"Empirical conjugacy residual ||Td_emp - P Tz_emp P^-1||_F = {conjugacy_emp:.3e}")

    eig_z_emp = spectral_info(Tz_emp)
    eig_d_emp = spectral_info(Td_emp)
    flush("\nEmpirical spectra:")
    flush("  eig(Tz_emp): " + " ; ".join(fmt_eigs(eig_z_emp)))
    flush("  eig(Td_emp): " + " ; ".join(fmt_eigs(eig_d_emp)))

    # Step-3B summary metrics
    flush("\n" + "=" * 78)
    flush("STEP3B VERDICT")
    flush("=" * 78)
    pass_exact = intertwiner_exact < 1e-12
    pass_emp = max(residual_z, residual_d, intertwiner_emp, conjugacy_emp) < 1e-12

    flush(f"Exact intertwiner: {'PASS' if pass_exact else 'FAIL'}")
    flush(f"Empirical bridge:  {'PASS' if pass_emp else 'FAIL'}")
    flush("Second eigenvalue magnitude (exact): "
          f"|lambda_2| = {abs(eig_z[1]):.8f} (expected 1/2)")

    if pass_exact and pass_emp:
        flush("\nConclusion:")
        flush("  Carry-side and Witt-correction-side operators are exactly intertwined")
        flush("  (in this base-2 addition model) by a permutation map.")
        flush("  This gives a concrete operator-level Step-3B bridge with spectral match.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

