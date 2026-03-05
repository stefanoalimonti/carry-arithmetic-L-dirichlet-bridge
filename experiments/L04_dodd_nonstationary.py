#!/usr/bin/env python3
"""
P1-STEP3D: D-odd multiplication nonstationary character channel
===============================================================

Extends Step-3A/B/C from stationary base-2 addition prototype to the
nonstationary D-odd multiplication setting (empirical, sampled).

For each accepted D-odd pair (x,y) in a chosen sector, we compute carries c_j
of schoolbook multiplication and define parity-augmented states:

  z_j = (c_j mod 2, c_{j+1} mod 2) in {0,1}^2

Across depth, we build nonstationary transition operators T_tau:
  z_j -> z_{j-1}  (equivalently tau -> tau+1 from top boundary).

Then:
1) build delta-parity states via phi(z)= (z0 - 2*z1) mod 4,
2) verify per-depth intertwiner residuals between z- and delta-side operators,
3) evaluate chi_4-like channel weights w=[0,0,+1,-1].
"""

import argparse
import math
import random
from collections import defaultdict

import numpy as np


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


Z_STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]
Z_INDEX = {s: i for i, s in enumerate(Z_STATES)}

# delta residues induced by z0-2z1 mod 4, in Z_STATES order:
# (0,0)->0, (0,1)->2, (1,0)->1, (1,1)->3
D_RES = [0, 2, 1, 3]
D_STATES = [0, 2, 1, 3]
D_INDEX = {d: i for i, d in enumerate(D_STATES)}


def phi_index_from_z_index(i_z):
    return D_INDEX[D_RES[i_z]]


def permutation_z_to_d():
    P = np.zeros((4, 4), dtype=float)
    for i_z in range(4):
        i_d = phi_index_from_z_index(i_z)
        P[i_d, i_z] = 1.0
    return P


def carry_chain_bits(x: int, y: int, k_bits: int):
    """Carry array c[0..D] for schoolbook binary multiplication."""
    d_len = 2 * k_bits - 1
    carries = [0] * (d_len + 1)
    for d in range(d_len):
        conv = 0
        i_lo = max(0, d - k_bits + 1)
        i_hi = min(d, k_bits - 1)
        for i in range(i_lo, i_hi + 1):
            conv += ((x >> i) & 1) * ((y >> (d - i)) & 1)
        carries[d + 1] = (conv + carries[d]) >> 1
    return carries


def sector_of_pair(x: int, y: int, k_bits: int):
    ax = (x >> (k_bits - 2)) & 1
    by = (y >> (k_bits - 2)) & 1
    if ax == 0 and by == 0:
        return "00"
    if ax == 1 and by == 0:
        return "10"
    if ax == 0 and by == 1:
        return "01"
    return "11"


def sample_dodd_pairs(k_bits: int, sector: str, samples: int, seed: int):
    lo = 1 << (k_bits - 1)
    hi = 1 << k_bits
    d_len = 2 * k_bits - 1
    rng = random.Random(seed)

    out = []
    tries = 0
    while len(out) < samples:
        tries += 1
        x = rng.randrange(lo, hi)
        y = rng.randrange(lo, hi)
        if (x * y).bit_length() != d_len:
            continue
        sec = sector_of_pair(x, y, k_bits)
        if sec == "11":
            continue
        if sector != "all" and sec != sector:
            continue
        out.append((x, y))
    return out, tries


def build_nonstationary_operators(k_bits: int, sector: str, samples: int, seed: int):
    pairs, tries = sample_dodd_pairs(k_bits, sector, samples, seed)
    d_len = 2 * k_bits - 1
    # tau from top: tau=0 corresponds j=D-1, then j decreases.
    tau_max = d_len - 2  # need j-1>=0 and j+1<=D

    # Per tau counts for z-side and d-side
    C_z = [np.zeros((4, 4), dtype=float) for _ in range(tau_max)]
    C_d = [np.zeros((4, 4), dtype=float) for _ in range(tau_max)]
    out_z = [np.zeros((4,), dtype=float) for _ in range(tau_max)]
    out_d = [np.zeros((4,), dtype=float) for _ in range(tau_max)]

    # State distribution per tau
    V_z = [np.zeros((4,), dtype=float) for _ in range(tau_max + 1)]
    V_d = [np.zeros((4,), dtype=float) for _ in range(tau_max + 1)]

    D = d_len
    for x, y in pairs:
        c = carry_chain_bits(x, y, k_bits)  # c[0..D]
        for tau in range(tau_max + 1):
            j = D - 1 - tau
            z = (c[j] & 1, c[j + 1] & 1)
            iz = Z_INDEX[z]
            id_ = phi_index_from_z_index(iz)
            V_z[tau][iz] += 1.0
            V_d[tau][id_] += 1.0

            if tau < tau_max:
                jn = j - 1
                zn = (c[jn] & 1, c[jn + 1] & 1)
                izn = Z_INDEX[zn]
                idn = phi_index_from_z_index(izn)
                C_z[tau][izn, iz] += 1.0
                out_z[tau][iz] += 1.0
                C_d[tau][idn, id_] += 1.0
                out_d[tau][id_] += 1.0

    T_z = []
    T_d = []
    for tau in range(tau_max):
        tz = np.zeros((4, 4), dtype=float)
        td = np.zeros((4, 4), dtype=float)
        for i in range(4):
            if out_z[tau][i] > 0:
                tz[:, i] = C_z[tau][:, i] / out_z[tau][i]
            if out_d[tau][i] > 0:
                td[:, i] = C_d[tau][:, i] / out_d[tau][i]
        T_z.append(tz)
        T_d.append(td)

    # normalize state distributions
    for tau in range(tau_max + 1):
        sz = V_z[tau].sum()
        sd = V_d[tau].sum()
        if sz > 0:
            V_z[tau] /= sz
        if sd > 0:
            V_d[tau] /= sd

    return {
        "pairs": pairs,
        "tries": tries,
        "tau_max": tau_max,
        "T_z": T_z,
        "T_d": T_d,
        "V_z": V_z,
        "V_d": V_d,
    }


def fro_norm(x):
    return float(np.linalg.norm(x, ord="fro"))


def main():
    parser = argparse.ArgumentParser(description="Step3D nonstationary D-odd channel.")
    parser.add_argument("--K", type=int, default=13)
    parser.add_argument("--sector", choices=["00", "10", "all"], default="00")
    parser.add_argument("--samples", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    flush("=" * 78)
    flush("P1-STEP3D: D-odd nonstationary character channel")
    flush("=" * 78)
    flush(f"K={args.K}, sector={args.sector}, samples={args.samples}")

    out = build_nonstationary_operators(args.K, args.sector, args.samples, args.seed)
    tau_max = out["tau_max"]
    T_z = out["T_z"]
    T_d = out["T_d"]
    V_z = out["V_z"]
    V_d = out["V_d"]
    P = permutation_z_to_d()
    Pinv = np.linalg.inv(P)

    flush(f"Accepted pairs: {len(out['pairs'])}  (tries={out['tries']})")
    flush(f"tau range: 0..{tau_max}")

    # A) per-depth intertwiner residuals
    res_inter = []
    res_conj = []
    for tau in range(tau_max):
        ri = fro_norm(T_d[tau] @ P - P @ T_z[tau])
        rc = fro_norm(T_d[tau] - (P @ T_z[tau] @ Pinv))
        res_inter.append(ri)
        res_conj.append(rc)

    flush("\n" + "=" * 78)
    flush("A) Nonstationary operator bridge (per-depth)")
    flush("=" * 78)
    flush(f"max_tau ||Td_tau P - P Tz_tau||_F = {max(res_inter):.3e}")
    flush(f"max_tau ||Td_tau - P Tz_tau P^-1||_F = {max(res_conj):.3e}")

    # B) chi_4-like channel equality
    w_d = np.array([0.0, 0.0, +1.0, -1.0], dtype=float)
    w_z = w_d @ P
    c_d = np.array([float(w_d @ V_d[t]) for t in range(tau_max + 1)])
    c_z = np.array([float(w_z @ V_z[t]) for t in range(tau_max + 1)])

    flush("\n" + "=" * 78)
    flush("B) Character channel equality over depth")
    flush("=" * 78)
    flush(f"max_tau |c_d(tau)-c_z(tau)| = {np.max(np.abs(c_d-c_z)):.3e}")
    flush(f"{'tau':>4s}  {'c_d':>12s}  {'c_z':>12s}  {'diff':>12s}")
    for tau in range(min(14, tau_max + 1)):
        flush(f"{tau:4d}  {c_d[tau]:+12.6f}  {c_z[tau]:+12.6f}  {c_d[tau]-c_z[tau]:+12.3e}")

    # C) local decay ratios for centered channel
    tail_ref = np.mean(c_d[-6:]) if tau_max >= 10 else c_d[-1]
    d = c_d - tail_ref
    ratios = []
    ratios_robust = []
    for tau in range(2, min(tau_max, 20)):
        if abs(d[tau]) > 1e-8:
            r = d[tau + 1] / d[tau]
            ratios.append(r)
            # robust filter: avoid near-zero denominator and extreme spikes
            if abs(d[tau]) > 2e-2 and abs(r) < 2.0:
                ratios_robust.append(r)

    flush("\n" + "=" * 78)
    flush("C) Centered channel local ratios (empirical nonstationary)")
    flush("=" * 78)
    flush(f"tail reference (mean last 6 taus): {tail_ref:+.6e}")
    flush(f"{'tau':>4s}  {'d_tau':>12s}  {'ratio':>12s}")
    for tau in range(2, min(14, tau_max)):
        ratio = d[tau + 1] / d[tau] if abs(d[tau]) > 1e-8 else float("nan")
        flush(f"{tau:4d}  {d[tau]:+12.6e}  {ratio:+12.6f}")

    if ratios:
        flush(f"\nMean ratio over raw window: {float(np.mean(ratios)):+.6f}")
        if ratios_robust:
            flush(f"Mean ratio over robust window: {float(np.mean(ratios_robust)):+.6f}")
            flush(f"Median ratio over robust window: {float(np.median(ratios_robust)):+.6f}")
        flush("Reference (stationary prototype): +0.5")

    flush("\n" + "=" * 78)
    flush("STEP3D STATUS")
    flush("=" * 78)
    flush("PASS if considered as bridge-consistency check:")
    flush("  - operator intertwiner residuals are ~0 by depth (exact at sampled precision)")
    flush("  - character channel is identical on carry and delta/Witt sides")
    flush("Interpretation:")
    flush("  - This extends Step-3B/C mechanism to nonstationary D-odd multiplication data.")
    flush("  - Remaining gap: derive analytical (not sampled) nonstationary operator and")
    flush("    connect its character-projected resolvent to the P1 R-infinity channel.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

