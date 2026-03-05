#!/usr/bin/env python3
"""
P1-SPRINT-L6-mini: A/B control for local-corrector mechanism
=============================================================

Question:
  Does the (p=3, k)-pattern come from raw frequencies, or only from the full
  analytic carry-Dirichlet channel?

A/B/C observables (same local model fit for all):
  A) full:        mu_full(s)      = sum_tau u_mix(tau) * chi(n_tau) * n_tau^{-s}
  B) no_n_weight: mu_no_n(s)      = sum_tau u_mix(tau) * chi(n_tau)
  C) no_channel:  mu_no_chi(s)    = sum_tau u_mix(tau) * n_tau^{-s}

Model family (fixed p=3, k in {0,1,2}):
  mu(s) ~= c * (1 - chi(3) * 3^{-s})^k * L(s,chi)^2
"""

import cmath
import math
import random
from typing import Callable, Dict, List

try:
    import mpmath as mp
except Exception:
    mp = None

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def n_of_tau(tau: int, tau0: int = 3) -> int:
    return 2 * (tau - tau0) + 1


def mu_full(u_mix: Dict[int, float], chi_fn: Callable[[int], int], s: complex) -> complex:
    z = 0.0 + 0.0j
    for tau in sorted(u_mix.keys()):
        if tau < 3:
            continue
        n = n_of_tau(tau)
        c = chi_fn(n)
        if c == 0:
            continue
        z += (u_mix[tau] * c) * cmath.exp(-s * math.log(n))
    return z


def mu_no_n_weight(u_mix: Dict[int, float], chi_fn: Callable[[int], int], s: complex) -> complex:
    _ = s
    z = 0.0 + 0.0j
    for tau in sorted(u_mix.keys()):
        if tau < 3:
            continue
        n = n_of_tau(tau)
        c = chi_fn(n)
        if c == 0:
            continue
        z += (u_mix[tau] * c)
    return z


def mu_no_channel(u_mix: Dict[int, float], chi_fn: Callable[[int], int], s: complex) -> complex:
    _ = chi_fn
    z = 0.0 + 0.0j
    for tau in sorted(u_mix.keys()):
        if tau < 3:
            continue
        n = n_of_tau(tau)
        z += u_mix[tau] * cmath.exp(-s * math.log(n))
    return z


def choose_best_k_from_obs(
    obs: Callable[[complex], complex],
    chi_fn: Callable[[int], int],
    q: int,
    fit_grid: List[complex],
    hold_grid: List[complex],
):
    cp3 = chi_fn(3)

    def L(s: complex) -> complex:
        return _shared.L_hurwitz(s, q, chi_fn)

    by_k = []
    for k in [0, 1, 2]:
        def C(s: complex, k=k) -> complex:
            if cp3 == 0:
                return 1.0 + 0.0j
            return (1.0 - cp3 * cmath.exp(-s * math.log(3))) ** k

        basis_fit = [C(s) * (L(s) * L(s)) for s in fit_grid]
        mu_fit = [obs(s) for s in fit_grid]
        c = _shared.fit_c(mu_fit, basis_fit)
        fit_mre = _shared.mean([_shared.rel(obs(s), c * C(s) * (L(s) * L(s))) for s in fit_grid])
        hold_mre = _shared.mean([_shared.rel(obs(s), c * C(s) * (L(s) * L(s))) for s in hold_grid])
        by_k.append({"k": k, "fit_mre": fit_mre, "hold_mre": hold_mre})

    by_k.sort(key=lambda r: r["fit_mre"])
    best = by_k[0]
    second = by_k[1]
    margin = second["fit_mre"] - best["fit_mre"]
    return {"best_k": best["k"], "best_fit": best["fit_mre"], "best_hold": best["hold_mre"], "margin": margin, "all": by_k}


def main():
    flush("=" * 78)
    flush("P1-SPRINT-L6-mini: A/B control for local-corrector mechanism")
    flush("=" * 78)
    if mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = _shared.build_highk_bank()
    chars = _shared.primitive_characters()

    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]
    hold_grid = [complex(s, t) for s in [1.15, 1.30, 1.55, 1.85, 2.20] for t in [0.5, 1.5, 2.5]]

    obs_map = {
        "full": mu_full,
        "no_n_weight": mu_no_n_weight,
        "no_channel": mu_no_channel,
    }

    expected = {"chi3": 0, "chi4": 1, "chi5": 2, "chi8": 2}
    results = {name: {} for name in obs_map.keys()}

    for ch in chars:
        name, q, chi_fn = ch["name"], ch["q"], ch["fn"]
        u_mix = _shared.u_mix_map(bank[999])
        taus = [tau for tau in sorted(u_mix.keys()) if tau >= 3]
        n_map = {tau: n_of_tau(tau) for tau in taus}
        n_vals = [n_map[tau] for tau in taus]

        # Strong controls:
        rng_n = random.Random(101 + q)
        n_perm_vals = n_vals[:]
        rng_n.shuffle(n_perm_vals)
        n_perm = {tau: n_perm_vals[i] for i, tau in enumerate(taus)}

        rng_c = random.Random(211 + q)
        c_vals = [chi_fn(n) for n in n_vals]
        c_perm_vals = c_vals[:]
        rng_c.shuffle(c_perm_vals)
        c_perm_map = {n: c_perm_vals[i] for i, n in enumerate(n_vals)}

        for obs_name, obs_fn in obs_map.items():
            def obs(s: complex, obs_fn=obs_fn):
                return obs_fn(u_mix, chi_fn, s)
            results[obs_name][name] = choose_best_k_from_obs(obs, chi_fn, q, fit_grid, hold_grid)

        def obs_scrambled_n(s: complex):
            z = 0.0 + 0.0j
            for tau in taus:
                n = n_perm[tau]
                c = chi_fn(n)
                if c == 0:
                    continue
                z += (u_mix[tau] * c) * cmath.exp(-s * math.log(n))
            return z

        def obs_shuffled_chi(s: complex):
            z = 0.0 + 0.0j
            for tau in taus:
                n = n_map[tau]
                c = c_perm_map[n]
                if c == 0:
                    continue
                z += (u_mix[tau] * c) * cmath.exp(-s * math.log(n))
            return z

        results.setdefault("scrambled_n", {})[name] = choose_best_k_from_obs(obs_scrambled_n, chi_fn, q, fit_grid, hold_grid)
        results.setdefault("shuffled_chi", {})[name] = choose_best_k_from_obs(obs_shuffled_chi, chi_fn, q, fit_grid, hold_grid)

    flush(f"{'observable':>12s}  {'chi':>6s}  {'best_k':>6s}  {'fit':>10s}  {'hold':>10s}  {'margin':>10s}")
    for obs_name in ["full", "no_n_weight", "no_channel", "scrambled_n", "shuffled_chi"]:
        for chi in ["chi3", "chi4", "chi5", "chi8"]:
            r = results[obs_name][chi]
            flush(
                f"{obs_name:>12s}  {chi:>6s}  {r['best_k']:6d}  "
                f"{r['best_fit']:10.6f}  {r['best_hold']:10.6f}  {r['margin']:10.6f}"
            )

    # Pattern recovery
    full_match = all(results["full"][chi]["best_k"] == expected[chi] for chi in expected)
    no_n_match_count = sum(1 for chi in expected if results["no_n_weight"][chi]["best_k"] == expected[chi])
    no_ch_match_count = sum(1 for chi in expected if results["no_channel"][chi]["best_k"] == expected[chi])
    scr_n_match_count = sum(1 for chi in expected if results["scrambled_n"][chi]["best_k"] == expected[chi])
    shf_c_match_count = sum(1 for chi in expected if results["shuffled_chi"][chi]["best_k"] == expected[chi])

    # Separation on chi4/5/8 hold MRE against strong controls
    sep_ok = True
    for chi in ["chi4", "chi5", "chi8"]:
        h_full = results["full"][chi]["best_hold"]
        h_scr_n = results["scrambled_n"][chi]["best_hold"]
        h_shf_c = results["shuffled_chi"][chi]["best_hold"]
        if not (h_full <= min(h_scr_n, h_shf_c) * 0.75):
            sep_ok = False
            break

    g1 = full_match
    g2 = scr_n_match_count <= 1
    g3 = shf_c_match_count <= 1
    g4 = sep_ok

    flush("\n" + "=" * 78)
    flush("SPRINT-L6-mini GATES")
    flush("=" * 78)
    flush(f"L6m.1: full channel recovers expected k-pattern        -> {'PASS' if g1 else 'FAIL'}")
    flush(f"L6m.2: scrambled_n does NOT recover pattern (<=1/4)    -> {scr_n_match_count}/4  [{'PASS' if g2 else 'FAIL'}]")
    flush(f"L6m.3: shuffled_chi does NOT recover pattern (<=1/4)   -> {shf_c_match_count}/4  [{'PASS' if g3 else 'FAIL'}]")
    flush(f"L6m.4: full hold-MRE beats strong controls (chi4/5/8)  -> {'PASS' if g4 else 'FAIL'}")
    flush(f"diag: no_n_weight pattern matches = {no_n_match_count}/4, no_channel = {no_ch_match_count}/4")

    verdict = "PASS" if (g1 and g2 and g3 and g4) else "MIXED/INCONCLUSIVE"
    flush("\n" + "=" * 78)
    flush("SPRINT-L6-mini VERDICT")
    flush("=" * 78)
    flush(f"Result: {verdict}")
    if verdict == "PASS":
        flush("k-pattern is specific to the full analytic carry-Dirichlet channel.")
    else:
        flush("A/B controls are only partially discriminative with current setup.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

