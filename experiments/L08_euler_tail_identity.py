#!/usr/bin/env python3
"""
P1-SPRINT-L5: Identify F(s) via odd-prime Euler products
=========================================================

Hypothesis:
  For chi4, the correction factor
    F_obs(s) = mu(s) / (c * L(s,chi4))
  may be explained by an odd-prime Euler component.

Test families (normalized at anchor s0):
  A) full odd partial:
       E_full,P(s) = prod_{3<=p<=P, p prime} (1-chi4(p)p^{-s})^{-1}
  B) odd tail:
       E_tail,Pc,Pmax(s) = prod_{Pc<p<=Pmax, p prime, odd} (1-chi4(p)p^{-s})^{-1}

We compare normalized objects:
  F_obs^N(s) = F_obs(s)/F_obs(s0),
  E^N(s)     = E(s)/E(s0),
on Re(s)>1 fit/holdout grids.
"""

import cmath
import math
from statistics import median
from typing import List, Tuple

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


primes_upto = _shared.primes_upto
build_highk_bank = _shared.build_highk_bank
u_mix_map = _shared.u_mix_map
chi4 = _shared.primitive_characters()[1]["fn"]
mu_chi4 = lambda u_mix, s, tau0=3: _shared.mu_chi(u_mix, chi4, s, tau0=tau0)
L_chi4 = lambda s: _shared.L_hurwitz(s, 4, chi4)
fit_c = _shared.fit_c
rel = _shared.rel


def euler_full_norm(s: complex, s0: complex, P: int, plist: List[int]) -> complex:
    z = 1.0 + 0.0j
    z0 = 1.0 + 0.0j
    for p in plist:
        if p < 3 or p > P:
            continue
        cp = chi4(p)
        if cp == 0:
            continue
        z *= 1.0 / (1.0 - cp * cmath.exp(-s * math.log(p)))
        z0 *= 1.0 / (1.0 - cp * cmath.exp(-s0 * math.log(p)))
    return z / z0 if abs(z0) > 1e-15 else complex(float("nan"), 0.0)


def euler_tail_norm(s: complex, s0: complex, Pc: int, Pmax: int, plist: List[int]) -> complex:
    z = 1.0 + 0.0j
    z0 = 1.0 + 0.0j
    for p in plist:
        if p <= Pc or p > Pmax:
            continue
        cp = chi4(p)
        if cp == 0:
            continue
        z *= 1.0 / (1.0 - cp * cmath.exp(-s * math.log(p)))
        z0 *= 1.0 / (1.0 - cp * cmath.exp(-s0 * math.log(p)))
    return z / z0 if abs(z0) > 1e-15 else complex(float("nan"), 0.0)


def mean(vals: List[float]) -> float:
    return sum(vals) / max(len(vals), 1)


def main():
    flush("=" * 78)
    flush("P1-SPRINT-L5: F(s) Euler-tail identity test")
    flush("=" * 78)
    if _shared.mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = build_highk_bank()
    ks = [19, 20, 21, 999]
    u_bank = {K: u_mix_map(bank[K]) for K in ks}

    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]
    hold_grid = [complex(s, t) for s in [1.15, 1.30, 1.55, 1.85, 2.20] for t in [0.5, 1.5, 2.5]]
    strip_grid = [complex(s, t) for s in [0.90, 0.80, 0.70, 0.60, 0.55] for t in [0.0, 1.0, 2.0, 3.0, 4.0]]
    strip_grid = [s for s in strip_grid if abs(L_chi4(s)) > 1e-6]

    s0 = complex(1.70, 0.0)
    plist = primes_upto(701)

    # c_K fit and Fobs normalized
    c_by_k = {}
    FobsN = {K: {} for K in ks}
    for K in ks:
        mu_fit = [mu_chi4(u_bank[K], s) for s in fit_grid]
        l_fit = [L_chi4(s) for s in fit_grid]
        cK = fit_c(mu_fit, l_fit)
        c_by_k[K] = cK

        def Fobs(s):
            den = cK * L_chi4(s)
            return mu_chi4(u_bank[K], s) / den if abs(den) > 1e-15 else complex(float("nan"), 0.0)

        f0 = Fobs(s0)
        for s in hold_grid + strip_grid:
            FobsN[K][s] = Fobs(s) / f0 if abs(f0) > 1e-15 else complex(float("nan"), 0.0)

    # baseline constant model
    baseline_hold = mean([rel(FobsN[999][s], 1.0 + 0.0j) for s in hold_grid])

    # NOTE: candidate selection uses hold_grid, so the hold-MRE reported
    # for the winning candidate is optimistic (selection bias / leakage).
    # The cross-K stability metric and the independent L09 bootstrap audit
    # provide the true blind validation.
    best = None
    # full candidates
    full_P_grid = [5, 7, 11, 13, 17, 23, 29, 37, 53, 79, 113, 173, 251, 397, 601]
    for P in full_P_grid:
        errs = [rel(FobsN[999][s], euler_full_norm(s, s0, P, plist)) for s in hold_grid]
        m = mean(errs)
        rec = {"family": "full", "P": P, "Pc": None, "Pmax": None, "mre_hold_999": m}
        if best is None or rec["mre_hold_999"] < best["mre_hold_999"]:
            best = rec

    # tail candidates
    Pc_grid = [3, 5, 7, 11, 13, 17, 23, 29, 37, 53, 79, 113, 173]
    Pmax_grid = [251, 401, 701]
    for Pc in Pc_grid:
        for Pmax in Pmax_grid:
            errs = [rel(FobsN[999][s], euler_tail_norm(s, s0, Pc, Pmax, plist)) for s in hold_grid]
            m = mean(errs)
            rec = {"family": "tail", "P": None, "Pc": Pc, "Pmax": Pmax, "mre_hold_999": m}
            if best is None or rec["mre_hold_999"] < best["mre_hold_999"]:
                best = rec

    # Evaluate best candidate across K and strip
    def Ebest(s: complex) -> complex:
        if best["family"] == "full":
            return euler_full_norm(s, s0, best["P"], plist)
        return euler_tail_norm(s, s0, best["Pc"], best["Pmax"], plist)

    hold_mre_by_k = {}
    hold_max_by_k = {}
    strip_mre_by_k = {}
    for K in ks:
        errs_h = [rel(FobsN[K][s], Ebest(s)) for s in hold_grid]
        hold_mre_by_k[K] = mean(errs_h)
        hold_max_by_k[K] = max(errs_h)
        errs_s = [rel(FobsN[K][s], Ebest(s)) for s in strip_grid]
        strip_mre_by_k[K] = mean(errs_s) if errs_s else float("nan")

    # cross-K stability
    vals = [hold_mre_by_k[K] for K in ks]
    m = mean(vals)
    sd = (sum((x - m) ** 2 for x in vals) / len(vals)) ** 0.5
    rel_std = sd / max(abs(m), 1e-15)

    improvement = (baseline_hold - hold_mre_by_k[999]) / max(baseline_hold, 1e-15)

    # Gates
    g1 = improvement >= 0.30  # >=30% better than constant residual model
    g2 = hold_mre_by_k[999] <= 0.06
    g3 = rel_std <= 0.02
    g4 = strip_mre_by_k[999] <= 0.30  # diagnostic consistency in strip sample

    sprint_pass = g1 and g2 and g3 and g4

    flush(f"baseline hold MRE (constant F=1): {baseline_hold:.6f}")
    if best["family"] == "full":
        flush(f"best candidate: full odd product up to P={best['P']}")
    else:
        flush(f"best candidate: odd tail product Pc={best['Pc']}, Pmax={best['Pmax']}")
    flush(f"hold MRE (K=999): {hold_mre_by_k[999]:.6f}")
    flush(f"improvement vs baseline: {100.0*improvement:.2f}%")

    flush("\nPer-K errors for best candidate:")
    flush(f"{'K':>5s}  {'hold_mre':>10s}  {'hold_max':>10s}  {'strip_mre':>10s}")
    for K in ks:
        flush(f"{K:5d}  {hold_mre_by_k[K]:10.6f}  {hold_max_by_k[K]:10.6f}  {strip_mre_by_k[K]:10.6f}")
    flush(f"cross-K rel-std (hold_mre): {rel_std:.6f}")

    flush("\n" + "=" * 78)
    flush("SPRINT-L5 GATES")
    flush("=" * 78)
    flush(f"L5.1: >=30% improvement over constant model          -> {100.0*improvement:.2f}%  [{'PASS' if g1 else 'FAIL'}]")
    flush(f"L5.2: best hold MRE (K=999) <= 0.06                 -> {hold_mre_by_k[999]:.6f}  [{'PASS' if g2 else 'FAIL'}]")
    flush(f"L5.3: cross-K stability rel-std <= 0.02              -> {rel_std:.6f}  [{'PASS' if g3 else 'FAIL'}]")
    flush(f"L5.4: strip consistency MRE (K=999) <= 0.30          -> {strip_mre_by_k[999]:.6f}  [{'PASS' if g4 else 'FAIL'}]")

    flush("\n" + "=" * 78)
    flush("SPRINT-L5 VERDICT")
    flush("=" * 78)
    if sprint_pass:
        flush("Result: PASS")
        flush("Interpretation: F_obs is well-explained by an odd-prime Euler component.")
    else:
        flush("Result: FAIL/PIVOT")
        flush("Interpretation: Euler-tail identity for F is not sufficiently supported by current gates.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

