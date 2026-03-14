#!/usr/bin/env python3
"""
P1-SPRINT-L5bis: full audit for residual identity models
=========================================================

Comprehensive tests on the chi4 residual factor:
  F_obs(s) = mu(s)/(c_K * L(s,chi4)), normalized at anchor s0.

Blocks:
  A) Model-suite comparison (M0/M1/M2/M3)
  B) Bootstrap comparison M1 vs M2
  C) Prime ablation (leave-one-out)
  D) Negative controls (sign-shuffled Euler factors)
  E) Cross-character L5-style summary

Model family (chi4):
  M0: constant 1
  M1: odd-prime tail Euler product (selected on fit grid only)
  M2: (1+3^{-s})*L(s,chi4) (normalized)
  M3: L(s,chi4) (normalized)
"""

import cmath
import math
import random
from typing import Callable, Dict, List, Tuple

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def chi4(n: int) -> int:
    r = n % 4
    if r % 2 == 0:
        return 0
    return 1 if r == 1 else -1


primes_upto = _shared.primes_upto
build_highk_bank = _shared.build_highk_bank
u_mix_map = _shared.u_mix_map
n_of_tau = _shared.n_of_tau
primitive_characters = _shared.primitive_characters
mu_chi = _shared.mu_chi
L_hurwitz = _shared.L_hurwitz
fit_c = _shared.fit_c
rel = _shared.rel
mean = _shared.mean


def euler_tail_norm_generic(
    s: complex,
    s0: complex,
    Pc: int,
    Pmax: int,
    plist: List[int],
    chi_fn: Callable[[int], int],
) -> complex:
    z = 1.0 + 0.0j
    z0 = 1.0 + 0.0j
    for p in plist:
        if p <= Pc or p > Pmax:
            continue
        cp = chi_fn(p)
        if cp == 0:
            continue
        z *= 1.0 / (1.0 - cp * cmath.exp(-s * math.log(p)))
        z0 *= 1.0 / (1.0 - cp * cmath.exp(-s0 * math.log(p)))
    return z / z0 if abs(z0) > 1e-15 else complex(float("nan"), 0.0)


def model_suite_chi4(
    FobsN_by_k: Dict[int, Dict[complex, complex]],
    fit_grid: List[complex],
    hold_grid: List[complex],
    strip_grid: List[complex],
    s0: complex,
    plist: List[int],
):
    # M1 selection on fit-grid only
    Pc_grid = [3, 5, 7, 11, 13, 17, 23, 29, 37, 53, 79, 113, 173]
    Pmax_grid = [251, 401, 701]
    best = None
    for Pc in Pc_grid:
        for Pmax in Pmax_grid:
            errs = [rel(FobsN_by_k[999][s], euler_tail_norm_generic(s, s0, Pc, Pmax, plist, chi4)) for s in fit_grid]
            m = mean(errs)
            rec = {"Pc": Pc, "Pmax": Pmax, "fit_mre": m}
            if best is None or rec["fit_mre"] < best["fit_mre"]:
                best = rec

    Pc_star, Pmax_star = best["Pc"], best["Pmax"]

    # model functions (normalized)
    def M0(s):
        return 1.0 + 0.0j

    def M1(s):
        return euler_tail_norm_generic(s, s0, Pc_star, Pmax_star, plist, chi4)

    def M2(s):
        # ((1+3^-s)L(s))/((1+3^-s0)L(s0))
        num = (1.0 + cmath.exp(-s * math.log(3))) * L_hurwitz(s, 4, chi4)
        den = (1.0 + cmath.exp(-s0 * math.log(3))) * L_hurwitz(s0, 4, chi4)
        return num / den if abs(den) > 1e-15 else complex(float("nan"), 0.0)

    def M3(s):
        num = L_hurwitz(s, 4, chi4)
        den = L_hurwitz(s0, 4, chi4)
        return num / den if abs(den) > 1e-15 else complex(float("nan"), 0.0)

    models = {"M0_const1": M0, "M1_tail": M1, "M2_(1+3^-s)L": M2, "M3_L": M3}
    metrics = {}
    for mname, mf in models.items():
        hold_mre = {}
        hold_max = {}
        strip_mre = {}
        for K in [19, 20, 21, 999]:
            eh = [rel(FobsN_by_k[K][s], mf(s)) for s in hold_grid]
            es = [rel(FobsN_by_k[K][s], mf(s)) for s in strip_grid]
            hold_mre[K] = mean(eh)
            hold_max[K] = max(eh)
            strip_mre[K] = mean(es)
        metrics[mname] = {"hold_mre": hold_mre, "hold_max": hold_max, "strip_mre": strip_mre}

    return {
        "Pc_star": Pc_star,
        "Pmax_star": Pmax_star,
        "fit_mre_star": best["fit_mre"],
        "models": models,
        "metrics": metrics,
    }


def direct_mu_model_compare_chi4(
    u_bank_chi4: Dict[int, Dict[int, float]],
    fit_grid: List[complex],
    hold_grid: List[complex],
    s0: complex,
    plist: List[int],
    Pc_star: int,
    Pmax_star: int,
):
    def basis_eval(name: str, s: complex) -> complex:
        Ls = L_hurwitz(s, 4, chi4)
        if name == "B1_L":
            return Ls
        if name == "B2_L2":
            return Ls * Ls
        if name == "B3_(1+3^-s)L2":
            return (1.0 + cmath.exp(-s * math.log(3))) * Ls * Ls
        if name == "B4_L_times_tail":
            return Ls * euler_tail_norm_generic(s, s0, Pc_star, Pmax_star, plist, chi4)
        raise ValueError(name)

    basis_names = ["B1_L", "B2_L2", "B3_(1+3^-s)L2", "B4_L_times_tail"]
    out = {}
    for bname in basis_names:
        fit_err = {}
        hold_err = {}
        for K in [19, 20, 21, 999]:
            mu_fit = [mu_chi(u_bank_chi4[K], chi4, s) for s in fit_grid]
            bas_fit = [basis_eval(bname, s) for s in fit_grid]
            c = fit_c(mu_fit, bas_fit)

            ef = [rel(mu_chi(u_bank_chi4[K], chi4, s), c * basis_eval(bname, s)) for s in fit_grid]
            eh = [rel(mu_chi(u_bank_chi4[K], chi4, s), c * basis_eval(bname, s)) for s in hold_grid]
            fit_err[K] = mean(ef)
            hold_err[K] = mean(eh)
        out[bname] = {"fit_mre": fit_err, "hold_mre": hold_err}
    return out


def bootstrap_delta(
    FobsN: Dict[complex, complex],
    model_a: Callable[[complex], complex],
    model_b: Callable[[complex], complex],
    hold_grid: List[complex],
    n_boot: int = 1000,
    seed: int = 7,
) -> Dict[str, float]:
    rng = random.Random(seed)
    deltas = []
    n = len(hold_grid)
    for _ in range(n_boot):
        sample = [hold_grid[rng.randrange(n)] for _ in range(n)]
        ea = [rel(FobsN[s], model_a(s)) for s in sample]
        eb = [rel(FobsN[s], model_b(s)) for s in sample]
        deltas.append(mean(eb) - mean(ea))  # positive => A better
    deltas.sort()

    def pct(p):
        idx = int(round((len(deltas) - 1) * p))
        return deltas[max(0, min(len(deltas) - 1, idx))]

    return {"lo": pct(0.025), "med": pct(0.5), "hi": pct(0.975)}


def ablation_leave_one_out(
    FobsN: Dict[complex, complex],
    hold_grid: List[complex],
    s0: complex,
    Pc: int,
    Pmax: int,
    plist: List[int],
) -> List[Tuple[int, float, float]]:
    # baseline M1 error
    def M1(s):
        return euler_tail_norm_generic(s, s0, Pc, Pmax, plist, chi4)

    base = mean([rel(FobsN[s], M1(s)) for s in hold_grid])

    support_primes = [p for p in plist if (p > Pc and p <= Pmax and chi4(p) != 0)]
    out = []
    for p_drop in support_primes:
        def Mdrop(s, p_drop=p_drop):
            z = 1.0 + 0.0j
            z0 = 1.0 + 0.0j
            for p in support_primes:
                if p == p_drop:
                    continue
                cp = chi4(p)
                z *= 1.0 / (1.0 - cp * cmath.exp(-s * math.log(p)))
                z0 *= 1.0 / (1.0 - cp * cmath.exp(-s0 * math.log(p)))
            return z / z0 if abs(z0) > 1e-15 else complex(float("nan"), 0.0)

        e = mean([rel(FobsN[s], Mdrop(s)) for s in hold_grid])
        out.append((p_drop, e, e - base))

    out.sort(key=lambda x: x[2], reverse=True)
    return out


def sign_shuffle_controls(
    FobsN: Dict[complex, complex],
    hold_grid: List[complex],
    s0: complex,
    Pc: int,
    Pmax: int,
    plist: List[int],
    n_draws: int = 1200,
    seed: int = 13,
) -> Dict[str, float]:
    rng = random.Random(seed)
    support = [p for p in plist if p > Pc and p <= Pmax and chi4(p) != 0]
    true_signs = [chi4(p) for p in support]
    n_pos = sum(1 for x in true_signs if x > 0)
    n_neg = len(true_signs) - n_pos

    def eval_with_signs(signs: List[int]) -> float:
        sign_map = {p: s for p, s in zip(support, signs)}

        def M(s):
            z = 1.0 + 0.0j
            z0 = 1.0 + 0.0j
            for p in support:
                cp = sign_map[p]
                z *= 1.0 / (1.0 - cp * cmath.exp(-s * math.log(p)))
                z0 *= 1.0 / (1.0 - cp * cmath.exp(-s0 * math.log(p)))
            return z / z0 if abs(z0) > 1e-15 else complex(float("nan"), 0.0)

        return mean([rel(FobsN[s], M(s)) for s in hold_grid])

    true_mre = eval_with_signs(true_signs)
    draws = []
    base_signs = [1] * n_pos + [-1] * n_neg
    for _ in range(n_draws):
        perm = base_signs[:]
        rng.shuffle(perm)
        draws.append(eval_with_signs(perm))
    draws.sort()
    better_or_equal = sum(1 for x in draws if x <= true_mre)
    p_emp = (better_or_equal + 1) / (len(draws) + 1)

    return {
        "true_mre": true_mre,
        "draw_med": draws[len(draws) // 2],
        "draw_lo": draws[int(round(0.025 * (len(draws) - 1)))],
        "draw_hi": draws[int(round(0.975 * (len(draws) - 1)))] ,
        "p_emp": p_emp,
    }


def cross_character_l5_summary(
    bank: Dict[int, dict],
    fit_grid: List[complex],
    hold_grid: List[complex],
    s0: complex,
    plist: List[int],
):
    chars = primitive_characters()
    rows = []
    for ch in chars:
        name, q, fn = ch["name"], ch["q"], ch["fn"]
        u_bank = {K: u_mix_map(bank[K]) for K in [19, 20, 21, 999]}

        # cK from fit
        cK = {}
        FobsN = {K: {} for K in [19, 20, 21, 999]}
        for K in [19, 20, 21, 999]:
            mu_fit = [mu_chi(u_bank[K], fn, s) for s in fit_grid]
            l_fit = [L_hurwitz(s, q, fn) for s in fit_grid]
            c = fit_c(mu_fit, l_fit)
            cK[K] = c
            f0 = mu_chi(u_bank[K], fn, s0) / (c * L_hurwitz(s0, q, fn))
            for s in fit_grid + hold_grid:
                den = c * L_hurwitz(s, q, fn)
                f = mu_chi(u_bank[K], fn, s) / den if abs(den) > 1e-15 else complex(float("nan"), 0.0)
                FobsN[K][s] = f / f0 if abs(f0) > 1e-15 else complex(float("nan"), 0.0)

        base = mean([rel(FobsN[999][s], 1.0 + 0.0j) for s in hold_grid])

        best = None
        for Pc in [3, 5, 7, 11, 13, 17, 23, 29, 37, 53]:
            for Pmax in [251, 401, 701]:
                def M(s, Pc=Pc, Pmax=Pmax):
                    return euler_tail_norm_generic(s, s0, Pc, Pmax, plist, fn)
                fit_mre = mean([rel(FobsN[999][s], M(s)) for s in fit_grid])
                hold_mre = mean([rel(FobsN[999][s], M(s)) for s in hold_grid])
                rec = {"Pc": Pc, "Pmax": Pmax, "fit_mre": fit_mre, "hold_mre": hold_mre}
                if best is None or rec["fit_mre"] < best["fit_mre"]:
                    best = rec

        imp = (base - best["hold_mre"]) / max(base, 1e-15)
        rows.append((name, base, best["hold_mre"], imp, best["Pc"], best["Pmax"]))
    return rows


def main():
    flush("=" * 78)
    flush("P1-SPRINT-L5bis: full audit")
    flush("=" * 78)
    if _shared.mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = build_highk_bank()
    plist = primes_upto(701)
    ks = [19, 20, 21, 999]

    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]
    hold_grid = [complex(s, t) for s in [1.15, 1.30, 1.55, 1.85, 2.20] for t in [0.5, 1.5, 2.5]]
    strip_grid = [complex(s, t) for s in [0.90, 0.80, 0.70, 0.60, 0.55] for t in [0.0, 1.0, 2.0, 3.0, 4.0]]
    strip_grid = [s for s in strip_grid if abs(L_hurwitz(s, 4, chi4)) > 1e-6]
    s0 = complex(1.70, 0.0)

    # chi4 FobsN
    u_bank_chi4 = {K: u_mix_map(bank[K]) for K in ks}
    c_by_k = {}
    FobsN = {K: {} for K in ks}
    for K in ks:
        mu_fit = [mu_chi(u_bank_chi4[K], chi4, s) for s in fit_grid]
        l_fit = [L_hurwitz(s, 4, chi4) for s in fit_grid]
        cK = fit_c(mu_fit, l_fit)
        c_by_k[K] = cK
        f0 = mu_chi(u_bank_chi4[K], chi4, s0) / (cK * L_hurwitz(s0, 4, chi4))
        for s in fit_grid + hold_grid + strip_grid:
            den = cK * L_hurwitz(s, 4, chi4)
            f = mu_chi(u_bank_chi4[K], chi4, s) / den if abs(den) > 1e-15 else complex(float("nan"), 0.0)
            FobsN[K][s] = f / f0 if abs(f0) > 1e-15 else complex(float("nan"), 0.0)

    suite = model_suite_chi4(FobsN, fit_grid, hold_grid, strip_grid, s0, plist)
    models = suite["models"]
    mets = suite["metrics"]
    Pc_star = suite["Pc_star"]
    Pmax_star = suite["Pmax_star"]

    flush(f"chi4 M1 selected on fit-grid: Pc={Pc_star}, Pmax={Pmax_star}, fit-MRE={suite['fit_mre_star']:.6f}")
    flush(f"{'model':>14s}  {'hold_mre@999':>12s}  {'hold_max@999':>12s}  {'strip_mre@999':>12s}")
    for mname in ["M0_const1", "M1_tail", "M2_(1+3^-s)L", "M3_L"]:
        flush(
            f"{mname:>14s}  {mets[mname]['hold_mre'][999]:12.6f}  "
            f"{mets[mname]['hold_max'][999]:12.6f}  {mets[mname]['strip_mre'][999]:12.6f}"
        )

    # bootstrap M1 vs M2 at K=999
    boot = bootstrap_delta(FobsN[999], models["M1_tail"], models["M2_(1+3^-s)L"], hold_grid, n_boot=1200, seed=11)
    flush(
        f"bootstrap delta (M2 - M1) on hold MRE: "
        f"median={boot['med']:.6f}, 95%CI=[{boot['lo']:.6f},{boot['hi']:.6f}]"
    )

    # ablation
    abl = ablation_leave_one_out(FobsN[999], hold_grid, s0, Pc_star, Pmax_star, plist)
    flush("top prime ablations (largest MRE increase):")
    for p, e, de in abl[:8]:
        flush(f"  p={p:3d}: hold_mre={e:.6f}, delta={de:+.6f}")

    # negative controls
    ctrl = sign_shuffle_controls(FobsN[999], hold_grid, s0, Pc_star, Pmax_star, plist, n_draws=1200, seed=13)
    flush(
        "shuffle controls: "
        f"true_mre={ctrl['true_mre']:.6f}, median_shuffle={ctrl['draw_med']:.6f}, "
        f"95% band=[{ctrl['draw_lo']:.6f},{ctrl['draw_hi']:.6f}], p_emp={ctrl['p_emp']:.4f}"
    )

    # direct mu-model comparison
    direct = direct_mu_model_compare_chi4(u_bank_chi4, fit_grid, hold_grid, s0, plist, Pc_star, Pmax_star)
    flush("\nchi4 direct mu-model comparison (hold MRE @ K=999):")
    for name in ["B1_L", "B2_L2", "B3_(1+3^-s)L2", "B4_L_times_tail"]:
        flush(f"  {name:18s}: {direct[name]['hold_mre'][999]:.6f}")

    # cross-character
    rows = cross_character_l5_summary(bank, fit_grid, hold_grid, s0, plist)
    flush("\ncross-character odd-tail summary:")
    flush(f"{'chi':>6s}  {'base':>10s}  {'best_tail':>10s}  {'impr%':>8s}  {'Pc':>4s}  {'Pmax':>5s}")
    for name, base, besth, imp, pc, pm in rows:
        flush(f"{name:>6s}  {base:10.6f}  {besth:10.6f}  {100.0*imp:8.2f}  {pc:4d}  {pm:5d}")

    # gates
    g1 = mets["M1_tail"]["hold_mre"][999] <= min(mets["M2_(1+3^-s)L"]["hold_mre"][999], mets["M3_L"]["hold_mre"][999]) + 1e-9
    g2 = boot["lo"] > 0.0  # M2-M1 >0 => M1 better with CI
    g3 = ctrl["p_emp"] <= 0.05
    g4 = sum(1 for _, base, besth, imp, _, _ in rows if imp >= 0.25) >= 3

    flush("\n" + "=" * 78)
    flush("SPRINT-L5bis GATES")
    flush("=" * 78)
    flush(f"L5bis.1: M1 <= min(M2,M3) on hold@999                -> {'PASS' if g1 else 'FAIL'}")
    flush(f"L5bis.2: bootstrap CI confirms M1 better than M2     -> {'PASS' if g2 else 'FAIL'}")
    flush(f"L5bis.3: shuffle controls reject random-sign model    -> p={ctrl['p_emp']:.4f}  [{'PASS' if g3 else 'FAIL'}]")
    flush(f"L5bis.4: >=3 characters with >=25% odd-tail gain      -> {'PASS' if g4 else 'FAIL'}")

    verdict = "PASS" if (g1 and g2 and g3 and g4) else "FAIL/PIVOT"
    flush("\n" + "=" * 78)
    flush("SPRINT-L5bis VERDICT")
    flush("=" * 78)
    flush(f"Result: {verdict}")
    if verdict == "PASS":
        flush("Residual structure is robust, non-random, and supports odd-prime Euler interpretation.")
    else:
        flush("Residual structure not fully robust under full-audit controls.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

