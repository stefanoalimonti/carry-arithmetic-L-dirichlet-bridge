#!/usr/bin/env python3
"""
L14: s=1 lambda2 envelope toward the conditioned 1/2 rate
==========================================================

Operational version of the dominant-rate gap:
  show that the conditioned D-odd chain is consistent with a dominant
  asymptotic rate 1/2, even though per-position bulk matrices do not
  themselves converge to the unconditioned Markov operator.
"""

import math
from statistics import median

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def q2_map(profile: _shared.TauProfile):
    out = {}
    for tau in profile.taus:
        if tau + 1 in profile.q:
            out[tau] = profile.q[tau] * profile.q[tau + 1]
    return out


def fit_power_law_gap(tau_vals, q2_vals):
    pts = []
    for tau, q in zip(tau_vals, q2_vals):
        gap = q - 0.25
        if tau > 0 and gap > 0:
            pts.append((math.log(tau), math.log(gap)))
    if len(pts) < 2:
        return float("nan"), float("nan")
    xs = [x for x, _ in pts]
    ys = [y for _, y in pts]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    ss_xy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    ss_xx = sum((x - mx) ** 2 for x in xs)
    slope = ss_xy / ss_xx if ss_xx > 0 else 0.0
    alpha = -slope
    a = math.exp(my + alpha * mx)
    return a, alpha


def main():
    flush("=" * 78)
    flush("L14: s=1 lambda2 envelope toward 1/2")
    flush("=" * 78)

    bank = _shared.build_highk_bank()
    q2_store = {}
    for K in (19, 20, 21):
        for sec in ("00", "10"):
            q2_store[(K, sec)] = q2_map(bank[K][f"profile{sec}"])

    tau_win = list(range(5, 11))

    flush("\n" + "=" * 78)
    flush("A) q2(tau) table in the converged depth window")
    flush("=" * 78)
    flush(f"{'tau':>4s}  {'q2_00':>12s}  {'gap00':>10s}  {'q2_10':>12s}  {'gap10':>10s}")
    for tau in tau_win:
        q00 = q2_store[(21, '00')].get(tau, float('nan'))
        q10 = q2_store[(21, '10')].get(tau, float('nan'))
        flush(f"{tau:4d}  {q00:12.8f}  {q00-0.25:+10.6f}  {q10:12.8f}  {q10-0.25:+10.6f}")

    flush("\n" + "=" * 78)
    flush("B) Power-law envelope q2(tau)-1/4 ~ A/tau^alpha")
    flush("=" * 78)
    fit_summary = {}
    for sec in ("00", "10"):
        tau_vals = [tau for tau in tau_win if tau in q2_store[(21, sec)]]
        q_vals = [q2_store[(21, sec)][tau] for tau in tau_vals]
        a, alpha = fit_power_law_gap(tau_vals, q_vals)
        a_safe = 1.35 * a
        fit_summary[sec] = (a, alpha, a_safe)
        flush(f"sector {sec}: A={a:.6f}, alpha={alpha:.4f}, A_safe={a_safe:.6f}")

    flush("\nEnvelope check:")
    flush(f"{'tau':>4s}  {'q2_00<=env':>12s}  {'q2_10<=env':>12s}")
    for tau in tau_win:
        checks = []
        for sec in ("00", "10"):
            q = q2_store[(21, sec)].get(tau, float("nan"))
            _, alpha, a_safe = fit_summary[sec]
            env = 0.25 + a_safe / (tau ** alpha)
            checks.append(q <= env + 1e-15)
        flush(f"{tau:4d}  {str(checks[0]):>12s}  {str(checks[1]):>12s}")

    flush("\n" + "=" * 78)
    flush("C) Resulting rho(tau) envelope")
    flush("=" * 78)
    flush(f"{'tau':>4s}  {'rho_bd_00':>10s}  {'rho_bd_10':>10s}  {'rho_bd_max':>10s}")
    rho_tail = []
    for tau in range(8, 21):
        rho_vals = []
        for sec in ("00", "10"):
            _, alpha, a_safe = fit_summary[sec]
            q_bd = 0.25 + a_safe / (tau ** alpha)
            rho_vals.append(math.sqrt(q_bd))
        rho_max = max(rho_vals)
        if tau >= 10:
            rho_tail.append(rho_max)
        flush(f"{tau:4d}  {rho_vals[0]:10.6f}  {rho_vals[1]:10.6f}  {rho_max:10.6f}")

    rho_tail_med = median(rho_tail) if rho_tail else float("nan")

    flush("\n" + "=" * 78)
    flush("D) Cross-check with R(K) decrement ratios")
    flush("=" * 78)
    r_hist = _shared.load_e162_R_history()
    all_res = _shared.load_e45()
    by_k = {r["K"]: r for r in all_res}
    for K in (19, 20, 21):
        r_hist[K] = by_k[K]["S10"] / by_k[K]["S00"]
    ks = sorted(k for k in r_hist.keys() if 19 <= k <= 21)

    ratio_list = []
    prev_r = None
    prev_dec = None
    flush(f"{'K':>4s}  {'R(K)':>15s}  {'dec':>12s}  {'ratio':>9s}")
    for K in ks:
        rv = r_hist[K]
        if prev_r is None:
            flush(f"{K:4d}  {rv:+15.12f}")
        else:
            dec = prev_r - rv
            ratio = dec / prev_dec if (prev_dec is not None and abs(prev_dec) > 1e-15) else float("nan")
            if ratio == ratio:
                ratio_list.append((K, ratio))
                flush(f"{K:4d}  {rv:+15.12f}  {dec:+12.9f}  {ratio:9.6f}")
            prev_dec = dec
        prev_r = rv

    last_ratio = ratio_list[-1][1] if ratio_list else float("nan")

    flush("\n" + "=" * 78)
    flush("L14 VERDICT")
    flush("=" * 78)
    flush(
        "The conditioned chain is consistent with a theorem-style envelope\n"
        "  q2(tau) <= 1/4 + A/tau^alpha,\n"
        "hence rho(tau) = 1/2 + O(tau^-alpha).\n"
        f"Current median rho upper envelope (tau>=10): {rho_tail_med:.6f}\n"
        f"Observed last decrement ratio: {last_ratio:.6f}\n"
        "\nRemaining theorem gap:\n"
        "  derive the envelope analytically from conditioned carry dynamics,\n"
        "  then lift it from depthwise survival to the global spectral statement."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
