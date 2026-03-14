#!/usr/bin/env python3
"""
L15: s=1 tail-majorant template for limit/sum exchange
======================================================

Operational version of the tail-control gap needed to justify
exchanging the K->inf limit with the stopping-time sum.
"""

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def envelope_B(bank, sec: str):
    taus = set()
    for K in (19, 20, 21, 999):
        taus.update(bank[K][f"profile{sec}"].u.keys())
    taus = sorted(taus)
    B = {}
    for tau in taus:
        B[tau] = max(abs(bank[K][f"profile{sec}"].u.get(tau, 0.0)) for K in (19, 20, 21, 999))
    return taus, B


def tail_bound_from_envelope(taus, B, T, tau_star, rho_cont):
    observed = 0.0
    for tau in taus:
        if T < tau <= tau_star:
            observed += B.get(tau, 0.0)
    continuation = 0.0
    if T < tau_star:
        continuation = B.get(tau_star, 0.0) * rho_cont / (1.0 - rho_cont)
    return observed + continuation, observed, continuation


def main():
    flush("=" * 78)
    flush("L15: s=1 tail-majorant and limit/sum exchange")
    flush("=" * 78)

    bank = _shared.build_highk_bank()
    rho_cont = 0.56

    for sec in ("00", "10"):
        taus, B = envelope_B(bank, sec)
        tau_star = max(taus)

        flush("\n" + "=" * 78)
        flush(f"Sector {sec}: uniform envelope B(tau)=sup_K |u(tau;K)|")
        flush("=" * 78)
        flush(f"tau range observed: {taus[0]}..{tau_star}")
        flush(f"continuation factor rho_cont = {rho_cont:.3f}")
        flush(f"{'tau':>4s}  {'B(tau)':>14s}")
        for tau in taus[:10]:
            flush(f"{tau:4d}  {B[tau]:14.6e}")
        flush("   ...")
        for tau in taus[-8:]:
            flush(f"{tau:4d}  {B[tau]:14.6e}")

        recent = []
        for tau in range(max(taus[0], tau_star - 7), tau_star):
            a = B.get(tau, 0.0)
            b = B.get(tau + 1, 0.0)
            if a > 1e-18 and b > 0:
                recent.append(b / a)
        if recent:
            recent_sorted = sorted(recent)
            flush(
                f"\nrecent B(t+1)/B(t): median={recent_sorted[len(recent_sorted)//2]:.4f}, "
                f"max={max(recent):.4f}"
            )

        flush("\nTail bounds:")
        flush(f"{'T':>4s}  {'tail_bound':>14s}  {'obs_part':>14s}  {'cont_part':>14s}")
        for T in (8, 10, 12, 14, 16):
            tb, obs, cont = tail_bound_from_envelope(taus, B, T, tau_star, rho_cont)
            flush(f"{T:4d}  {tb:14.6e}  {obs:14.6e}  {cont:14.6e}")

    flush("\n" + "=" * 78)
    flush("Gap-(b) lemma template")
    flush("=" * 78)
    flush(
        "Assume for each sector ab in {00,10}:\n"
        "  (H1) for K in {19,20,21,inf-proxy}, |u_ab(tau;K)| <= B_ab(tau),\n"
        "  (H2) for tau>tau_*: B_ab(tau+1) <= rho_cont * B_ab(tau), rho_cont<1.\n"
        "Then, for every T<tau_*:\n"
        "  sup_K sum_{tau>T} |u_ab(tau;K)|\n"
        "  <= sum_{tau=T+1}^{tau_*} B_ab(tau) + B_ab(tau_*)*rho_cont/(1-rho_cont).\n"
        "Hence tails are uniformly controlled; limit/sum exchange reduces to proving\n"
        "H1-H2 analytically."
    )

    flush("\nEpsilon targets (max over sectors):")
    flush(f"{'T':>4s}  {'epsilon(T)':>14s}")
    for T in (8, 10, 12, 14, 16):
        vals = []
        for sec in ("00", "10"):
            taus, B = envelope_B(bank, sec)
            tau_star = max(taus)
            tb, _, _ = tail_bound_from_envelope(taus, B, T, tau_star, rho_cont)
            vals.append(tb)
        flush(f"{T:4d}  {max(vals):14.6e}")

    flush("\n" + "=" * 78)
    flush("L15 VERDICT")
    flush("=" * 78)
    flush(
        "The limit/sum exchange gap is reduced to an explicit two-hypothesis\n"
        "tail-majorant lemma with numeric epsilon(T) targets. Remaining gap:\n"
        "prove the continuation inequality H2 from conditioned carry dynamics."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
