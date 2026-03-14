#!/usr/bin/env python3
"""
L12: shared operator-core regression checks
==========================================

Regression script for the unified operator/stopping-time core in `_shared.py`.

Checks:
  1) exact low-K stopping-time profiles satisfy
       sum_tau u(tau) == mu == x0^T (I-Q)^(-1) r
  2) parsed high-K profiles satisfy the same identity
  3) record-level channel helper matches explicit u_mix + mu_chi evaluation
  4) unified core bank contains both low-K exact and high-K parsed records
"""

import argparse

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def profile_identity_error(profile: _shared.TauProfile) -> float:
    mu_sum = sum(profile.u.values())
    mu_res = _shared.mu_from_resolvent(profile)
    return max(abs(mu_sum - profile.mu_direct), abs(mu_res - profile.mu_direct))


def main():
    parser = argparse.ArgumentParser(description="Operator-core regression checks.")
    parser.add_argument("--k-exact-min", type=int, default=7)
    parser.add_argument("--k-exact-max", type=int, default=9)
    args = parser.parse_args()

    flush("=" * 78)
    flush("L12: operator-core regression")
    flush("=" * 78)

    chars = {rec["name"]: rec for rec in _shared.primitive_characters()}
    chi4 = chars["chi4"]["fn"]
    s_probe = complex(1.4, 1.0)

    exact_bank = _shared.build_exact_bank(args.k_exact_min, args.k_exact_max)
    high_bank = _shared.build_highk_bank()
    core_bank = _shared.build_core_bank(args.k_exact_min, args.k_exact_max)

    flush("\nExact low-K profile identities")
    flush("-" * 78)
    exact_ok = True
    for K in sorted(exact_bank.keys()):
        for sec in ("00", "10"):
            profile = exact_bank[K][f"profile{sec}"]
            err = profile_identity_error(profile)
            exact_ok = exact_ok and (err <= 1e-12)
            flush(
                f"K={K:2d} sec={sec}: "
                f"mu_direct={profile.mu_direct:+.12e}  "
                f"|identity_err|={err:.3e}"
            )

    flush("\nHigh-K parsed profile identities")
    flush("-" * 78)
    high_ok = True
    for K in (19, 20, 21, 999):
        for sec in ("00", "10"):
            profile = high_bank[K][f"profile{sec}"]
            err = profile_identity_error(profile)
            high_ok = high_ok and (err <= 1e-12)
            flush(
                f"K={K:3d} sec={sec}: "
                f"mu_direct={profile.mu_direct:+.12e}  "
                f"|identity_err|={err:.3e}"
            )

    flush("\nCharacter-channel record helper")
    flush("-" * 78)
    channel_ok = True
    for K in sorted(core_bank.keys()):
        rec = core_bank[K]
        direct = _shared.mu_chi_from_record(rec, chi4, s_probe)
        via_umix = _shared.mu_chi(_shared.u_mix_map(rec), chi4, s_probe)
        err = abs(direct - via_umix)
        channel_ok = channel_ok and (err <= 1e-14)
        flush(f"K={K:3d}: |helper-explicit|={err:.3e}")

    flush("\nUnified bank coverage")
    flush("-" * 78)
    expected_exact = set(range(args.k_exact_min, args.k_exact_max + 1))
    expected_high = {19, 20, 21, 999}
    bank_keys = set(core_bank.keys())
    cover_ok = expected_exact.issubset(bank_keys) and expected_high.issubset(bank_keys)
    flush(f"Exact keys present : {sorted(expected_exact & bank_keys)}")
    flush(f"High keys present  : {sorted(expected_high & bank_keys)}")

    overall = exact_ok and high_ok and channel_ok and cover_ok

    flush("\n" + "=" * 78)
    flush("L12 VERDICT")
    flush("=" * 78)
    flush(f"Exact identities      : {'PASS' if exact_ok else 'FAIL'}")
    flush(f"High-K identities     : {'PASS' if high_ok else 'FAIL'}")
    flush(f"Channel helper        : {'PASS' if channel_ok else 'FAIL'}")
    flush(f"Unified bank coverage : {'PASS' if cover_ok else 'FAIL'}")
    flush(f"\nOVERALL: {'PASS' if overall else 'FAIL'}")
    flush("=" * 78)


if __name__ == "__main__":
    main()
