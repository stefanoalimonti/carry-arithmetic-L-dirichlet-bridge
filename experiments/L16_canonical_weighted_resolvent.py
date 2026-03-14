#!/usr/bin/env python3
"""
L16: canonical analytic object candidate via weighted first-return resolvent
===========================================================================

Define the candidate carry-side analytic object

    M_chi(K, s) = omega(K) * M_10(K, s) - M_00(K, s)

with

    M_ab(K, s) = x0^T (I - Q_ab(K))^(-1) r_ab,chi(s),
    r_ab,chi(s)[tau] = a_ab(tau) E[val|tau,ab] chi(n(tau)) n(tau)^(-s).

This is the intrinsic operator realization of the mixed character channel.
The script verifies that it reproduces the direct Dirichlet-series definition
exactly on the current finite / extrapolated profiles.
"""

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def main():
    flush("=" * 78)
    flush("L16: canonical weighted resolvent object")
    flush("=" * 78)

    bank = _shared.build_core_bank(7, 9)
    chars = _shared.primitive_characters()
    s_grid = [
        complex(1.20, 0.0),
        complex(1.40, 1.0),
        complex(1.70, 2.0),
        complex(0.80, 0.5),
        complex(0.60, 3.0),
    ]

    overall_ok = True
    for ch in chars:
        flush("\n" + "=" * 78)
        flush(f"Character: {ch['name']}")
        flush("=" * 78)
        for K in sorted(bank.keys()):
            errs = []
            for s in s_grid:
                direct = _shared.mu_chi_from_record(bank[K], ch["fn"], s)
                canonical = _shared.mu_chi_from_record_resolvent(bank[K], ch["fn"], s)
                errs.append(abs(direct - canonical))
            max_err = max(errs) if errs else 0.0
            overall_ok = overall_ok and (max_err <= 1e-12)
            flush(f"K={K:3d}: max|series - resolvent| = {max_err:.3e}")

    flush("\n" + "=" * 78)
    flush("L16 VERDICT")
    flush("=" * 78)
    flush(
        "Candidate canonical object defined:\n"
        "  M_chi(K,s) = omega(K) M_10(K,s) - M_00(K,s),\n"
        "  with M_ab(K,s) realized as a weighted first-return resolvent.\n"
        f"Direct series equivalence on tested grids: {'PASS' if overall_ok else 'FAIL'}.\n"
        "\nInterpretation:\n"
        "  the character channel is no longer only a post hoc Dirichlet sum;\n"
        "  it is realized intrinsically as an s-dependent operator quantity."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
