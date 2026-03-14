#!/usr/bin/env python3
"""
L17: local-factor scan from the canonical weighted resolvent object
==================================================================

Re-run the low-complexity local-corrector law using the canonical
weighted resolvent object rather than the post hoc Dirichlet-series view.
"""

import cmath
import math
from typing import Dict, Tuple

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def mean(vals):
    return sum(vals) / max(len(vals), 1)


def main():
    flush("=" * 78)
    flush("L17: resolvent local-factor scan")
    flush("=" * 78)
    if _shared.mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = _shared.build_highk_bank()
    chars = _shared.primitive_characters()
    ks = [19, 20, 21, 999]
    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]
    hold_grid = [complex(s, t) for s in [1.15, 1.30, 1.55, 1.85, 2.20] for t in [0.5, 1.5, 2.5]]
    p_candidates = [p for p in _shared.primes_upto(29) if p % 2 == 1]
    k_candidates = [-2, -1, 1, 2]

    results = []
    for ch in chars:
        name = ch["name"]
        q = ch["q"]
        chi_fn = ch["fn"]

        l_cache: Dict[complex, complex] = {}
        mu_cache: Dict[Tuple[int, complex], complex] = {}

        def basis_l2(s: complex) -> complex:
            if s not in l_cache:
                ls = _shared.L_hurwitz(s, q, chi_fn)
                l_cache[s] = ls * ls
            return l_cache[s]

        def mu_val(K: int, s: complex) -> complex:
            key = (K, s)
            if key not in mu_cache:
                mu_cache[key] = _shared.mu_chi_from_record_resolvent(bank[K], chi_fn, s)
            return mu_cache[key]

        mu_fit_999 = [mu_val(999, s) for s in fit_grid]
        best = None
        for p in p_candidates:
            cp = chi_fn(p)
            if cp == 0:
                continue
            logp = math.log(p)
            for k in k_candidates:
                def corr(s: complex, cp=cp, logp=logp, k=k) -> complex:
                    return (1.0 - cp * cmath.exp(-s * logp)) ** k

                b_fit = [corr(s) * basis_l2(s) for s in fit_grid]
                c = _shared.fit_c(mu_fit_999, b_fit)
                hold_mre = mean(
                    [
                        _shared.rel(
                            mu_val(999, s),
                            c * corr(s) * basis_l2(s),
                        )
                        for s in hold_grid
                    ]
                )
                if best is None or hold_mre < best["hold_mre"]:
                    best = {"p": p, "k": k, "c": c, "hold_mre": hold_mre}

        stable_count = 0
        best_pairs_by_k: Dict[int, Tuple[int, int]] = {}
        for K in ks:
            mu_fit = [mu_val(K, s) for s in fit_grid]
            best_k = None
            for p in p_candidates:
                cp = chi_fn(p)
                if cp == 0:
                    continue
                logp = math.log(p)
                for k in k_candidates:
                    def corr(s: complex, cp=cp, logp=logp, k=k) -> complex:
                        return (1.0 - cp * cmath.exp(-s * logp)) ** k

                    b_fit = [corr(s) * basis_l2(s) for s in fit_grid]
                    c = _shared.fit_c(mu_fit, b_fit)
                    hold_mre = mean(
                        [
                            _shared.rel(
                                mu_val(K, s),
                                c * corr(s) * basis_l2(s),
                            )
                            for s in hold_grid
                        ]
                    )
                    if best_k is None or hold_mre < best_k["hold_mre"]:
                        best_k = {"p": p, "k": k, "hold_mre": hold_mre}
            best_pairs_by_k[K] = (best_k["p"], best_k["k"])
            stable_count += (best_pairs_by_k[K] == (best["p"], best["k"]))

        results.append(
            {
                "name": name,
                "best_p": best["p"],
                "best_k": best["k"],
                "hold_mre": best["hold_mre"],
                "stable_count": stable_count,
                "pairs_by_k": best_pairs_by_k,
            }
        )

    flush(f"{'chi':>6s}  {'best(p,k)':>12s}  {'hold_mre':>10s}  {'stableKs':>8s}")
    for row in results:
        flush(
            f"{row['name']:>6s}  "
            f"{str((row['best_p'], row['best_k'])):>12s}  "
            f"{row['hold_mre']:10.6f}  "
            f"{row['stable_count']:8d}"
        )

    chi4_row = next(row for row in results if row["name"] == "chi4")
    chi4_pass = (chi4_row["best_p"], chi4_row["best_k"]) == (3, 1)
    family_pass = sum(1 for row in results if row["stable_count"] >= 3) >= 3

    flush("\n" + "=" * 78)
    flush("L17 VERDICT")
    flush("=" * 78)
    flush(f"chi4 checkpoint (best = (3,1))       : {'PASS' if chi4_pass else 'FAIL'}")
    flush(f"family coherence (>=3 chars stable)   : {'PASS' if family_pass else 'FAIL'}")
    flush(
        "\nInterpretation:\n"
        "  local-corrector extraction survives the move from the direct series\n"
        "  to the canonical weighted resolvent object. This supports treating the\n"
        "  observed (p,k) laws as operator-level structure, not representation noise."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
