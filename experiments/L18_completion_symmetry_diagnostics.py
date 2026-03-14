#!/usr/bin/env python3
"""
L18: completion / functional-equation diagnostics for the canonical object
==========================================================================

Search for simple gamma-like completions of the canonical weighted resolvent
object and test whether any of them exhibit approximate symmetry

    Xi(s) ~= eps * Xi(1-s)

on a mirrored strip sample.
"""

import math

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def chi_parity(chi_fn, q: int) -> int:
    val = chi_fn(q - 1)
    return 0 if val == 1 else 1


def gamma_completion_factor(s: complex, q: int, a: int) -> complex:
    ss = _shared.mp.mpc(s.real, s.imag)
    return complex(_shared.mp.power(_shared.mp.pi / q, -(ss + a) / 2) * _shared.mp.gamma((ss + a) / 2))


def fit_epsilon(lhs, rhs):
    num = 0.0 + 0.0j
    den = 0.0
    for x, y in zip(lhs, rhs):
        num += x * y.conjugate()
        den += (y.real * y.real + y.imag * y.imag)
    return num / den if den > 1e-30 else complex(float("nan"), float("nan"))


def mean(vals):
    return sum(vals) / max(len(vals), 1)


def main():
    flush("=" * 78)
    flush("L18: completion symmetry diagnostics")
    flush("=" * 78)
    if _shared.mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = _shared.build_highk_bank()
    chars = _shared.primitive_characters()
    sample = [complex(sig, t) for sig in [0.60, 0.75, 1.20, 1.50] for t in [0.5, 1.5, 3.0]]

    for ch in chars:
        name = ch["name"]
        q = ch["q"]
        chi_fn = ch["fn"]
        parity = chi_parity(chi_fn, q)

        def M(s: complex) -> complex:
            return _shared.mu_chi_from_record_resolvent(bank[999], chi_fn, s)

        candidates = {
            "raw": lambda s: M(s),
            "gamma_even": lambda s: gamma_completion_factor(s, q, 0) * M(s),
            "gamma_odd": lambda s: gamma_completion_factor(s, q, 1) * M(s),
            "gamma_parity": lambda s, parity=parity: gamma_completion_factor(s, q, parity) * M(s),
        }

        best_name = None
        best_score = None
        best_eps = None
        rows = []
        for cname, Xi in candidates.items():
            lhs = [Xi(s) for s in sample]
            rhs = [Xi(1.0 - s) for s in sample]
            eps = fit_epsilon(lhs, rhs)
            errs = [_shared.rel(a, eps * b) for a, b in zip(lhs, rhs)]
            score = mean(errs)
            rows.append((cname, eps, score))
            if best_score is None or score < best_score:
                best_name = cname
                best_score = score
                best_eps = eps

        flush("\n" + "=" * 78)
        flush(f"Character: {name} (parity a={parity})")
        flush("=" * 78)
        for cname, eps, score in rows:
            flush(f"{cname:>12s}: mean symmetry defect = {score:.6f}   eps={eps.real:+.4f}{eps.imag:+.4f}i")
        flush(f"best completion: {best_name}   defect={best_score:.6f}")
        flush(f"simple FE candidate found (<0.10 defect): {'YES' if best_score < 0.10 else 'NO'}")

    flush("\n" + "=" * 78)
    flush("L18 VERDICT")
    flush("=" * 78)
    flush(
        "This diagnostic tests whether a simple classical gamma-completion already\n"
        "puts the canonical weighted resolvent object near a functional equation.\n"
        "If all defects stay large, the current carry object captures amplitude\n"
        "structure but still lacks the right phase/completion mechanism."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
