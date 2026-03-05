#!/usr/bin/env python3
"""
P1-STEP3A: Base-2 carry <-> Witt correction identity
====================================================

Goal (Step 3A):
  Build a clean, machine-checkable algebraic core for base 2:

  1) Canonical truncated Witt embedding is a ring homomorphism modulo 2^n.
  2) The classical n-step carry vector of integer addition is in bijection
     with the Witt correction vector Delta = s - x - y (coordinate-wise).

This is the precise "carry = Witt correction terms" statement used as Step-3A
foundation before any Frobenius/trace claims.
"""

import argparse
import math
from collections import defaultdict


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def ghost_component(a, p, k):
    return sum((p ** i) * (a[i] ** (p ** (k - i))) for i in range(k + 1))


def from_ghost(ghosts, p, n):
    """Inverse ghost map for truncated length n."""
    s = []
    for k in range(n):
        rhs = ghosts[k] - sum((p ** i) * (s[i] ** (p ** (k - i))) for i in range(k))
        s.append((rhs // (p ** k)) % p)
    return tuple(s)


def int_to_witt(z, p, n):
    """
    Canonical truncated embedding Z/(p^n) -> W_n(F_p) via constant ghosts.
    """
    z_mod = z % (p ** n)
    return from_ghost([z_mod] * n, p, n)


def witt_add(x, y, p, n):
    gs = [ghost_component(x, p, k) + ghost_component(y, p, k) for k in range(n)]
    return from_ghost(gs, p, n)


def witt_mul(x, y, p, n):
    gs = [ghost_component(x, p, k) * ghost_component(y, p, k) for k in range(n)]
    return from_ghost(gs, p, n)


def carry_vector_add(a, b, n, p=2):
    """
    Classical schoolbook carry vector c_1..c_n for n digits (incoming carry=0).
    c_j is carry entering position j (1-indexed in this doc; 0-index in tuple).
    """
    c = 0
    out = []
    for j in range(n):
        aj = (a >> j) & 1
        bj = (b >> j) & 1
        s = aj + bj + c
        c = s // p
        out.append(c)
    return tuple(out)


def delta_vector_witt_add(a, b, n, p=2):
    """
    Delta = s - x - y in Witt coordinates, where x=iota(a), y=iota(b), s=x+y.
    """
    x = int_to_witt(a, p, n)
    y = int_to_witt(b, p, n)
    s = witt_add(x, y, p, n)
    return tuple(s[k] - x[k] - y[k] for k in range(n))


def run_for_n(n):
    p = 2
    mod = p ** n
    all_vals = range(mod)

    # A) Ring-homomorphism checks modulo 2^n
    add_fail = 0
    mul_fail = 0

    # B) Carry <-> Delta map checks
    c_to_d = defaultdict(set)
    d_to_c = defaultdict(set)

    pairs = 0
    for a in all_vals:
        x = int_to_witt(a, p, n)
        for b in all_vals:
            y = int_to_witt(b, p, n)
            pairs += 1

            # Homomorphism add/mul modulo 2^n
            s_w = witt_add(x, y, p, n)
            s_i = int_to_witt((a + b) % mod, p, n)
            if s_w != s_i:
                add_fail += 1

            m_w = witt_mul(x, y, p, n)
            m_i = int_to_witt((a * b) % mod, p, n)
            if m_w != m_i:
                mul_fail += 1

            c = carry_vector_add(a, b, n, p)
            d = tuple(s_w[k] - x[k] - y[k] for k in range(n))
            c_to_d[c].add(d)
            d_to_c[d].add(c)

    max_d_per_c = max(len(v) for v in c_to_d.values())
    max_c_per_d = max(len(v) for v in d_to_c.values())
    bijection = (max_d_per_c == 1 and max_c_per_d == 1 and len(c_to_d) == len(d_to_c))

    # Full coverage of carry states (all binary n-tuples)
    expected_states = 2 ** n
    carry_coverage_ok = len(c_to_d) == expected_states

    return {
        "n": n,
        "pairs": pairs,
        "add_fail": add_fail,
        "mul_fail": mul_fail,
        "carry_states": len(c_to_d),
        "delta_states": len(d_to_c),
        "expected_states": expected_states,
        "carry_coverage_ok": carry_coverage_ok,
        "max_d_per_c": max_d_per_c,
        "max_c_per_d": max_c_per_d,
        "bijection": bijection,
        "c_to_d": c_to_d,
    }


def main():
    parser = argparse.ArgumentParser(description="Step3A Witt/carry identity checks.")
    parser.add_argument("--n-min", type=int, default=4)
    parser.add_argument("--n-max", type=int, default=9)
    parser.add_argument("--show-map-n", type=int, default=4,
                        help="Print explicit carry->Delta map for this n.")
    args = parser.parse_args()

    flush("=" * 78)
    flush("P1-STEP3A: BASE-2 CARRY <-> WITT CORRECTION IDENTITY")
    flush("=" * 78)
    flush(f"Range: n={args.n_min}..{args.n_max} (modulus 2^n)")

    results = []
    for n in range(args.n_min, args.n_max + 1):
        res = run_for_n(n)
        results.append(res)
        flush("\n" + "-" * 78)
        flush(f"n={n}, pairs={res['pairs']:,}")
        flush(f"  Homomorphism add failures: {res['add_fail']}")
        flush(f"  Homomorphism mul failures: {res['mul_fail']}")
        flush(f"  Carry states seen: {res['carry_states']} / {res['expected_states']}")
        flush(f"  Delta states seen: {res['delta_states']} / {res['expected_states']}")
        flush(f"  max |Delta-set per carry|: {res['max_d_per_c']}")
        flush(f"  max |carry-set per Delta|: {res['max_c_per_d']}")
        flush(f"  Carry coverage complete: {'YES' if res['carry_coverage_ok'] else 'NO'}")
        flush(f"  Carry <-> Delta bijection: {'YES' if res['bijection'] else 'NO'}")

    # Print explicit map for one small n
    target = None
    for r in results:
        if r["n"] == args.show_map_n:
            target = r
            break
    if target is not None:
        flush("\n" + "=" * 78)
        flush(f"Explicit carry -> Delta map (n={target['n']})")
        flush("=" * 78)
        items = sorted(target["c_to_d"].items(), key=lambda kv: kv[0])
        for c, dset in items[:32]:
            d = next(iter(dset))
            flush(f"  c={c}  ->  Delta={d}")
        if len(items) > 32:
            flush(f"  ... ({len(items) - 32} more states)")

    # Global verdict
    all_add_ok = all(r["add_fail"] == 0 for r in results)
    all_mul_ok = all(r["mul_fail"] == 0 for r in results)
    all_bij_ok = all(r["bijection"] for r in results)
    all_cov_ok = all(r["carry_coverage_ok"] for r in results)

    flush("\n" + "=" * 78)
    flush("STEP3A VERDICT")
    flush("=" * 78)
    flush(f"Ring homomorphism modulo 2^n (add): {'PASS' if all_add_ok else 'FAIL'}")
    flush(f"Ring homomorphism modulo 2^n (mul): {'PASS' if all_mul_ok else 'FAIL'}")
    flush(f"Carry-state full coverage: {'PASS' if all_cov_ok else 'FAIL'}")
    flush(f"Carry <-> Delta bijection: {'PASS' if all_bij_ok else 'FAIL'}")

    if all_add_ok and all_mul_ok and all_cov_ok and all_bij_ok:
        flush("\nConclusion:")
        flush("  In base 2, classical addition carries are exactly encoded by")
        flush("  Witt-coordinate correction vectors under the canonical embedding.")
        flush("  This provides a clean algebraic Step-3A bridge.")

    flush("\n" + "=" * 78)
    flush("DONE")
    flush("=" * 78)


if __name__ == "__main__":
    main()

