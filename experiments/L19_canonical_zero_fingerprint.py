#!/usr/bin/env python3
"""
L19: zero fingerprint for the canonical weighted resolvent object
=================================================================

Re-run the critical-strip zero test on the canonical weighted resolvent
object for chi4 and compare its in-box zero structure against L(s, chi4).
"""

import cmath
import math

import _shared


def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def boundary_points(s1: float, s2: float, t1: float, t2: float, n: int = 250):
    pts = []
    for i in range(n):
        x = s1 + (s2 - s1) * i / (n - 1)
        pts.append(complex(x, t1))
    for i in range(1, n):
        y = t1 + (t2 - t1) * i / (n - 1)
        pts.append(complex(s2, y))
    for i in range(1, n):
        x = s2 - (s2 - s1) * i / (n - 1)
        pts.append(complex(x, t2))
    for i in range(1, n - 1):
        y = t2 - (t2 - t1) * i / (n - 1)
        pts.append(complex(s1, y))
    return pts


def count_box_winding(f, s1: float, s2: float, t1: float, t2: float, n_edge: int = 250):
    pts = boundary_points(s1, s2, t1, t2, n=n_edge)
    vals = [f(z) for z in pts]
    min_mod = min(abs(v) for v in vals)
    total = 0.0
    prev = cmath.phase(vals[0])
    for v in vals[1:]:
        cur = cmath.phase(v)
        d = cur - prev
        while d <= -math.pi:
            d += 2.0 * math.pi
        while d > math.pi:
            d -= 2.0 * math.pi
        total += d
        prev = cur
    w = total / (2.0 * math.pi)
    return int(round(w)), abs(w - round(w)), min_mod


def mean(vals):
    return sum(vals) / max(len(vals), 1)


def main():
    flush("=" * 78)
    flush("L19: canonical zero fingerprint")
    flush("=" * 78)
    if _shared.mp is None:
        flush("ERROR: mpmath required.")
        return

    bank = _shared.build_highk_bank()
    chi4 = next(ch for ch in _shared.primitive_characters() if ch["name"] == "chi4")
    chi_fn = chi4["fn"]

    sigma1, sigma2, t1, t2 = 0.3, 0.7, 0.0, 15.0
    fit_grid = [complex(s, t) for s in [1.20, 1.40, 1.70, 2.10] for t in [0.0, 1.0, 2.0]]

    def Lchi(s: complex) -> complex:
        return _shared.L_hurwitz(s, 4, chi_fn)

    mu_fit = [_shared.mu_chi_from_record_resolvent(bank[999], chi_fn, s) for s in fit_grid]
    l_fit = [Lchi(s) for s in fit_grid]
    c_ref = _shared.fit_c(mu_fit, l_fit)

    flush(f"Reference c_ref on Re(s)>1: {c_ref.real:+.6f}{c_ref.imag:+.6f}i")

    flush("\nIn-box windings for the canonical object")
    flush("-" * 78)
    for K in (19, 20, 21, 999):
        fK = lambda z, K=K: _shared.mu_chi_from_record_resolvent(bank[K], chi_fn, z)
        n0, drift, minbd = count_box_winding(fK, sigma1, sigma2, t1, t2, n_edge=320)
        flush(f"K={K:3d}: winding={n0:+d}  drift={drift:.3e}  min|bdry|={minbd:.3e}")

    fL = lambda z: Lchi(z)
    nL, driftL, minL = count_box_winding(fL, sigma1, sigma2, t1, t2, n_edge=320)
    flush(f"\nL(s,chi4): winding={nL:+d}  drift={driftL:.3e}  min|bdry|={minL:.3e}")

    fF = lambda z: _shared.mu_chi_from_record_resolvent(bank[999], chi_fn, z) / (c_ref * Lchi(z)) if abs(c_ref * Lchi(z)) > 1e-12 else complex(1e12, 0.0)
    nF, driftF, minF = count_box_winding(fF, sigma1, sigma2, t1, t2, n_edge=320)
    flush(f"F(s)=M/(cL): index={nF:+d}  drift={driftF:.3e}  min|bdry|={minF:.3e}")

    zero_heights = [6.020949, 10.243770, 12.988098]
    flush("\nCritical-line minima near L-zero heights")
    flush("-" * 78)
    for t0 in zero_heights:
        ts = [t0 + (j - 20) * 0.01 for j in range(41)]
        vals = [
            abs(_shared.mu_chi_from_record_resolvent(bank[999], chi_fn, complex(0.5, t)))
            for t in ts
        ]
        idx = min(range(len(vals)), key=lambda i: vals[i])
        flush(f"target t={t0:.6f}  min at t={ts[idx]:.6f}  |M|={vals[idx]:.6e}")

    flush("\n" + "=" * 78)
    flush("L19 VERDICT")
    flush("=" * 78)
    flush(
        "The canonical weighted resolvent object can now be tested directly for\n"
        "critical-strip zero structure. If the in-box winding remains zero while\n"
        "L(s,chi4) has nonzero winding, the zero-transfer obstruction persists\n"
        "even after passing to the canonical operator realization."
    )
    flush("=" * 78)


if __name__ == "__main__":
    main()
