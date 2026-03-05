"""
Shared utilities for Paper L experiments.

Contains data-loading, profile-building, character definitions,
Dirichlet L-function evaluation, and D-odd enumeration routines.
All data files are loaded from the local data/ directory.
"""

import cmath
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

try:
    import mpmath as mp
except ImportError:
    mp = None

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# ── E45 parsing and profile building ────────────────────────────

def parse_e45_file(path: str):
    results = []
    current = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            m = re.match(r">>> COMPLETED K=(\d+)", line)
            if m:
                current = {
                    "K": int(m.group(1)),
                    "cas_00": {},
                    "cas_10": {},
                    "stop_00": {},
                    "stop_10": {},
                }
                results.append(current)
                continue
            if current is None:
                continue
            m = re.match(r"\s+D\s+=\s+(\d+)", line)
            if m:
                current["D"] = int(m.group(1))
            m = re.match(r"\s+S00\s+=\s+(-?\d+)\s+\(n00\s*=\s*(\d+)\)", line)
            if m:
                current["S00"] = int(m.group(1))
                current["n00"] = int(m.group(2))
            m = re.match(r"\s+S10\s+=\s+(-?\d+)\s+\(n10\s*=\s*(\d+)\)", line)
            if m:
                current["S10"] = int(m.group(1))
                current["n10"] = int(m.group(2))
            m = re.match(
                r"\s+(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$", line
            )
            if m:
                d = int(m.group(1))
                current["cas_00"][d] = int(m.group(2))
                current["cas_10"][d] = int(m.group(3))
                current["stop_00"][d] = int(m.group(4))
                current["stop_10"][d] = int(m.group(5))
    return results


def load_e45():
    """Load E45 high-K profiles from local data/ directory."""
    out = []
    for fname in ["E45_K19_K20.txt", "E45_K21.txt"]:
        fpath = os.path.join(_DATA_DIR, fname)
        if os.path.exists(fpath):
            out.extend(parse_e45_file(fpath))
    if not out:
        raise FileNotFoundError(
            f"E45 data not found in {_DATA_DIR}. "
            "Ensure data/E45_K19_K20.txt and data/E45_K21.txt are present."
        )
    return sorted(out, key=lambda x: x["K"])


def load_e162_R_history() -> Dict[int, float]:
    """Load R(K) = sigma_10/sigma_00 from local data/ directory."""
    fpath = os.path.join(_DATA_DIR, "E162_K19_analysis_output.txt")
    R_hist: Dict[int, float] = {}
    if not os.path.exists(fpath):
        return R_hist
    with open(fpath) as f:
        for line in f:
            m = re.search(r"K\s*=\s*(\d+).*?=\s*([+-]?\d+\.\d+)\s+Q/Nt", line)
            if m:
                R_hist[int(m.group(1))] = float(m.group(2))
    return R_hist


# ── TauProfile ──────────────────────────────────────────────────

@dataclass
class TauProfile:
    K: int
    sec: str
    taus: List[int]
    p: Dict[int, float]
    e: Dict[int, float]
    u: Dict[int, float]
    q: Dict[int, float]
    a: Dict[int, float]
    ccdf: Dict[int, float]
    mu_direct: float


def build_profile(res: dict, sec: str) -> TauProfile:
    if sec not in ("00", "10"):
        raise ValueError("sector must be 00 or 10")

    D = res["D"]
    n = res["n00"] if sec == "00" else res["n10"]
    stop_map = res["stop_00"] if sec == "00" else res["stop_10"]
    cas_map = res["cas_00"] if sec == "00" else res["cas_10"]

    p: Dict[int, float] = {}
    e: Dict[int, float] = {}
    u: Dict[int, float] = {}
    taus: List[int] = []
    for d, s in stop_map.items():
        tau = D - d
        if s <= 0:
            continue
        taus.append(tau)
        p_tau = s / n
        e_tau = cas_map[d] / s
        p[tau] = p_tau
        e[tau] = e_tau
        u[tau] = p_tau * e_tau

    taus = sorted(taus)
    tau_min = taus[0]
    tau_max = taus[-1]

    ccdf: Dict[int, float] = {}
    cdf = 0.0
    ccdf[tau_min - 1] = 1.0
    for t in range(tau_min, tau_max + 1):
        cdf += p.get(t, 0.0)
        ccdf[t] = max(0.0, 1.0 - cdf)

    q: Dict[int, float] = {}
    a: Dict[int, float] = {}
    for t in taus:
        s_prev = ccdf[t - 1]
        if s_prev <= 1e-18:
            q[t] = 0.0
            a[t] = 0.0
            continue
        q[t] = min(max(ccdf[t] / s_prev, 0.0), 1.0)
        a[t] = p[t] / s_prev

    mu_direct = sum(u.values())
    return TauProfile(
        K=res["K"], sec=sec, taus=taus,
        p=p, e=e, u=u, q=q, a=a, ccdf=ccdf,
        mu_direct=mu_direct,
    )


def richardson_map(v20: Dict[int, float], v21: Dict[int, float]) -> Dict[int, float]:
    out: Dict[int, float] = {}
    taus = sorted(set(v20.keys()) | set(v21.keys()))
    for t in taus:
        if t in v20 and t in v21:
            out[t] = 2.0 * v21[t] - v20[t]
        elif t in v21:
            out[t] = v21[t]
    return out


def normalize_prob(p_map: Dict[int, float]) -> Dict[int, float]:
    cleaned = {t: max(0.0, p) for t, p in p_map.items()}
    s = sum(cleaned.values())
    if s <= 1e-18:
        return cleaned
    return {t: p / s for t, p in cleaned.items()}


def build_extrapolated_profile(
    prof20: TauProfile, prof21: TauProfile, sec: str
) -> TauProfile:
    p_inf = normalize_prob(richardson_map(prof20.p, prof21.p))
    e_inf = richardson_map(prof20.e, prof21.e)

    taus = sorted(t for t in p_inf.keys() if p_inf[t] > 0.0)
    if not taus:
        raise RuntimeError("empty extrapolated profile")

    u_inf: Dict[int, float] = {}
    for t in taus:
        if t not in e_inf:
            e_inf[t] = prof21.e.get(t, 0.0)
        u_inf[t] = p_inf[t] * e_inf[t]

    tau_min = taus[0]
    tau_max = taus[-1]
    ccdf: Dict[int, float] = {}
    cdf = 0.0
    ccdf[tau_min - 1] = 1.0
    for t in range(tau_min, tau_max + 1):
        cdf += p_inf.get(t, 0.0)
        ccdf[t] = max(0.0, 1.0 - cdf)

    q: Dict[int, float] = {}
    a: Dict[int, float] = {}
    for t in taus:
        s_prev = ccdf[t - 1]
        if s_prev <= 1e-18:
            q[t] = 0.0
            a[t] = 0.0
        else:
            q[t] = min(max(ccdf[t] / s_prev, 0.0), 1.0)
            a[t] = p_inf[t] / s_prev

    return TauProfile(
        K=999, sec=sec, taus=taus,
        p=p_inf, e=e_inf, u=u_inf, q=q, a=a, ccdf=ccdf,
        mu_direct=sum(u_inf.values()),
    )


# ── High-K bank builder ────────────────────────────────────────

def build_highk_bank() -> dict:
    """Load E45, build profiles for K=19,20,21 + Richardson K=999."""
    all_res = load_e45()
    by_k = {r["K"]: r for r in all_res}

    prof = {}
    for k in (19, 20, 21):
        prof[(k, "00")] = build_profile(by_k[k], "00")
        prof[(k, "10")] = build_profile(by_k[k], "10")
    prof[(999, "00")] = build_extrapolated_profile(
        prof[(20, "00")], prof[(21, "00")], "00"
    )
    prof[(999, "10")] = build_extrapolated_profile(
        prof[(20, "10")], prof[(21, "10")], "10"
    )

    omega = {k: by_k[k]["n10"] / by_k[k]["n00"] for k in (19, 20, 21)}
    omega[999] = 2.0 * omega[21] - omega[20]

    bank = {}
    for k in (19, 20, 21, 999):
        bank[k] = {
            "omega": omega[k],
            "u00": prof[(k, "00")].u,
            "u10": prof[(k, "10")].u,
        }
    return bank


# ── u_mix and Dirichlet series ──────────────────────────────────

def u_mix_map(rec: dict) -> Dict[int, float]:
    u00 = rec["u00"]
    u10 = rec["u10"]
    omega = rec["omega"]
    taus = sorted(set(u00.keys()) | set(u10.keys()))
    return {t: omega * u10.get(t, 0.0) - u00.get(t, 0.0) for t in taus}


def n_of_tau(tau: int, tau0: int = 3) -> int:
    return 2 * (tau - tau0) + 1


def mu_chi(
    u_mix: Dict[int, float],
    chi_fn: Callable[[int], int],
    s: complex,
    tau0: int = 3,
) -> complex:
    total = 0.0 + 0.0j
    for tau in sorted(u_mix.keys()):
        if tau < tau0:
            continue
        n = n_of_tau(tau, tau0=tau0)
        c = chi_fn(n)
        if c == 0:
            continue
        total += (u_mix[tau] * c) * cmath.exp(-s * math.log(n))
    return total


# ── Characters ──────────────────────────────────────────────────

def primitive_characters():
    def chi3(n: int) -> int:
        r = n % 3
        if r == 0:
            return 0
        return 1 if r == 1 else -1

    def chi4(n: int) -> int:
        r = n % 4
        if r % 2 == 0:
            return 0
        return 1 if r == 1 else -1

    def chi5(n: int) -> int:
        r = n % 5
        if r == 0:
            return 0
        return 1 if r in (1, 4) else -1

    def chi8(n: int) -> int:
        r = n % 8
        if r % 2 == 0:
            return 0
        return 1 if r in (1, 7) else -1

    return [
        {"name": "chi3", "q": 3, "fn": chi3},
        {"name": "chi4", "q": 4, "fn": chi4},
        {"name": "chi5", "q": 5, "fn": chi5},
        {"name": "chi8", "q": 8, "fn": chi8},
    ]


# ── L-function evaluation ──────────────────────────────────────

def L_hurwitz(s: complex, q: int, chi_fn: Callable[[int], int]) -> complex:
    """Evaluate L(s, chi) via Hurwitz zeta (requires mpmath)."""
    if mp is None:
        raise RuntimeError("mpmath is required for L_hurwitz but not installed.")
    mp.mp.dps = 70
    ss = mp.mpc(s.real, s.imag)
    total = mp.mpc(0)
    qpow = mp.power(q, ss)
    for a in range(1, q + 1):
        ca = chi_fn(a)
        if ca == 0:
            continue
        total += ca * mp.hurwitz(ss, mp.mpf(a) / q)
    return complex(total / qpow)


# ── Fitting utilities ───────────────────────────────────────────

def fit_c(mu_vals: List[complex], basis_vals: List[complex]) -> complex:
    num = 0.0 + 0.0j
    den = 0.0
    for m, b in zip(mu_vals, basis_vals):
        num += m * b.conjugate()
        den += b.real * b.real + b.imag * b.imag
    if den <= 1e-30:
        return complex(float("nan"), float("nan"))
    return num / den


def rel(a: complex, b: complex) -> float:
    return abs(a - b) / max(abs(a), 1e-15)


def mean(vals: List[float]) -> float:
    return sum(vals) / max(len(vals), 1)


def primes_upto(n: int) -> List[int]:
    out = []
    for x in range(2, n + 1):
        ok = True
        d = 2
        while d * d <= x:
            if x % d == 0:
                ok = False
                break
            d += 1
        if ok:
            out.append(x)
    return out


# ── D-odd enumeration ──────────────────────────────────────────

def enumerate_dodd_sector(K: int, sector: str):
    """Yield (x, y, carries[0..D]) for all D-odd pairs in given sector."""
    D = 2 * K - 1
    lo = 1 << (K - 1)
    hi = 1 << K

    for x in range(lo, hi):
        a1 = (x >> (K - 2)) & 1
        for y in range(lo, hi):
            prod = x * y
            if prod.bit_length() != D:
                continue
            b1 = (y >> (K - 2)) & 1
            sec_bits = a1 * 2 + b1
            if sector == "00" and sec_bits != 0:
                continue
            if sector == "10" and sec_bits != 2:
                continue

            carries = [0] * (D + 1)
            for d in range(D):
                cv = 0
                for i in range(max(0, d - K + 1), min(d, K - 1) + 1):
                    cv += ((x >> i) & 1) * ((y >> (d - i)) & 1)
                carries[d + 1] = (cv + carries[d]) >> 1
            yield x, y, carries
