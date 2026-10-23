"""Numerical ceiling of the truncation idea for the squarefree comparison (NOT a certificate).

For a schedule (C0, bands, deltas) compute, in floating point, the comparison stage bounds
    untruncated:  E (S - delta)_+ / (1 - delta),            S = 2^N / p
    truncated:    E (min(S, 1) - delta)_+ / (1 - delta)
for every stage p > C0 up to P_EXACT, where N is the Poisson-binomial count of present earlier
primes (prime q present with probability nu_q / q, nu_q = 1/(1 - delta_q)).  Stages p <= C0 are
taken from a verifier output (their S never exceeds 1 at C0 = 8).  If even the truncated total
is >= 1 the truncation direction is dead for this C0.
Usage: python truncation_ceiling.py schedule.json verifier_out.json [P_EXACT]
"""
import json, sys, bisect
import numpy as np


def primes_to(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]


def main():
    w = json.load(open(sys.argv[1])); v = json.load(open(sys.argv[2]))
    P = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    C0 = w["C0"]; edges = w["edges"]; params = w["params"]
    def delta(k):  # k = 1-based stage index
        i = bisect.bisect_left(edges, k); return params[min(i, len(params) - 1)]
    ps = primes_to(P)
    rows = {r["p"]: r for r in (v.get("rows") or [])}
    NMAX = 120
    dist = np.zeros(NMAX + 1); dist[0] = 1.0   # law of N (count of present earlier primes)
    tot_u = tot_t = 0.0; small = 0.0
    per = []
    for k, p in enumerate(ps, start=1):
        d = delta(k)
        if p <= C0:
            loss = float(rows[p]["loss"]) if p in rows else 0.0
            small += loss; per.append((p, loss, loss))
        else:
            n = np.arange(NMAX + 1)
            S = 2.0 ** n / p
            u = np.sum(dist * np.maximum(S - d, 0)) / (1 - d)
            t = np.sum(dist * np.maximum(np.minimum(S, 1) - d, 0)) / (1 - d)
            tot_u += u; tot_t += t; per.append((p, u, t))
        pr = (1 / (1 - d)) / p
        dist = dist * (1 - pr) + np.concatenate(([0.0], dist[:-1])) * pr
    print(f"C0={C0} P_EXACT={P}: stages p<=C0 {small:.6f}; p>C0: untruncated {tot_u:.6f}, truncated {tot_t:.6f}, saving {tot_u - tot_t:.6f}")
    print("far+tail from verifier:", v["parts"].get("far"), v["parts"].get("tail"))
    print("per stage (p, untruncated, truncated, saving) for p < 60:")
    for p, u, t in per:
        if p < 60: print(f"  {p:3d} {u:.6f} {t:.6f} {u - t:.6f}")
    far = float(v["parts"]["far"]) + float(v["parts"]["tail"])
    print(f"TOTAL untruncated {small + tot_u + far:.6f}   truncated ceiling {small + tot_t + far:.6f}")


if __name__ == "__main__":
    main()
