"""Floating-point optimiser for the squarefree convex stop-loss schedule (NOT a certificate).

Model (squarefree, cutoff C0: all moduli > C0).  Stage k, prime p = p_k, earlier primes q present
independently with probability nu_q/q, nu_q = 1/(1-delta_q).  S = (1/p) #{m | D : m > C0/p}.
Stage bound = min( E S , E(S-d)_+/(1-d) , E S^2/(4d(1-d)) ).
Exact range: primes <= P_EXACT (Poisson-binomial on the count N of present primes, truncated at NMAX;
for p <= C0 condition on the primes <= C0/p).  Stages above P_EXACT up to stage N_STAGES: closed
forms E S = C/p, E S^2 = T/p^2.  Tail beyond N_STAGES: CFT Lemma 4.7 with delta = 1/2.
Output: JSON with the schedule (band upper stage indices + deltas) and the totals.
Usage: python sqfree_float_opt.py C0 out.json [P_EXACT] [sweeps] [init.json]
"""
import json, math, sys, time
import numpy as np

N_STAGES = 10 ** 6
NMAX = 120


def primes_to(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


PR = primes_to(16_000_000)
assert len(PR) >= N_STAGES
POW2 = 2.0 ** np.arange(NMAX + 1)


def divisor_counts_le(small, b):
    n = len(small); out = [0] * (1 << n)
    for mask in range(1 << n):
        divs = [1]
        for i in range(n):
            if mask >> i & 1: divs = divs + [d * small[i] for d in divs]
        out[mask] = sum(1 for d in divs if d <= b)
    return out


def pb_dist(prs):
    dist = np.zeros(NMAX + 1); dist[0] = 1.0
    for pr in prs:
        new = dist * (1 - pr); new[1:] += dist[:-1] * pr; dist = new
    return dist


def stage_bounds(S, w, d):
    ES = float(np.dot(w, S)); ES2 = float(np.dot(w, S * S))
    conv = float(np.dot(w, np.maximum(S - d, 0.0))) / (1 - d)
    quad = ES2 / (4 * d * (1 - d)) if d > 0 else math.inf
    return ES, ES2, conv, quad


def tail_cft(pN, M0):
    """CFT Lemma 4.7: sum_{k>N} w_k(B_k) <= 2 c2 M0 e^{6 c1} poly(log pN)/pN, delta=1/2 beyond N."""
    L = math.log(pN)
    c1 = -math.log(L) + 1 / L ** 2
    c2 = 1 + 3 / (2 * L)
    poly = L ** 5 + 5 * L ** 4 + 20 * L ** 3 + 60 * L ** 2 + 120 * L + 120
    return 2 * c2 * M0 * math.exp(6 * c1) * poly / pN


class Model:
    def __init__(self, C0, P_EXACT):
        self.C0 = C0
        self.ps_exact = [int(x) for x in PR[PR <= P_EXACT]]
        self.n_exact = len(self.ps_exact)
        self.ps_far = PR[self.n_exact:N_STAGES].astype(float)

    def total(self, delta_of_k, want_rows=False):
        C0 = self.C0; nu = {}; rows = []; tot = 0.0
        dist = None  # Poisson-binomial over all earlier primes (used for p > C0)
        C = 1.0; T = 1.0; E3 = 1.0
        for k, p in enumerate(self.ps_exact, start=1):
            d = delta_of_k(k)
            b = C0 / p
            if p > C0:
                if dist is None:
                    dist = pb_dist([nu[q] / q for q in nu if q < p])
                S = POW2 / p; w = dist
                ES, ES2, conv, quad = stage_bounds(S, w, d)
                # truncation remainder (upper bound) E[2^N 1{N>NMAX}]/p <= (2/3)^NMAX E 3^N / p
                rem = (2 / 3) ** NMAX * E3 / p
                conv += rem / (1 - d); ES += rem
            else:
                qs = [q for q in nu if q < p]
                small = [q for q in qs if q <= b]; large = [q for q in qs if q > b]
                dl = pb_dist([nu[q] / q for q in large])
                cnt = divisor_counts_le(small, b)
                ES = ES2 = conv = quad = 0.0
                for mask in range(1 << len(small)):
                    pm = 1.0; kk = 0
                    for i, q in enumerate(small):
                        pr = nu[q] / q
                        if mask >> i & 1: pm *= pr; kk += 1
                        else: pm *= 1 - pr
                    S = (POW2 * 2 ** kk - cnt[mask]) / p
                    a, b2, c, q2 = stage_bounds(S, dl, d)
                    ES += pm * a; ES2 += pm * b2; conv += pm * c
                quad = ES2 / (4 * d * (1 - d)) if d > 0 else math.inf
            loss = min(ES, conv, quad)
            tot += loss
            if want_rows: rows.append(dict(k=k, p=p, delta=d, ES=ES, ES2=ES2, convex=conv, quad=quad, loss=loss))
            nu[p] = 1 / (1 - d); pr = nu[p] / p
            C *= 1 + pr; T *= 1 + 3 * pr; E3 *= 1 + 2 * pr
            if dist is not None:
                new = dist * (1 - pr); new[1:] += dist[:-1] * pr; dist = new
        # far stages: closed forms with band deltas
        ks = np.arange(self.n_exact + 1, N_STAGES + 1)
        ds = delta_of_k.vec(ks)
        nus = 1 / (1 - ds); ps = self.ps_far
        fT = np.concatenate(([T], T * np.cumprod(1 + 3 * nus[:-1] / ps[:-1])))
        fC = np.concatenate(([C], C * np.cumprod(1 + nus[:-1] / ps[:-1])))
        quad = fT / ps ** 2 / (4 * ds * (1 - ds)); first = fC / ps
        far = np.minimum(quad, first); far_sum = float(far.sum())
        M0 = fT[-1] * (1 + 3 * nus[-1] / ps[-1])
        tail = tail_cft(float(ps[-1]), M0)
        return tot, far_sum, tail, rows


def make_bands(per_stage=30):
    """Band upper stage indices.  Per-stage up to k=per_stage, then growing bands."""
    edges = list(range(1, per_stage + 1))
    for e in [35, 40, 50, 60, 75, 100, 130, 170, 230, 300, 400, 500, 700, 1000, 1400, 2000, 3000, 5000, 7000, 10000,
              15000, 20000, 30000, 50000, 80000, 100000, 150000, 200000, 300000, 500000, 700000, 1000000]:
        if e > edges[-1]: edges.append(e)
    return edges


class banded_fn:
    def __init__(self, params, edges):
        self.arr = np.array(edges); self.params = np.array(params, dtype=float)
    def __call__(self, k):
        return float(self.params[int(np.searchsorted(self.arr, k))])
    def vec(self, ks):
        return self.params[np.searchsorted(self.arr, ks)]


def main():
    C0 = int(sys.argv[1]); out = sys.argv[2]
    P_EXACT = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 5
    sweeps = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    init_file = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != "-" else None
    per_stage = int(sys.argv[6]) if len(sys.argv) > 6 else 30
    t0 = time.time()
    model = Model(C0, P_EXACT)
    edges = make_bands(per_stage)
    nb = len(edges)
    params = []
    for i, e in enumerate(edges):
        p = int(PR[e - 1])
        if p <= 13: params.append(0.0)
        elif p < 30: params.append(0.2)
        elif p < 100: params.append(0.27)
        elif p < 1000: params.append(0.31)
        elif p < 10000: params.append(0.35)
        else: params.append(0.40)
    if init_file:
        prev = json.load(open(init_file))
        pf = banded_fn(prev["params"], prev["edges"])
        params = [pf(e) for e in edges]

    def ev(ps):
        a, b, c, _ = model.total(banded_fn(ps, edges)); return a + b + c
    cur = ev(params)
    print(f"C0={C0} P_EXACT={P_EXACT} n_exact={model.n_exact} bands={nb} initial total={cur:.6f} ({time.time()-t0:.1f}s)", flush=True)
    step = 0.06
    for sweep in range(sweeps):
        improved = False
        for i in range(nb):
            for sgn in (+1, -1):
                trial = list(params); trial[i] = min(0.5, max(0.0, trial[i] + sgn * step))
                if trial[i] == params[i]: continue
                v = ev(trial)
                if v < cur - 1e-10:
                    params, cur = trial, v; improved = True
                    while True:
                        t2 = list(params); t2[i] = min(0.5, max(0.0, t2[i] + sgn * step))
                        if t2[i] == params[i]: break
                        v2 = ev(t2)
                        if v2 < cur - 1e-10: params, cur = t2, v2
                        else: break
                    break
        print(f"  sweep {sweep} step {step:.4f} total={cur:.6f} ({time.time()-t0:.0f}s)", flush=True)
        if not improved: step /= 2
        a, b, c, rows = model.total(banded_fn(params, edges), want_rows=True)
        json.dump(dict(C0=C0, P_EXACT=P_EXACT, n_exact=model.n_exact, edges=edges, params=params,
                       exact_sum=a, far_sum=b, tail=c, total=a + b + c,
                       rows_first_60=rows[:60]), open(out, "w"), indent=1)
    a, b, c, rows = model.total(banded_fn(params, edges), want_rows=True)
    json.dump(dict(C0=C0, P_EXACT=P_EXACT, n_exact=model.n_exact, edges=edges, params=params,
                   exact_sum=a, far_sum=b, tail=c, total=a + b + c,
                   rows_first_60=rows[:60]), open(out, "w"), indent=1)
    print(f"FINAL C0={C0}: exact={a:.6f} far={b:.6f} tail={c:.6f} total={a+b+c:.6f}")
    for r in rows[:25]:
        print(f"   k={r['k']} p={r['p']} d={r['delta']:.3f} ES={r['ES']:.5f} conv={r['convex']:.5f} quad={r['quad']:.5f} loss={r['loss']:.5f}")


if __name__ == "__main__":
    main()
