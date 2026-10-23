"""Minimum number of uncovered residues modulo a squarefree N by distinct progressions with moduli
d | N, d >= M (CP-SAT optimisation).  Uses the lossless normal form of cpsat_cover.solve (all moduli
used, class of each odd prime p at 0, class of 2p at 1 mod p, parity of the smallest 2p fixed), which
is valid for the minimum as well: adding unused moduli cannot increase the uncovered set, and the
normalisation acts by bijections of Z_N.  An OPTIMAL status proves the lower bound on the uncovered
count; FEASIBLE gives only an upper bound (the reported configuration).
Usage: python cpsat_min_uncovered.py --N 2310 --M 3 [--time 600] [--workers 8] [--output out.json]
"""
import argparse, json, math, time
from ortools.sat.python import cp_model
from cpsat_cover import divisors


def solve(N, M, time_limit=600.0, workers=8):
    D = [d for d in divisors(N) if d >= M]
    t = time.time()
    m = cp_model.CpModel()
    x = {(d, r): m.NewBoolVar(f"x{d}_{r}") for d in D for r in range(d)}
    u = [m.NewBoolVar(f"u{n}") for n in range(N)]
    for d in D:
        m.AddExactlyOne(x[(d, r)] for r in range(d))
    for n in range(N):
        m.AddBoolOr([x[(d, n % d)] for d in D] + [u[n]])
    odd_primes = [q for q in range(3, N + 1) if N % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    assert all(N % (q * q) for q in range(2, N + 1) if N % q == 0 and q * q <= N), "N must be squarefree"
    for q in odd_primes:
        if q in D: m.Add(x[(q, 0)] == 1)
        if 2 * q in D:
            for r in range(2 * q):
                if r % q != 1: m.Add(x[(2 * q, r)] == 0)
    q0 = odd_primes[0]
    if 2 * q0 in D: m.Add(x[(2 * q0, 1)] == 1)
    m.Minimize(sum(u))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    st = s.Solve(m)
    res = dict(N=N, M=M, status=s.StatusName(st), seconds=round(time.time() - t, 1), moduli=len(D))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        res["uncovered"] = int(round(s.ObjectiveValue())); res["lower_bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-9))
        sysm = sorted((d, r) for (d, r), v in x.items() if s.Value(v))
        unc = [n for n in range(N) if not any(n % d == r for d, r in sysm)]
        assert len(unc) == res["uncovered"], (len(unc), res["uncovered"])
        res["system"] = [dict(modulus=d, residue=r) for d, r in sysm]; res["uncovered_residues"] = unc
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True); ap.add_argument("--M", type=int, default=3)
    ap.add_argument("--time", type=float, default=600.0); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output")
    a = ap.parse_args()
    r = solve(a.N, a.M, a.time, a.workers)
    print(json.dumps({k: v for k, v in r.items() if k not in ("system", "uncovered_residues")}))
    if "uncovered_residues" in r: print("uncovered:", r["uncovered_residues"][:60])
    if a.output: json.dump(r, open(a.output, "w"), indent=1)
