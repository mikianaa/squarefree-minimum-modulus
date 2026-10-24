"""Slice-capacity relaxation for the extension L -> L q (least modulus 3, squarefree).

Over all families F over Z_L (one progression per divisor d >= 3 of L; the normal form of
cpsat_cover is used) with uncovered set R, maximise
      cap(F) = Sum_{d | L, d > 1} c_d(R) - (q - 1) |R|,     c_d(R) = max_a |R cap (a + dZ)|.
A cover of Z_{Lq} needs a slice assignment, which needs cap(F) >= 0 for its q-free part (each of
the q-1 slices not covered by the modulus q must cover R by progressions mod d, one per shape d).
So an OPTIMAL value < 0 proves that no cover with lcm Lq and least modulus 3 exists.
Linearisation: c_d <= Sum_{n = b mod d} u_n + L (1 - z_{d,b}), Sum_b z_{d,b} = 1.
Usage: python slice_capacity.py L q [--time 1200] [--workers 8] [--output out.json]
"""
import argparse, json, time
from ortools.sat.python import cp_model
from cpsat_cover import divisors


def solve(L, q, time_limit, workers):
    D = [d for d in divisors(L) if d >= 3]
    shapes = [d for d in divisors(L) if d > 1]
    m = cp_model.CpModel()
    x = {(d, r): m.NewBoolVar(f"x{d}_{r}") for d in D for r in range(d)}
    u = [m.NewBoolVar(f"u{n}") for n in range(L)]
    for d in D: m.AddExactlyOne(x[(d, r)] for r in range(d))
    for n in range(L):
        m.AddBoolOr([x[(d, n % d)] for d in D] + [u[n]])
        for d in D: m.AddImplication(x[(d, n % d)], u[n].Not())      # u_n = 1 only if uncovered
    odd_primes = [p for p in range(3, L + 1) if L % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    for p in odd_primes:
        m.Add(x[(p, 0)] == 1)
        for r in range(2 * p):
            if r % p != 1: m.Add(x[(2 * p, r)] == 0)
    m.Add(x[(2 * odd_primes[0], 1)] == 1)
    R = sum(u)
    c = {}
    for d in shapes:
        c[d] = m.NewIntVar(0, L // d, f"c{d}")
        z = [m.NewBoolVar(f"z{d}_{b}") for b in range(d)]
        m.AddExactlyOne(z)
        for b in range(d):
            m.Add(c[d] <= sum(u[n] for n in range(b, L, d)) + (L // d) * (1 - z[b]))
    m.Maximize(sum(c.values()) - (q - 1) * R)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    t = time.time(); st = s.Solve(m)
    res = dict(L=L, q=q, status=s.StatusName(st), seconds=round(time.time() - t, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        res["objective"] = int(round(s.ObjectiveValue())); res["bound"] = s.BestObjectiveBound()
        res["uncovered"] = int(round(s.Value(R))); res["c_d"] = {d: s.Value(c[d]) for d in shapes}
        res["system"] = [dict(modulus=d, residue=r) for (d, r), v in x.items() if s.Value(v)]
        # BestObjectiveBound is a valid upper bound on the maximum: negative => no family reaches cap >= 0
        res["verdict"] = ("NO COVER at L*q (proved: bound < 0)" if res["bound"] < 0 else
                          "not excluded (relaxation feasible)" if res["objective"] >= 0 else "undecided")
    else:
        res["bound"] = s.BestObjectiveBound(); res["verdict"] = "NO COVER at L*q (proved: bound < 0)" if res["bound"] < 0 else "undecided"
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("L", type=int); ap.add_argument("q", type=int)
    ap.add_argument("--time", type=float, default=1200); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output")
    a = ap.parse_args()
    r = solve(a.L, a.q, a.time, a.workers)
    print(json.dumps({k: v for k, v in r.items() if k not in ("system", "c_d")}))
    if "c_d" in r: print("c_d:", r["c_d"])
    if a.output: json.dump(r, open(a.output, "w"), indent=1)
