"""Slice-extension probe.  For squarefree L and a divisor m >= 3 of L, decide whether Z_L can be covered
by one progression for every divisor d >= 3 of L plus a SECOND progression with modulus m (CP-SAT).
If so the residual of the distinct family (moduli d >= 3) can be made a single class mod m, which is
the most concentrated residual with L/m points; the slice argument then asks whether the shapes
q*d (d | L) can cover that class in each of the q-1 slices of Z_{Lq}.
Only translation normalisation is used (class of each odd prime at 0, parity of 2*p_min fixed).
Usage: python duplicate_test.py --L 2310 [--time 60] [--workers 8] [--m 105]
"""
import argparse, json, time
from ortools.sat.python import cp_model
from cpsat_cover import divisors


def solve(L, m_dup, time_limit, workers):
    D = [d for d in divisors(L) if d >= 3]
    copies = [(d, i) for d in D for i in range(2 if d == m_dup else 1)]
    mdl = cp_model.CpModel()
    x = {(d, i, r): mdl.NewBoolVar(f"x{d}_{i}_{r}") for d, i in copies for r in range(d)}
    for d, i in copies:
        mdl.AddExactlyOne(x[(d, i, r)] for r in range(d))
    for n in range(L):
        mdl.AddBoolOr([x[(d, i, n % d)] for d, i in copies])
    odd_primes = [q for q in range(3, L + 1) if L % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    for q in odd_primes:
        mdl.Add(x[(q, 0, 0)] == 1)
    q0 = odd_primes[0]
    for r in range(2 * q0):
        if r % 2 == 0: mdl.Add(x[(2 * q0, 0, r)] == 0)
    if m_dup in D:   # break the symmetry between the two copies
        mdl.Add(sum(r * x[(m_dup, 0, r)] for r in range(m_dup)) <= sum(r * x[(m_dup, 1, r)] for r in range(m_dup)))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    t = time.time(); st = s.Solve(mdl)
    res = dict(L=L, m_dup=m_dup, status=s.StatusName(st), seconds=round(time.time() - t, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sysm = sorted((d, i, r) for (d, i, r), v in x.items() if s.Value(v))
        assert all(any(n % d == r for d, i, r in sysm) for n in range(L))
        res["status"] = "SAT"; res["system"] = [dict(modulus=d, copy=i, residue=r) for d, i, r in sysm]
    elif st == cp_model.INFEASIBLE:
        res["status"] = "UNSAT"
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=2310); ap.add_argument("--time", type=float, default=60)
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--m", type=int); ap.add_argument("--output")
    a = ap.parse_args()
    ms = [a.m] if a.m else [d for d in divisors(a.L) if d >= 3]
    out = []
    for m in ms:
        r = solve(a.L, m, a.time, a.workers)
        print(json.dumps({k: v for k, v in r.items() if k != "system"}), flush=True)
        out.append(r)
    if a.output: json.dump(out, open(a.output, "w"), indent=1)
