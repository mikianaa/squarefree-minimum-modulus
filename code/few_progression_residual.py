"""Pigeonhole test for the extension L -> L q (least modulus 3, squarefree).

The tau(L) - 1 shapes d > 1 are spread over the q - 1 non-free slices, so some slice holds at most
k = floor((tau(L) - 1) / (q - 1)) shapes, and the q-free residual R satisfies R subset of a union of
at most k progressions a_i + d_i Z (d_i | L, d_i > 1, distinct).  This script asks CP-SAT for a
q-free family over Z_L (normal form of cpsat_cover) whose residual lies in such a union, with
|R| >= rlb.  INFEASIBLE proves that no cover with lcm L q exists.  When a family is found, the exact
slice assignment for its residual is tested (slice_extend.assign_slices); if that fails, that exact
residual R is forbidden (a smaller residual inside the same union could still be assignable, so the
union itself must NOT be forbidden) and the search continues (lazy loop; INFEASIBLE after valid cuts
is still a proof).
Usage: python few_progression_residual.py L q k rlb [--time 600] [--workers 8] [--iters 50] [--output out.json]
"""
import argparse, json, time
from ortools.sat.python import cp_model
from cpsat_cover import divisors
from slice_extend import assign_slices


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("L", type=int); ap.add_argument("q", type=int); ap.add_argument("k", type=int); ap.add_argument("rlb", type=int)
    ap.add_argument("--time", type=float, default=600); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--iters", type=int, default=50); ap.add_argument("--output")
    a = ap.parse_args()
    L, q, k = a.L, a.q, a.k
    D = [d for d in divisors(L) if d >= 3]; shapes = [d for d in divisors(L) if d > 1]
    m = cp_model.CpModel()
    x = {(d, r): m.NewBoolVar(f"x{d}_{r}") for d in D for r in range(d)}
    u = [m.NewBoolVar(f"u{n}") for n in range(L)]
    for d in D: m.AddExactlyOne(x[(d, r)] for r in range(d))
    for n in range(L):
        m.AddBoolOr([x[(d, n % d)] for d in D] + [u[n]])
        for d in D: m.AddImplication(x[(d, n % d)], u[n].Not())
    odd_primes = [p for p in range(3, L + 1) if L % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    for p in odd_primes:
        m.Add(x[(p, 0)] == 1)
        for r in range(2 * p):
            if r % p != 1: m.Add(x[(2 * p, r)] == 0)
    m.Add(x[(2 * odd_primes[0], 1)] == 1)
    m.Add(sum(u) >= a.rlb)
    z = [{(d, b): m.NewBoolVar(f"z{i}_{d}_{b}") for d in shapes for b in range(d)} for i in range(k)]
    for i in range(k): m.AddExactlyOne(z[i].values())
    # distinct shapes, ordered to break symmetry
    sh = [m.NewIntVar(2, L, f"sh{i}") for i in range(k)]
    for i in range(k):
        m.Add(sh[i] == sum(d * z[i][(d, b)] for (d, b) in z[i]))
    for i in range(k - 1): m.Add(sh[i] < sh[i + 1])
    for n in range(L):
        m.Add(sum(z[i][(d, n % d)] for i in range(k) for d in shapes) >= 1).OnlyEnforceIf(u[n])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.time; s.parameters.num_workers = a.workers
    log = []
    for it in range(a.iters):
        t = time.time(); st = s.Solve(m); name = s.StatusName(st)
        entry = dict(iteration=it, status=name, seconds=round(time.time() - t, 1))
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            union = [next((d, b) for (d, b), v in z[i].items() if s.Value(v)) for i in range(k)]
            fam = sorted((d, r) for (d, r), v in x.items() if s.Value(v))
            R = [n for n in range(L) if s.Value(u[n])]
            assert all(not any(n % d == r for d, r in fam) for n in R)
            assert all(any(n % d == b for d, b in union) for n in R)
            asg, st2 = assign_slices(L, R, q, time_limit=300, workers=a.workers)
            entry.update(union=union, uncovered=len(R), slice_assignment=st2)
            print(json.dumps(entry), flush=True)
            if asg is not None:
                entry["family"] = fam; entry["R"] = R; entry["assignment"] = {d: list(v) for d, v in asg.items()}
                log.append(entry); print("COVER at L*q EXISTS: family + assignment recorded"); break
            # forbid exactly this residual (valid cut); forbidding the union would be unsound
            m.AddBoolOr([u[n].Not() for n in R] + [u[n] for n in range(L) if n not in set(R)])
        else:
            entry["verdict"] = "NO COVER at L*q (proved)" if st == cp_model.INFEASIBLE else "undecided (time limit)"
            print(json.dumps(entry), flush=True)
        log.append(entry)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
    if a.output: json.dump(log, open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
