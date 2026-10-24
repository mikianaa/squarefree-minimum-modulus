"""CP-SAT (OR-Tools) search for distinct covering systems with moduli dividing N
and minimum modulus >= M.  Same exact model as code/fable_sat_cover.py:
one Boolean per (modulus d | N, d >= M, residue r), at most one residue per
modulus, every residue class mod N covered.  Symmetry: modulus N, if used,
has residue 0 (translation).  UNSAT (status INFEASIBLE) is a proof of nonexistence.

Usage: python code/fable_cpsat_cover.py --M 4 --N 360 [--time 600] [--workers 8]
Originally code/fable_cpsat_cover.py of the minimum-modulus repository.
"""
import argparse, json, math, time
from fractions import Fraction
from ortools.sat.python import cp_model


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n: ds.append(n // i)
        i += 1
    return sorted(ds)


def solve(N, M, time_limit=600.0, workers=8, require=(), forbid=(), redundant=False, linearization=1, normalize=False, slice_cuts=0):
    """normalize=True (squarefree N, M >= 3): lossless normal form for existence.  Any cover extends
    to one using every divisor >= M of N (add unused moduli anywhere), so all moduli are required;
    translating each CRT coordinate Z_p puts the class of modulus p at residue 0 for every odd prime
    p | N; a class 2p whose residue is 0 mod p is redundant and can be moved, and the units of Z_p
    then map its residue mod p to 1 (parity free); the free translation of Z_2 (modulus 2 absent)
    fixes the parity of the class 2p for the smallest odd prime to 1.  The modulus-N translation
    symmetry constraint is dropped in this mode."""
    D = [d for d in divisors(N) if d >= M and d not in forbid]
    if sum(Fraction(1, d) for d in D) < 1:
        return dict(N=N, M=M, status="INFEASIBLE", reason="sum of 1/d below 1", seconds=0.0)
    t = time.time()
    m = cp_model.CpModel()
    x = {(d, r): m.NewBoolVar(f"x{d}_{r}") for d in D for r in range(d)}
    for d in D:
        m.AddAtMostOne(x[(d, r)] for r in range(d))
    for n in range(N):
        m.AddBoolOr([x[(d, n % d)] for d in D])
    if normalize:
        assert M >= 3 and all(N % (q * q) for q in range(2, 200) if N % q == 0), "normalize needs squarefree N and M >= 3"
        odd_primes = [q for q in range(3, N + 1) if N % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
        for d in D:
            m.AddBoolOr([x[(d, r)] for r in range(d)])          # every modulus used
        for q in odd_primes:
            if q in D:
                m.Add(x[(q, 0)] == 1)
            if 2 * q in D:
                for r in range(2 * q):
                    if r % q != 1: m.Add(x[(2 * q, r)] == 0)
        q0 = odd_primes[0]
        if 2 * q0 in D: m.Add(x[(2 * q0, 1)] == 1)
    else:
        for r in range(1, N):
            m.Add(x[(N, r)] == 0)
    for d in require:
        m.AddBoolOr([x[(d, r)] for r in range(d)])
    if slice_cuts:
        # Slice cuts for the largest prime q of N (squarefree N, M >= 3, normalize=True so the class of q is 0):
        # u_n = 1 iff n in Z_{N/q} is uncovered by the q-free progressions; |R| >= slice_cuts (a proven lower
        # bound on the minimum uncovered count of the q-free level); every slice c != 0 covers R with its
        # shapes, so the total size of its progressions is at least |R|.
        q = max(p for p in range(2, 200) if N % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1)))
        Lq = N // q
        Dfree = [d for d in D if Lq % d == 0]
        u = [m.NewBoolVar(f"u{n}") for n in range(Lq)]
        for n in range(Lq):
            m.AddBoolOr([x[(d, n % d)] for d in Dfree] + [u[n]])
            for d in Dfree: m.AddImplication(x[(d, n % d)], u[n].Not())
        R = sum(u)
        m.Add(R >= slice_cuts)
        for c in range(1, q):
            m.Add(sum((Lq // d) * x[(q * d, r)] for d in [e for e in divisors(Lq) if e > 1] for r in range(c, q * d, q)) >= R)
    if redundant:
        # implied counting constraints: for g=N/p every class c mod g needs total weight >= 1,
        # where modulus d contributes gcd(d,g)/d if its residue is = c mod gcd(d,g).
        primes = [q for q in range(2, 200) if N % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
        for q in primes:
            g = N // q
            agg = {}
            for d in D:
                h = math.gcd(d, g)
                for s_ in range(h):
                    agg[(d, s_)] = sum(x[(d, r)] for r in range(s_, d, h))  # residue r = s_ mod h
            for c in range(g):
                m.Add(sum((N * math.gcd(d, g) // d) * agg[(d, c % math.gcd(d, g))] for d in D) >= N)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    s.parameters.linearization_level = linearization
    st = s.Solve(m)
    name = s.StatusName(st)
    res = dict(N=N, M=M, status=name, seconds=round(time.time() - t, 1), moduli_available=len(D), normalize=normalize)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        system = sorted((d, r) for (d, r), v in x.items() if s.Value(v))
        assert len({d for d, _ in system}) == len(system)
        assert all(any(n % d == r for d, r in system) for n in range(N)), "coverage check failed"
        res["status"] = "SAT"
        res["system"] = [dict(modulus=d, residue=r) for d, r in system]
        res["moduli_used"] = len(system); res["min_modulus"] = min(d for d, _ in system)
    elif st == cp_model.INFEASIBLE:
        res["status"] = "UNSAT"
    else:
        res["status"] = "UNKNOWN"
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", type=int, required=True); ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--time", type=float, default=600.0); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output"); ap.add_argument("--redundant", action="store_true"); ap.add_argument("--normalize", action="store_true"); ap.add_argument("--slice-cuts", type=int, default=0, help="lower bound on the q-free residual; enables the slice cuts"); ap.add_argument("--lin", type=int, default=1)
    a = ap.parse_args()
    r = solve(a.N, a.M, a.time, a.workers, redundant=a.redundant, linearization=a.lin, normalize=a.normalize, slice_cuts=a.slice_cuts)
    print(json.dumps({k: v for k, v in r.items() if k != "system"}))
    if "system" in r: print("system:", [(e["residue"], e["modulus"]) for e in r["system"]])
    if a.output: json.dump(r, open(a.output, "w"), indent=1)
