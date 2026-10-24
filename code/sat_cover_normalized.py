"""Pure SAT (CaDiCaL via python-sat) version of the normalized squarefree cover search.

Decides whether a covering system with distinct moduli d | N, d >= M exists, in the lossless normal
form of cpsat_cover.solve (squarefree N, M >= 3): every divisor >= M of N is used, the class of each
odd prime p is at residue 0, the class of 2p has residue 1 mod p, the parity of the class of 2 p_min
is fixed.  Memory: one clause per residue class mod N (length = number of moduli) plus at-most-one
constraints per modulus.  UNSAT is a proof of nonexistence.
Usage: python sat_cover_normalized.py --N 30030 --M 3 [--solver cadical153] [--output out.json]
"""
import argparse, json, time
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from cpsat_cover import divisors


def build(N, M):
    D = [d for d in divisors(N) if d >= M]
    var = {}; nv = 0
    for d in D:
        for r in range(d):
            nv += 1; var[(d, r)] = nv
    odd_primes = [q for q in range(3, N + 1) if N % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    clauses = []
    for d in D:
        clauses.append([var[(d, r)] for r in range(d)])                       # used
        lits = [var[(d, r)] for r in range(d)]
        enc = CardEnc.atmost(lits=lits, bound=1, top_id=nv, encoding=EncType.seqcounter)
        nv = max(nv, enc.nv); clauses.extend(enc.clauses)
    for n in range(N):
        clauses.append([var[(d, n % d)] for d in D])
    for q in odd_primes:
        if q in D: clauses.append([var[(q, 0)]])
        if 2 * q in D:
            for r in range(2 * q):
                if r % q != 1: clauses.append([-var[(2 * q, r)]])
    q0 = odd_primes[0]
    if 2 * q0 in D: clauses.append([var[(2 * q0, 1)]])
    return D, var, clauses, nv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True); ap.add_argument("--M", type=int, default=3)
    ap.add_argument("--solver", default="cadical153"); ap.add_argument("--output")
    a = ap.parse_args()
    t = time.time()
    D, var, clauses, nv = build(a.N, a.M)
    tb = time.time() - t
    with Solver(name=a.solver, bootstrap_with=clauses) as s:
        sat = s.solve()
        model = s.get_model() if sat else None
    res = dict(N=a.N, M=a.M, moduli=len(D), variables=nv, clauses=len(clauses), build_seconds=round(tb, 1),
               seconds=round(time.time() - t, 1), status="SAT" if sat else "UNSAT", solver=a.solver)
    if sat:
        pos = set(l for l in model if l > 0)
        system = sorted((d, r) for (d, r), v in var.items() if v in pos)
        assert all(any(n % d == r for d, r in system) for n in range(a.N))
        res["system"] = [dict(modulus=d, residue=r) for d, r in system]
    print(json.dumps({k: v for k, v in res.items() if k != "system"}))
    if a.output: json.dump(res, open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
