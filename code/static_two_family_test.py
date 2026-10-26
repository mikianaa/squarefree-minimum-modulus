"""Direction 2 of notes/direction-analysis-2026-10-25.md: the static two-family
model, tested directly (not just its "G-residues zero" aligned special case
from notes/aligned-construction-test-2026-10-26.md).

Given a period L, a next prime q, and an EXPLICIT residual R (the actual
uncovered points of a known family F over Z_L, moduli >= 3), decide: can a
single residue be assigned to every divisor d > 1 of L (the "G" family, no
alignment restriction -- any residue, not just 0) so that every point of R
gets covered by at least q-1 of these 31 (for L=2310) progressions?

This is a small, exact CP-SAT feasibility check (|R| points, divisors(L)-1
boolean "which residue" variables), decidable regardless of how the
aligned-only test fares. It is a genuine test of direction 2's extra power
over direction 3, not just a restatement of it.

Usage: python code/static_two_family_test.py --source data/min-uncovered-lcm2310-feasible121.json --q 13
"""
import argparse
import json
from pathlib import Path

from ortools.sat.python import cp_model


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n:
                ds.append(n // i)
        i += 1
    return sorted(ds)


def run(source_path, q, time_limit=120.0):
    source = json.loads(Path(source_path).read_text())
    L = source["N"]
    R = source["uncovered_residues"]
    assert len(R) == source["uncovered"]
    ds = [d for d in divisors(L) if d > 1]
    need = q - 1

    m = cp_model.CpModel()
    y = {}
    for d in ds:
        y[d] = [m.NewBoolVar(f"y_{d}_{a}") for a in range(d)]
        m.AddExactlyOne(y[d])
    for x in R:
        m.Add(sum(y[d][x % d] for d in ds) >= need)

    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = 8
    status = s.Solve(m)
    name = s.StatusName(status)
    result = dict(L=L, q=q, need=need, R_size=len(R), num_divisors=len(ds), status=name)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        assignment = {d: int(next(a for a in range(d) if s.Value(y[d][a]))) for d in ds}
        # independent check: recompute multiplicities from the assignment, no solver trust
        for x in R:
            mult = sum(1 for d in ds if x % d == assignment[d])
            assert mult >= need, (x, mult, need)
        result["status"] = "SAT"
        result["G"] = [dict(modulus=d, residue=a) for d, a in sorted(assignment.items())]
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--time", type=float, default=120.0)
    ap.add_argument("--output")
    args = ap.parse_args()
    out = run(args.source, args.q, args.time)
    print(json.dumps({k: v for k, v in out.items() if k != "G"}, indent=2))
    if args.output:
        Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
