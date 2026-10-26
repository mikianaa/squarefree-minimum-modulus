"""Exact squarefree covering search in CRT boxes after prime/twice-prime anchors.

For minimum modulus 3 and N = 2 * product(odd primes), normalize p to
label 0, 2p to label 1 in parity b_p, and b_3 (or the first prime) to 0.
The uncovered anchor boxes have coordinate alphabets 1..p-1, excluding 1
in parity b_p. Remaining resources are d and 2d for odd composite d|N.

All permutations of labels 2..p-1 are still symmetries. First-use ordering
of these labels, across a fixed resource order, removes these symmetries.
This is stronger than affine normalization and remains lossless by CRT.
Every solution is verified on the entire integer period, independently of
the reduced coverage constraints. A fixed anchor pattern is only ONE case.
"""
import argparse
import itertools
import json
import math
import time
from ortools.sat.python import cp_model


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            n //= p
            if n % p == 0:
                raise ValueError("N must be squarefree")
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def crt(labels, primes):
    n = math.prod(primes)
    return sum(a * (n // p) * pow(n // p, -1, p)
               for a, p in zip(labels, primes)) % n


def verify(n, system):
    assert len({d for d, _ in system}) == len(system)
    assert all(d >= 3 and n % d == 0 for d, _ in system)
    cov = bytearray(n)
    for d, a in system:
        for x in range(a, n, d):
            cov[x] = 1
    return [x for x, v in enumerate(cov) if not v]


def build(n, pattern, minimize=False, symmetry=True, residual_lb=0, pair_relaxation=False):
    primes = primes_of(n)
    assert primes[0] == 2 and len(primes) >= 2
    ps = primes[1:]
    assert len(pattern) == len(ps) and pattern[0] == 0
    assert all(b in (0, 1) for b in pattern)
    mdl = cp_model.CpModel()
    boxes = [[tuple(v for v in range(1, p) if v != 1 or pattern[i] != e)
              for i, p in enumerate(ps)] for e in (0, 1)]
    points = [(e, a) for e in (0, 1) for a in itertools.product(*boxes[e])]
    covering = [[] for _ in points]
    resources = []
    omitted_capacity = 0
    for k in range(2, len(ps) + 1):
        for subset in itertools.combinations(range(len(ps)), k):
            if pair_relaxation and k > 2:
                capacities = [math.prod(len(boxes[e][i]) for i in range(len(ps)) if i not in subset) for e in (0,1)]
                omitted_capacity += sum(capacities) + max(capacities)
                continue
            for even in (False, True):
                d = math.prod(ps[i] for i in subset) * (2 if even else 1)
                resources.append((d, subset, even))
    resources.sort()
    by_coordinate = [[] for _ in ps]
    choices = []
    for d, subset, even in resources:
        placements = []
        lookup = {}
        for e in ((0, 1) if even else (None,)):
            alphabets = [boxes[e][i] if even else range(1, ps[i]) for i in subset]
            for labels in itertools.product(*alphabets):
                if not even and not any(all(a in boxes[j][i] for i, a in zip(subset, labels))
                                        for j in (0, 1)):
                    continue  # This placement lies entirely in anchor classes.
                v = mdl.NewBoolVar(f"x{d}_{e}_{'_'.join(map(str, labels))}")
                lookup[(e, labels)] = v
                placements.append((e, labels, v))
        mdl.AddExactlyOne(v for _, _, v in placements)
        for pos, i in enumerate(subset):
            label = mdl.NewIntVar(1, ps[i] - 1, f"label{d}_{ps[i]}")
            mdl.Add(label == sum(a[pos] * v for _, a, v in placements))
            by_coordinate[i].append(label)
        for j, (e, a) in enumerate(points):
            key = (e if even else None, tuple(a[i] for i in subset))
            covering[j].append(lookup[key])
        choices.append((d, subset, even, placements))
    if symmetry:
        for i, labels in enumerate(by_coordinate):
            previous = 1
            for j, label in enumerate(labels):
                mdl.Add(label <= previous + 1)
                current = mdl.NewIntVar(1, ps[i] - 1, f"maximum{ps[i]}_{j}")
                mdl.AddMaxEquality(current, [previous, label])
                previous = current
    missed = []
    for j, lits in enumerate(covering):
        if minimize or pair_relaxation:
            u = mdl.NewBoolVar(f"u{j}")
            mdl.AddBoolOr(lits + [u])
            missed.append(u)
        else:
            mdl.AddBoolOr(lits)
    if minimize:
        mdl.Add(sum(missed) >= residual_lb)
        mdl.Minimize(sum(missed))
    elif pair_relaxation:
        mdl.Add(sum(missed) <= omitted_capacity)
    return mdl, ps, points, choices, omitted_capacity


def solve(n, pattern, seconds=60, workers=8, minimize=False, symmetry=True, residual_lb=0, pair_relaxation=False):
    start = time.monotonic()
    mdl, ps, points, choices, omitted_capacity = build(n, pattern, minimize, symmetry, residual_lb, pair_relaxation)
    built = time.monotonic()
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    status = solver.Solve(mdl)
    result = dict(N=n, pattern=list(pattern), status=solver.StatusName(status),
                  minimize=minimize, label_symmetry=symmetry, points=len(points),
                  pair_relaxation=pair_relaxation, omitted_capacity=omitted_capacity,
                  build_seconds=round(built-start, 3),
                  seconds=round(time.monotonic()-start, 3))
    if minimize:
        result["lower_bound"] = math.ceil(solver.BestObjectiveBound() - 1e-7)
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        system = [(p, 0) for p in ps]
        system += [(2*p, crt([b, 1], [2, p])) for p, b in zip(ps, pattern)]
        for d, subset, even, placements in choices:
            e, labels, _ = next(t for t in placements if solver.Value(t[2]))
            qs = [ps[i] for i in subset]
            a = crt(([e] + list(labels)) if even else labels, ([2] + qs) if even else qs)
            system.append((d, a))
        uncovered = verify(n, system)
        if minimize:
            # Auxiliary missed indicators need not all be tight at a nonoptimal
            # incumbent. Report independently verified coverage as the upper bound.
            assert len(uncovered) <= round(solver.ObjectiveValue())
            if status == cp_model.OPTIMAL:
                assert len(uncovered) == round(solver.ObjectiveValue())
        elif not pair_relaxation:
            assert not uncovered
        else:
            assert len(uncovered) <= omitted_capacity
        result["uncovered"] = len(uncovered)
        result["system"] = [dict(modulus=d, residue=a) for d, a in sorted(system)]
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("N", type=int)
    ap.add_argument("--pattern", help="anchor parity bits, starting with 0; default: every pattern")
    ap.add_argument("--time", type=float, default=60)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--minimize", action="store_true")
    ap.add_argument("--no-symmetry", action="store_true")
    ap.add_argument("--residual-lb", type=int, default=0)
    ap.add_argument("--pair-relaxation", action="store_true", help="replace higher-order classes by their total capacity; feasible does NOT imply a cover")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    k = len(primes_of(args.N)) - 1
    patterns = [tuple(map(int, args.pattern))] if args.pattern else [(0,)+b for b in itertools.product((0, 1), repeat=k-1)]
    with open(args.output, "w") as out:
        for pattern in patterns:
            result = solve(args.N, pattern, args.time, args.workers, args.minimize, not args.no_symmetry, args.residual_lb, args.pair_relaxation)
            out.write(json.dumps(result) + "\n")
            out.flush()
            print(json.dumps({k: v for k, v in result.items() if k != "system"}), flush=True)
