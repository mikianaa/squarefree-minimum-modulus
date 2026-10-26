"""Search simultaneous actual core projections for the local F<=0 conditions.

No high-prime labels are chosen here. A feasible result is a relaxation
witness, not a covering. Failure of this heuristic is not an obstruction.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np


class JointProjections:
    def __init__(self, source):
        self.source = source
        self.points = np.asarray(source['core_residual'], dtype=np.int64)
        self.n = len(source['high_primes'])
        self.groups = [s for k in range(1, self.n+1) for s in itertools.combinations(range(self.n), k)]
        self.shapes = [[d for d in range(3 if h < self.n else 1, 211) if 210 % d == 0]
                       for h in range(len(self.groups))]
        self.D = np.asarray([[q-1-int(b == x % 2) for x in self.points]
                             for q, b in zip(source['high_primes'], source['high_pattern'])])
        self.masks = [sum(1 << i for i in s) for s in self.groups]
        self.full = (1 << self.n)-1
        self.complements = [self.full ^ m for m in self.masks[self.n:]]
        self.resources = []
        for h, shape in enumerate(self.shapes):
            for j, d in enumerate(shape):
                residues = sorted(set(map(int, self.points % d)))
                hits = np.asarray([self.points % d == a for a in residues], dtype=np.int64)
                self.resources.append((h, j, d, residues, hits))

    def counts(self, assignments):
        return np.asarray([[sum(int(x % d == a) for d, a in zip(shape, row)) for x in self.points]
                           for shape, row in zip(self.shapes, assignments)], dtype=np.int64)

    def monomials(self, counts):
        free = np.maximum(0, self.D-counts[:self.n])
        products = np.ones((1 << self.n, len(self.points)), dtype=np.int64)
        for mask in range(1, 1 << self.n):
            bit = mask & -mask
            products[mask] = products[mask ^ bit]*free[bit.bit_length()-1]
        return products

    def evaluate(self, counts):
        products = self.monomials(counts)
        return products[self.full]-np.sum(counts[self.n:]*products[self.complements], axis=0)

    def increment(self, counts, h):
        """Change in F from increasing group h's count by one at each point."""
        products = self.monomials(counts)
        if h >= self.n:
            return -products[self.full ^ self.masks[h]]
        bit = 1 << h
        derivative = products[self.full ^ bit].copy()
        for j, mask in enumerate(self.complements, self.n):
            if mask & bit:
                derivative -= counts[j]*products[mask ^ bit]
        return np.where(counts[h] < self.D[h], -derivative, 0)

    def group_costs(self, counts, h):
        costs = []
        state = counts.copy()
        for k in range(len(self.shapes[h])+1):
            state[h] = k
            costs.append(np.maximum(self.evaluate(state), 0))
        table = np.asarray(costs)
        # With other groups fixed, positive F is nonincreasing in a group's
        # count. This allows the existing nonnegative-event optimizer.
        assert np.all(table[1:] <= table[:-1])
        return table


def search(source, seconds=120, random_seed=1, initial=None):
    model = JointProjections(source)
    rng = np.random.default_rng(random_seed)
    assignments = [[int(rng.choice(model.points % d)) for d in shape] for shape in model.shapes]
    if initial is not None:
        assignments = [list(row) for row in initial]
    counts = model.counts(assignments)
    values = model.evaluate(counts)
    # Raw positive volume is integer and is the reported objective. Search
    # weights can emphasize persistent violations without changing that test.
    best = int(np.maximum(values, 0).sum())
    best_assignments = [row[:] for row in assignments]
    weights = np.ones(len(model.points))
    start = time.monotonic()
    passes = 0
    history = [dict(pass_index=0, positive_volume=best, violated_points=int((values > 0).sum()))]
    stale = 0
    while best > 0 and time.monotonic()-start < seconds:
        changed = 0
        for rid in rng.permutation(len(model.resources)):
            h, j, d, residues, hits = model.resources[rid]
            old = (model.points % d == assignments[h][j]).astype(np.int64)
            counts[h] -= old
            base = model.evaluate(counts)
            delta = model.increment(counts, h)
            # Each option changes independent point terms, so scores are dot
            # products of its incidence vector with the exact one-step gain.
            gains = weights*(np.maximum(base, 0)-np.maximum(base+delta, 0))
            scores = hits @ gains
            candidates = np.flatnonzero(scores >= scores.max()-1e-9)
            pick = int(rng.choice(candidates))
            changed += int(assignments[h][j] != residues[pick])
            assignments[h][j] = residues[pick]
            counts[h] += hits[pick]
            values = base+hits[pick]*delta
        passes += 1
        assert np.array_equal(values, model.evaluate(counts))
        value = int(np.maximum(values, 0).sum())
        if value < best:
            best = value
            best_assignments = [row[:] for row in assignments]
            row = dict(pass_index=passes, positive_volume=best, violated_points=int((values > 0).sum()),
                       maximum_F=int(values.max()), seconds=time.monotonic()-start)
            history.append(row)
            print(json.dumps(row), flush=True)
            stale = 0
        else:
            stale += 1
        if best == 0:
            break
        # Penalize persistent positive fibres, then perturb a few choices to
        # leave a coordinatewise minimum. Occasionally restart from the best.
        weights = .97*weights+.03+0.1*(values > 0)
        if stale >= 3:
            if passes % 7 == 0:
                assignments = [row[:] for row in best_assignments]
                counts = model.counts(assignments)
                weights[:] = 1
            for rid in rng.choice(len(model.resources), size=min(3+stale, 60), replace=False):
                h, j, d, residues, hits = model.resources[rid]
                counts[h] -= (model.points % d == assignments[h][j])
                pick = int(rng.integers(len(residues)))
                assignments[h][j] = residues[pick]
                counts[h] += hits[pick]
            values = model.evaluate(counts)
    counts = model.counts(best_assignments)
    values = model.evaluate(counts)
    return dict(input=source, group_assignments=best_assignments, local_F=list(map(int, values)),
                positive_volume=int(np.maximum(values, 0).sum()), passes=passes, random_seed=random_seed,
                seconds=time.monotonic()-start, history=history, verified=False,
                scope='Simultaneous actual core projections only; high-prime labels and covering untested. Heuristic failure proves no impossibility.')


def polish(data, seconds=60, seconds_per_group=1):
    """Replace entire group assignments, always retaining an explicit incumbent."""
    from group_event_coherence import maximize
    model = JointProjections(data['input'])
    assignments = [row[:] for row in data['group_assignments']]
    counts = model.counts(assignments)
    initial = int(np.maximum(model.evaluate(counts), 0).sum())
    best = initial
    start = time.monotonic()
    diagnostics = []
    while best > 0 and time.monotonic()-start < seconds:
        before = best
        for h, shape in enumerate(model.shapes):
            remaining = seconds-(time.monotonic()-start)
            if remaining <= 0:
                break
            table = model.group_costs(counts, h)
            gains = table[:-1]-table[1:]
            events = [(x, t+1) for t, row in enumerate(gains) for x, w in enumerate(row) if w]
            weights = [int(gains[t-1, x]) for x, t in events]
            if not events:
                continue
            result = maximize(list(map(int, model.points)), shape, events, weights,
                              min(seconds_per_group, remaining), hint=assignments[h])
            proposal = result.get('residues', assignments[h])
            old = assignments[h]
            assignments[h] = list(proposal)
            candidate = model.counts(assignments)
            value = int(np.maximum(model.evaluate(candidate), 0).sum())
            if value <= best:
                counts = candidate
            else:
                assignments[h] = old
                value = best
            row = dict(axis=h, before=best, after=value, solver_status=result['status'])
            diagnostics.append(row)
            if value < best:
                print(json.dumps(row), flush=True)
            best = value
            if best == 0:
                break
        if best >= before:
            break
    values = model.evaluate(counts)
    return dict(data, group_assignments=assignments, local_F=list(map(int, values)),
                positive_volume=best, verified=False, group_polish=dict(initial_positive_volume=initial,
                final_positive_volume=best, seconds=time.monotonic()-start, diagnostics=diagnostics,
                solver_bounds_used_as_proofs=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--initial')
    ap.add_argument('--time', type=float, default=120)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--output', required=True)
    ap.add_argument('--polish-time', type=float, default=0)
    ap.add_argument('--group-time', type=float, default=1)
    args = ap.parse_args()
    source = json.loads(Path(args.source).read_text())['input']
    initial = json.loads(Path(args.initial).read_text())['group_assignments'] if args.initial else None
    result = search(source, args.time, args.seed, initial)
    if args.polish_time:
        result = polish(result, args.polish_time, args.group_time)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input','group_assignments','history')}))
