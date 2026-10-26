"""Numerical, per-prime stop-loss optimization for cutoffs C0=6 or C0=8.

This is a candidate generator, not a proof verifier. The exported distortions
are rounded to denominator 10**6 and clipped to [0, 1/2]. Verify the resulting
schedule with BOTH verify_squarefree.py and audit_verify.py.

For a fixed prefix and fixed later controls, put r=1/(p*(1-delta)).
The current charge is E[max(0, 1+r*(A-p))], with A=2**N-1 for
5<=p<=C0 and A=2**N thereafter. Future cost is affine in r. Thus the
one-coordinate minimum occurs at a breakpoint or an interval endpoint.
Forward coordinate sweeps use backwards-computed continuation values.
This does NOT establish global optimality.

Reproduce the search producing the upper bound 7:
  python code/optimize_stoploss_breakpoints.py data/schedule8-float-only.json OUT.json
Then verify OUT.json with --p-exact 3000000 and --n-stages 1000000.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from verify_squarefree import primes_to


def candidate_probabilities(p, modulus_cutoff=8):
    """Search the full Bernoulli-admissible interval; export clips delta<=1/2."""
    shift = 1 if p <= modulus_cutoff else 0
    values = (2.0 ** np.arange(32) - shift) / p
    deltas = np.unique(np.r_[0.0, values[(values >= 0) & (values <= 1 - 1 / p)], 1 - 1 / p])
    return 1 / (p * (1 - deltas))


def stage_costs(p, probabilities, states, modulus_cutoff=8):
    counts = 2.0 ** np.arange(states) - (1 if p <= modulus_cutoff else 0)
    return np.maximum(0.0, 1 + probabilities[:, None] * (counts - p))


def quadratic_continuation(primes, cutoff):
    """Numerical optimal quadratic continuation, ending in the existing CFT tail."""
    last = primes[-1]
    log = math.log(last)
    polynomial = log**5 + 5*log**4 + 20*log**3 + 60*log**2 + 120*log + 120
    value = 2*(1+3/(2*log))*math.exp(6*(-math.log(log)+1/log**2))*polynomial/last
    later = [p for p in primes if p > cutoff]
    deltas = []
    for p in reversed(later):
        root = math.sqrt(1 + 12*p*value)
        deltas.append(1 / (1 + root))
        value += (1 + root)**2 / (4*p*p)
    return value, deltas[::-1]


def compress_rounded_schedule(deltas):
    """Only the exported rational schedule is used by the proof verifiers."""
    integers = [min(500000, max(0, round(float(d)*1000000))) for d in deltas]
    edges, params = [], []
    for i, value in enumerate(integers):
        if i == len(integers)-1 or integers[i+1] != value:
            edges.append(i+1)
            params.append(value / 1000000)
    return edges, params


def optimize(initial_path, output_path, cutoff=3000000, sweeps=4, states=50, modulus_cutoff=8):
    stages = 1000000
    all_primes = primes_to(16000000)[:stages]
    assert len(all_primes) == stages and 11 <= cutoff < all_primes[-1]
    primes = [p for p in all_primes if 5 <= p <= cutoff]
    initial_bytes = Path(initial_path).read_bytes()
    initial = json.loads(initial_bytes)
    assert initial['C0'] in (6, 8) and modulus_cutoff in (6, 8)
    ds = np.array(initial['params'])[np.searchsorted(initial['edges'], np.arange(3, len(primes)+3))]
    probabilities = 1 / (np.array(primes)*(1-ds))
    coefficient, later_deltas = quadratic_continuation(all_primes, cutoff)
    terminal = coefficient * 4.0**np.arange(states)
    start = np.zeros(states)
    start[:3] = [1/3, 1/2, 1/6]  # Bernoulli(1/2)+Bernoulli(1/3)
    total = None
    for sweep in range(sweeps):
        continuation = np.zeros((len(primes)+1, states))
        continuation[-1] = terminal
        for i in range(len(primes)-1, -1, -1):
            v = continuation[i+1]
            r = probabilities[i]
            continuation[i] = stage_costs(primes[i], np.array([r]), states, modulus_cutoff)[0] + (1-r)*v + r*np.r_[v[1:], 0.0]
        print(f'sweep {sweep}: starting numerical total {start @ continuation[0]:.15f}', flush=True)
        distribution = start.copy()
        total = 0.0
        for i, p in enumerate(primes):
            options = candidate_probabilities(p, modulus_cutoff)
            v = continuation[i+1]
            objectives = stage_costs(p, options, states, modulus_cutoff) + (1-options[:, None])*v + options[:, None]*np.r_[v[1:], 0.0]
            r = options[np.argmin(objectives @ distribution)]
            probabilities[i] = r
            total += float(distribution @ stage_costs(p, np.array([r]), states, modulus_cutoff)[0])
            distribution = (1-r)*distribution + r*np.r_[0.0, distribution[:-1]]
        total += float(distribution @ terminal)
    raw_deltas = [0.0, 0.0] + (1-1/(np.array(primes)*probabilities)).tolist() + later_deltas
    edges, params = compress_rounded_schedule(raw_deltas)
    result = dict(C0=modulus_cutoff, edges=edges, params=params,
                  origin='per-prime breakpoint optimization; export clipped to delta<=1/2 and rounded to 1e-6',
                  numerical_total_before_clipping=total,
                  search_parameters=dict(p_exact=cutoff, n_stages=stages, states=states, sweeps=sweeps,
                                         initial_file=str(initial_path), initial_sha256=hashlib.sha256(initial_bytes).hexdigest()),
                  clipped_parameters=sum(d > .5 for d in raw_deltas))
    Path(output_path).write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('edges', 'params')}, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('initial')
    parser.add_argument('output')
    parser.add_argument('--p-exact', type=int, default=3000000)
    parser.add_argument('--sweeps', type=int, default=4)
    parser.add_argument('--C0', type=int, choices=(6, 8), default=8)
    args = parser.parse_args()
    optimize(args.initial, args.output, args.p_exact, args.sweeps, modulus_cutoff=args.C0)
