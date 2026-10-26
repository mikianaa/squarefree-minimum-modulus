"""Lower numerical case counts until every incompatible-event clique fits."""
import argparse
import copy
import json
from pathlib import Path
from clique_case_constraints import strongest


def repair_pool(data):
    result = copy.deepcopy(data)
    R = data['input']['core_residual']
    n = len(data['input']['high_primes'])
    total = changed = 0
    largest = 0
    for group in result['group_scenarios']:
        h = group['axis']
        shape = [d for d in range(3 if h < n else 1, 211) if 210 % d == 0]
        rewards = [[-sum(row[:k+1]) for k in range(len(shape)+1)]
                   for row in data['baseline_event_multipliers'][h]]
        for scenario in group['scenarios']:
            remaining = [d for d in shape if d not in dict(scenario['fixed'])]
            rows = scenario['count_distributions']
            old, _ = strongest(R, remaining, rows)
            factor = min(1., (1-1e-12)/old) if old else 1.
            if factor < 1:
                scenario['count_distributions'] = [[1-factor+factor*row[0]]+[factor*v for v in row[1:]] for row in rows]
                changed += 1
            new, _ = strongest(R, remaining, scenario['count_distributions'])
            assert new <= 1+1e-9
            base = int(1 in remaining)
            hits = [sum(p % d == a for d, a in scenario['fixed']) for p in R]
            scenario['payoffs'] = [sum(v*rewards[x][hits[x]+base+j] for j, v in enumerate(row))
                                   for x, row in enumerate(scenario['count_distributions'])]
            for key in ('objective', 'numerical_dual'):
                scenario.pop(key, None)
            scenario.update(clique_compatibility=True, clique_repair_factor=factor,
                            clique_check_maximum=new,
                            scope='Numerical clique-compatible count reduction; resource means and prior pair inequalities are preserved. Not an actual assignment or rational certificate.')
            largest = max(largest, old)
            total += 1
        print(json.dumps(dict(axis=h, processed=total, repaired=changed)), flush=True)
    result['clique_repair'] = dict(scenarios=total, repaired=changed,
                                   maximum_previous_clique_sum=largest,
                                   verified=False, covering_obstruction=False)
    result['scope'] = 'Numerically clique-compatible scenarios, preserving case order and resource distributions.'
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    result = repair_pool(json.loads(Path(args.source).read_text()))
    result['clique_repair']['source'] = args.source
    Path(args.output).write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(json.dumps(result['clique_repair']))
