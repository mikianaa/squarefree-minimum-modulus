"""Reproduce the 30030 count refinement using a saved 2310 solver optimum.

The lower-bound input remains solver-dependent. This script does not replace
the solver's optimality proof with an integer-only proof of that input.
"""
import argparse
import json
from pathlib import Path
from parity_capacity import pareto, extend, can_cover_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pattern',default='0111')
    ap.add_argument('--new-parity',type=int,choices=(0,1),default=1)
    ap.add_argument('--output',required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    pattern = list(map(int,a.pattern))
    entries = json.loads((root/'data/min-uncovered-2310-five-patterns-optimal.json').read_text())
    source = next(e for e in entries if e['pattern']==pattern)
    assert source['status']=='OPTIMAL' and source['uncovered']==source['lower_bound']
    minimum = source['lower_bound']
    rows = json.loads((root/'data/parity-capacity-30030-certificate.json').read_text())
    old = next(e['frontier'] for e in rows if e['pattern']==pattern)
    refined = pareto((max(x,t),max(y,minimum-t)) for x,y in old for t in range(minimum+1))
    result_frontier = extend((3,5,7,11),pattern,refined,13,a.new_parity)
    extended_pattern = pattern+[a.new_parity]
    survivors = [sum(can_cover_counts((3,5,7,11,13),extended_pattern,r,17,b) for r in result_frontier) for b in (0,1)]
    result = dict(pattern=extended_pattern,prefix_minimum=minimum,
                  source='CP-SAT pattern optimum plus audited parity-count propagation',
                  uncovered_lower=min(map(sum,result_frontier)),frontier=result_frontier,
                  extension17_survivors=survivors,
                  verdict='not excluded' if any(survivors) else 'both extension patterns excluded by capacity')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='frontier'}),flush=True)


if __name__ == '__main__':
    main()
