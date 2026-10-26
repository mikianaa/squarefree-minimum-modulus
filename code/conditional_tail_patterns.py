"""Generate or independently audit all 16 high-prime patterns for one core."""
import argparse
import copy
import json
from pathlib import Path
import time


def generate(seconds):
    from conditional_tail_lp import solve
    read = lambda name: json.loads(Path('data',name).read_text())
    previous = read('conditional-tail-core-all-patterns.json')
    base = previous['base_certificate']
    concentration = read('core-hit-concentration-7-15.json')
    incumbent = read('core-fibre-eight-plane.json')
    cases = []
    for row in previous['patterns']:
        pattern = row['pattern']
        if row['integer_lower']>0:
            result = copy.deepcopy(base)
            result['high_pattern'] = pattern
            result['certificate']['point_multipliers'] = row['point_multipliers']
            result['certificate']['numerator'] = row['numerator']
            result['certificate']['integer_lower'] = row['integer_lower']
            result['converged_to_tolerance'] = False
            result['generation'] = 'reused dual with all point columns repriced'
        else:
            incumbent['high_pattern'] = pattern
            result = solve(incumbent,concentration,seconds,aggregate=True,stop_after_positive=True)
        cases.append(result)
        Path('data/conditional-tail-core-pattern-certificates.json').write_text(json.dumps(dict(cases=cases),indent=2)+'\n')
        print(json.dumps(dict(pattern=pattern,integer_lower=result['certificate']['integer_lower'],cases_done=len(cases))),flush=True)


def audit_all():
    from audit_conditional_tail import audit
    from audit_fibre_concentration import gray_capacities
    start = time.monotonic()
    read = lambda name: json.loads(Path('data',name).read_text())
    data = read('conditional-tail-core-pattern-certificates.json')
    concentration = read('core-hit-concentration-7-15.json')
    incumbent = read('core-fibre-eight-plane.json')
    assert gray_capacities(concentration['core_residual'],concentration['moduli'])==concentration['two_layer_capacities']
    assert len(data['cases'])==16 and len({tuple(d['high_pattern']) for d in data['cases']})==16
    results = []
    for case in data['cases']:
        result = audit(case,concentration,incumbent,check_subsets=False)
        results.append({k:v for k,v in result.items() if k not in ('core_system','scope','core_subsets_checked')})
        print(json.dumps(dict(pattern=case['high_pattern'],lower=result['conditional_residual_lower'],verified=True)),flush=True)
    out = dict(verified=True,all_core_subsets=1<<len(concentration['core_residual']),cases=results,
               minimum_lower=min(r['conditional_residual_lower'] for r in results),
               all_16_patterns_excluded=all(r['excludes_this_core_and_pattern'] for r in results),
               seconds=round(time.monotonic()-start,3),
               scope='one specified 210-core only; not all cores and not a global eight-prime obstruction')
    Path('data/audit-conditional-tail-core-patterns.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cases'}))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit',action='store_true')
    ap.add_argument('--time',type=float,default=45)
    a = ap.parse_args()
    audit_all() if a.audit else generate(a.time)
