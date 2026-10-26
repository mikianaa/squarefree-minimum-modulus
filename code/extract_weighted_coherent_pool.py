"""Evaluate saved weighted duals exactly on an explicitly named finite pool."""
import argparse
import copy
from fractions import Fraction
import json
from pathlib import Path


def extract(data,pool):
    assert data['input']==pool['input']
    source = data['input'];R = source['core_residual'];n = len(source['high_primes'])
    configurations = pool['result']['group_columns']
    counts = []
    for h,rows in enumerate(configurations):
        shape = [d for d in range(3 if h<n else 1,211) if 210%d==0]
        counts.append([[sum(x%d==a for a,d in zip(residues,shape)) for x in R] for residues in rows])
    best,best_value = None,None
    for index,original in enumerate(data['result']['pricing_history']):
        c = copy.deepcopy(original)
        maxima = [max(sum(e['weight'] for e in priced['events'] if row[e['point']]>=e['threshold']) for row in rows)
                  for rows,priced in zip(counts,c['group_pricing'])]
        numerator = sum(c['point_multipliers'])-c['point_weight_denominator']*sum(maxima)
        value = Fraction(numerator,c['denominator'])
        if best_value is None or value>best_value:
            c.update(pool_upper_bounds=maxima,pool_sizes=list(map(len,configurations)),pool_numerator=numerator,
                     source_pricing_index=index)
            best,best_value = c,value
    assert best is not None
    return dict(input=source,result=dict(group_columns=configurations,pricing_history=[best],
                scope='dual evaluated on the explicitly named finite pool; no unrestricted exclusion asserted'))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--pool',required=True);ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = extract(json.loads(Path(args.source).read_text()),json.loads(Path(args.pool).read_text()))
    result['provenance'] = dict(dual_source=args.source,pool_source=args.pool)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    c = result['result']['pricing_history'][0]
    print(json.dumps(dict(pool_numerator=c['pool_numerator'],denominator=c['denominator'],
                          capacity_numerator=c['capacity_numerator'],source_pricing_index=c['source_pricing_index'])))
