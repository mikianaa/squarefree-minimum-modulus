"""Fit a fixed point witness with explicit mixtures of whole-group choices.

Configuration pools are shared only between groups with identical resources.
Each small LP maximizes uniform domination of that group's event marginals.
A feasible domination of all target events gives a directly checkable coherent
witness. Failure using a finite pool does not exclude omitted configurations.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
from ortools.linear_solver import pywraplp
from audit_coherent_groups import audit
from two_configuration_mixtures import structure,cost


def fit(source,witness,pool,pricing_seconds=0):
    start = time.monotonic()
    R,n,scale = source['core_residual'],len(source['high_primes']),witness['denominator']
    groups,shapes = structure(source)
    caches = []
    for lo,hi in ((0,n),(n,len(groups))):
        distinct = {}
        for rows in pool['group_columns'][lo:hi]:
            for residues in rows:
                counts = tuple(sum(x%d==a for d,a in zip(shapes[lo],residues)) for x in R)
                distinct.setdefault(counts,residues)
        caches.append(list(distinct.items()))
    hist = [[[0]*(len(shape)+1) for x in R] for shape in shapes]
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[h][row['point_index']][v] += row['weight']
    saved_groups,reports,priced_columns = [],[],[[] for h in groups]
    for h,shape in enumerate(shapes):
        configurations = caches[h>=n]
        events = [(x,t,sum(line[t:])/scale) for x,line in enumerate(hist[h])
                  for t in range(1,len(line)) if line[t]>0]
        solver = pywraplp.Solver.CreateSolver('GLOP')
        gamma = solver.NumVar(0,1,'gamma')
        variables = [solver.NumVar(0,1,str(j)) for j in range(len(configurations))]
        solver.Add(sum(variables)==1)
        constraints = []
        for x,t,target in events:
            constraints.append(solver.Add(target*gamma<=sum(var for (counts,residues),var in zip(configurations,variables) if counts[x]>=t)))
        solver.Maximize(gamma)
        status = solver.Solve()
        assert status==pywraplp.Solver.OPTIMAL
        values = [max(Fraction(0),Fraction(str(v.solution_value()))) for v in variables]
        total = sum(values)
        weights = [int(v*scale//total) for v in values]
        weights[values.index(max(values))] += scale-sum(weights)
        rows = [dict(residues=residues,weight=w) for (counts,residues),w in zip(configurations,weights) if w]
        saved_groups.append(dict(axis=h,assignments=rows))
        reports.append(dict(axis=h,domination=gamma.solution_value(),configurations=len(configurations),
                            active_configurations=len(rows),events=len(events)))
        if pricing_seconds and gamma.solution_value()<1-1e-8:
            from group_event_coherence import maximize
            weights = [max(0,round(c.dual_value()*10**6)) for c in constraints]
            priced = maximize(R,shape,[(x,t) for x,t,target in events],weights,pricing_seconds,collect=True)
            if 'counts' in priced:
                priced_columns[h] = [row['residues'] for row in priced.get('candidate_assignments',[])]+[priced['residues']]
                required = sum(w*sum(hist[h][x][t:]) for w,(x,t,target) in zip(weights,events))
                reports[-1]['pricing'] = dict(objective=priced['objective'],upper_bound=priced['upper_bound'],
                                              target_numerator=required,denominator=scale,
                                              configurations_found=len(priced_columns[h]))
                if required>priced['upper_bound']*scale:
                    reports[-1]['candidate_cut'] = dict(axis=h,events=[dict(point=x,threshold=t,weight=w)
                          for (x,t,target),w in zip(events,weights) if w],bound=priced['upper_bound'],
                          violation_numerator=required-priced['upper_bound']*scale,
                          witness_denominator=scale,weight_scale=10**6,pricing=priced)
        print(json.dumps(reports[-1]),flush=True)
    supply = [[[0]*(len(shape)+1) for x in R] for shape in shapes]
    for group,shape in zip(saved_groups,shapes):
        for row in group['assignments']:
            for x,v in enumerate(R):
                hits = sum(v%d==a for d,a in zip(shape,row['residues']))
                for t in range(1,hits+1):
                    supply[group['axis']][x][t] += row['weight']
    support,repairs = [],[]
    for x in range(len(R)):
        selected = [r for r in witness['support'] if r['point_index']==x]
        ratio = Fraction(1)
        for h,shape in enumerate(shapes):
            for t in range(1,len(shape)+1):
                used = sum(r['weight'] for r in selected if r['counts'][h]>=t)
                if used:
                    ratio = min(ratio,Fraction(supply[h][x][t],used))
        mass = 0
        for row in selected:
            w = row['weight']*ratio.numerator//ratio.denominator
            if w:
                support.append(dict(row,weight=w))
                mass += w
        if mass<scale:
            support.append(dict(point_index=x,counts=[0]*len(groups),weight=scale-mass))
        repairs.append([ratio.numerator,ratio.denominator])
    exact = dict(denominator=scale,groups=saved_groups,support=support,
                 objective_numerator=sum(r['weight']*cost(source,R[r['point_index']],r['counts']) for r in support))
    result = dict(coherent_primal_witness=exact,group_reports=reports,point_repair_ratios=repairs,
                  group_columns=[[row['residues'] for row in group['assignments']]+new for group,new in zip(saved_groups,priced_columns)],
                  seconds=time.monotonic()-start,
                  scope='finite-pool marginal fitting; a positive repaired witness is not a lower bound')
    result['audit'] = audit(dict(input=source,result=result))
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--pool',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--pricing-time',type=float,default=0)
    args = ap.parse_args()
    data,pool = [json.loads(Path(p).read_text()) for p in (args.source,args.pool)]
    assert data['input']['core_residual']==pool['input']['core_residual']
    result = fit(data['input'],data['result']['point_primal_witness'],pool['result'],args.pricing_time)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(seconds=result['seconds'],audit=result['audit'])))
