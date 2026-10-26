"""Solver bounds on weighted residual counts in the two parity layers.

Unlike minimizing the total alone, these constrain the realizable tradeoff
between layers. Lower bounds depend on CP-SAT, not independently checked
proof logs. Every incumbent is checked using ordinary integer congruences.
"""
import argparse
import json
import math
from pathlib import Path
import time
from ortools.sat.python import cp_model
from parity_cover import build,crt,verify


def solve_weighted(n,pattern,weights,seconds=60,workers=4):
    assert len(weights)==2 and all(isinstance(w,int) and w>0 for w in weights)
    start = time.monotonic()
    model,primes,points,choices,_ = build(n,pattern,minimize=True)
    indicators = {v.name:model.GetBoolVarFromProtoIndex(i)
                  for i,v in enumerate(model.proto.variables) if v.name.startswith('u')}
    assert len(indicators)==len(points)
    model.Minimize(sum(weights[e]*indicators[f'u{j}'] for j,(e,a) in enumerate(points)))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.linearization_level = 2
    status = solver.Solve(model)
    result = dict(N=n,pattern=list(pattern),weights=list(weights),status=solver.StatusName(status),
                  seconds=round(time.monotonic()-start,3),
                  scope='weighted lower bound depends on CP-SAT; incumbent verified over full period')
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE,cp_model.UNKNOWN):
        result['weighted_lower'] = max(0,math.ceil(solver.BestObjectiveBound()-1e-7))
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        system = [(p,0) for p in primes]+[(2*p,crt([b,1],[2,p])) for p,b in zip(primes,pattern)]
        for d,subset,even,placements in choices:
            e,labels,_ = next(t for t in placements if solver.Value(t[2]))
            qs = [primes[i] for i in subset]
            a = crt(([e]+list(labels)) if even else labels,([2]+qs) if even else qs)
            system.append((d,a))
        missed = verify(n,system)
        counts = [sum(v%2==e for v in missed) for e in (0,1)]
        value = sum(w*c for w,c in zip(weights,counts))
        assert result['weighted_lower']<=value<=round(solver.ObjectiveValue())
        if status==cp_model.OPTIMAL:
            assert result['weighted_lower']==value==round(solver.ObjectiveValue())
        result.update(residual_counts=counts,weighted_upper=value,
                      system=[dict(modulus=d,residue=a) for d,a in sorted(system)])
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--period',type=int,default=2310)
    ap.add_argument('--pattern',default='0111')
    ap.add_argument('--weights',nargs=2,type=int,required=True)
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    result = solve_weighted(args.period,tuple(map(int,args.pattern)),args.weights,args.time,args.workers)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='system'}),flush=True)


if __name__=='__main__':
    main()
