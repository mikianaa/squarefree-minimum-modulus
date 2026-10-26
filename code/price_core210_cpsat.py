"""CP-SAT pricing of a single core stop-loss envelope; solver bounds are diagnostic.

This does not emit a proof log and is not an independent proof certificate.
Weights are nonnegative integers, normalized by their actual sum.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from core210_joint_probe import DIVS


def solve(measure_file,threshold,output,seconds=60):
    data=json.loads(Path(measure_file).read_text());points=data['support']
    weights=[max(0,round(w*1000000)) for w in data['weights']];den=sum(weights)
    model=cp_model.CpModel();z={(d,r):model.NewBoolVar(f'z{d}_{r}') for d in DIVS[1:] for r in range(d)}
    for d in DIVS[1:]:model.AddExactlyOne(z[d,r] for r in range(d))
    counts=[];losses={}
    for x in range(210):
        count=model.NewIntVar(1,16,f'n{x}');model.Add(count==1+sum(z[d,x%d] for d in DIVS[1:]));counts.append(count)
    # Valid global concentration cuts from the uniform product comparison.
    full_law={1:48,2:92,4:56,8:13,16:1}
    for t in sorted(set([1,2,4,8,threshold])):
        ls=[]
        for x in range(210):
            loss=model.NewIntVar(0,16-t,f'h{t}_{x}');model.AddMaxEquality(loss,[0,counts[x]-t]);ls.append(loss)
            if t==threshold:losses[x]=loss
        model.Add(sum(ls)<=sum(n*max(0,f-t) for f,n in full_law.items()))
    model.Maximize(sum(w*losses[x] for x,w in zip(points,weights)))
    pool=[a for a in data['pool'] if a['threshold']==threshold]
    values=[]
    for a in pool:
        counts_hint=[sum(x%d==r for d,r in zip(DIVS,a['residues'])) for x in points]
        values.append(sum(w*max(0,n-threshold) for w,n in zip(weights,counts_hint)))
    best=pool[int(np.argmax(values))]['residues']
    for i,d in enumerate(DIVS[1:],1):
        for r in range(d):model.AddHint(z[d,r],int(r==best[i]))
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=seconds;solver.parameters.num_workers=4
    start=time.time();status=solver.Solve(model)
    result=dict(threshold=threshold,status=solver.StatusName(status),denominator=den,
                solver_upper=solver.BestObjectiveBound(),normalized_solver_upper=solver.BestObjectiveBound()/den,
                seconds=time.time()-start,scope='Solver bound only, not an audited proof log')
    if status in [cp_model.OPTIMAL,cp_model.FEASIBLE]:
        residues=[0]+[next(r for r in range(d) if solver.Value(z[d,r])) for d in DIVS[1:]]
        counts_actual=[sum(x%d==r for d,r in zip(DIVS,residues)) for x in points]
        score=sum(w*max(0,n-threshold) for w,n in zip(weights,counts_actual))
        assert abs(score-solver.ObjectiveValue())<.1
        result.update(residues=residues,integer_witness=score,normalized_witness=score/den)
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('threshold',type=int);ap.add_argument('output');ap.add_argument('--seconds',type=int,default=60)
    a=ap.parse_args();solve(a.measure,a.threshold,a.output,a.seconds)
