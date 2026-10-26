"""Numerical solver probe for the robust capped free-class catalogue.

At a base point x, at most min(s_x,n_x) fresh lifts can be covered. This
bound covers arbitrary 11-colours. Solver results here are diagnostics,
not independently checked proof certificates.
"""
import json,time
from pathlib import Path
from ortools.sat.python import cp_model


def run():
    t=time.time();core=json.loads(Path('data/core2310-hard-projection-fibre-envelope.json').read_text());data=json.loads(Path(core['measure_file']).read_text())
    support=set(data['support']);base=data['base'];weights=core['weights_numerators'];den=core['weights_denominator'];w={};sizes={}
    for x in base:
        fibre=[y for y in support if y%210==x];values={weights[y] for y in fibre};assert len(values)==1;w[x]=values.pop();sizes[x]=len(fibre)
    mods=[d//11 for d in data['free_moduli']];model=cp_model.CpModel();choices={};linear=0
    for m in mods:
        residues=sorted({x%m for x in base});choices[m]={r:model.NewBoolVar(f'c{m}_{r}') for r in residues};model.AddExactlyOne(list(choices[m].values()))
        linear+=max(sum(w[x] for x in base if x%m==r) for r in residues)
    terms=[]
    for x in base:
        n=model.NewIntVar(0,len(mods),f'n{x}');model.Add(n==sum(choices[m][x%m] for m in mods))
        hit=model.NewIntVar(0,min(sizes[x],len(mods)),f'h{x}');model.AddMinEquality(hit,[n,sizes[x]]);terms.append(w[x]*hit)
    model.Maximize(sum(terms));solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=45;solver.parameters.num_search_workers=4
    status=solver.Solve(model);result=dict(status=solver.StatusName(status),incumbent=solver.ObjectiveValue(),upper=solver.BestObjectiveBound(),
          linear_bound=linear,denominator=den,possible_saving=(linear-solver.BestObjectiveBound())/den,seconds=time.time()-t,
          scope='CP-SAT diagnostic only; not an independently audited upper or impossibility certificate')
    if status in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        residues={m:next(r for r,v in choices[m].items() if solver.Value(v)) for m in mods}
        result['base_residue_witness']=residues
        active={m:{x for x in base if x%m==residues[m] and w[x]>0} for m in mods}
        crt={(y%210,y%11):y for y in range(2310)}
        allowed={m:{c for c in range(11) if all(crt[x,c] in support for x in active[m])} for m in mods}
        graph={m:{n for n in mods if n!=m and active[m]&active[n]} for m in mods}
        def colour(assigned):
            if len(assigned)==len(mods):return assigned
            opts={m:allowed[m]-{assigned[n] for n in graph[m] if n in assigned} for m in mods if m not in assigned}
            m=min(opts,key=lambda m:(len(opts[m]),-len(graph[m])))
            for c in sorted(opts[m]):
                found=colour({**assigned,m:c})
                if found is not None:return found
            return None
        colours=colour({})
        if colours is not None:
            result['actual_residue_witness']={11*m:next(r for r in range(11*m) if r%m==residues[m] and r%11==colours[m]) for m in mods}
    Path('data/core2310-capped-free-union-probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':run()
