"""Shared exhaustive cut generation, followed by optional 165 case splitting.

Every final decision uses exact integer core scores and audited rational
tail coefficients. Numerical LP bounds only decide which search to try next.
"""
import argparse,contextlib,io,json,time
from fractions import Fraction
from pathlib import Path
from density30_branch_lp import optimize
from batch_core2310_power import build


def read(p):return json.loads(Path(p).read_text())


def score(core,G):
    w=core['weights_numerators'];den=core['weights_denominator'];free=core['initial_free_moduli']
    cap=sum(max(sum(w[x] for x in range(r,2310,d)) for r in range(d)) for d in free)
    return Fraction(cap,den)+sum(Fraction(v,core['power_seed_denominator'])*g for v,g in zip(core['power_seed_numerators'],G))


def run(rounds=6):
    started=time.time();baseline=read('data/core2310-allproj-primary-scores.json');schedule=read('data/schedule-core2310-conditional-final.json')
    p=read('data/core2310-common-tail-basis.json');a=read('data/audit-core2310-common-tail-basis.json')
    G=[max(Fraction(int(x),int(p['scale'])),Fraction(int(y),int(a['scale']))) for x,y in zip(p['upper_integers'],a['upper_integers'])]
    accepted=[r for r in baseline if Fraction(r['total_upper'])<1];pending=[r for r in baseline if Fraction(r['total_upper'])>=1]
    cuts=set();hard=[]
    for row in pending:
        data=read(read(row['core_file'])['measure_file']);cuts.update(map(tuple,data['pricing_cuts']))
    for iteration in range(rounds):
        candidates=[]
        for row in pending:
            data=read(read(row['core_file'])['measure_file']);data['pricing_cuts']=[list(c) for c in sorted(cuts)]
            key=row['projections'];name='-'.join(map(str,key));out=f'data/core2310-refined-{iteration}-{name}-measure.json'
            result=optimize(data,out,schedule,quiet=True,price=False)
            if result['numerical_restricted_lp']>=1.00002:
                hard.append(dict(projections=key,measure_file=out,restricted_lp=result['numerical_restricted_lp']));continue
            candidates.append(dict(projections=key,measure_file=out))
        source=f'data/core2310-refinement-{iteration}-input.json';manifest=f'data/core2310-refinement-{iteration}-manifest.json'
        if not candidates:pending=[];break
        Path(source).write_text(json.dumps(candidates,indent=2))
        with contextlib.redirect_stdout(io.StringIO()):build(source,manifest)
        pending=[]
        for row,entry in zip(candidates,read(manifest)):
            core=read(entry['core_file']);cuts.update(map(tuple,core['maximizing_catalogues']));total=score(core,G)
            record=dict(projections=row['projections'],core_file=entry['core_file'],total_upper=str(total),value=float(total))
            (accepted if total<1 else pending).append(record)
        Path('data/core2310-refinement-state.json').write_text(json.dumps(dict(accepted=accepted,pending=pending,hard=hard,cuts=[list(c) for c in sorted(cuts)]),indent=2))
        print('REFINEMENT',iteration,'accepted',len(accepted),'pending',len(pending),'hard',len(hard),'seconds',round(time.time()-started,1),flush=True)
        if not pending:break
    for row in pending:hard.append(dict(projections=row['projections'],measure_file=read(row['core_file'])['measure_file']))
    Path('data/core2310-refinement-state.json').write_text(json.dumps(dict(accepted=accepted,pending=[],hard=hard,cuts=[list(c) for c in sorted(cuts)]),indent=2))
    assert len(accepted)+len(hard)==1800
    print('REFINEMENT END',len(accepted),'accepted',len(hard),'need branching',flush=True)


def branch165():
    started=time.time();state=read('data/core2310-refinement-state.json');schedule=read('data/schedule-core2310-conditional-final.json');rows=[];cuts=state['cuts']
    for entry in state['hard']:
        data=read(entry['measure_file']);data['free_moduli']=[d for d in data['free_moduli'] if d!=165];data['pricing_cuts']=cuts
        for a in range(15):
            data['family']['165']=next(r for r in range(165) if r%15==a and r%11==4)
            fixed={int(d):r for d,r in data['family'].items() if int(d) not in data['free_moduli']};support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
            data.update(support=support,base=sorted({x%210 for x in support}),fixed_family=fixed)
            name='-'.join(map(str,entry['projections']));out=f'data/core2310-branch165-{name}-{a}-measure.json'
            result=optimize(data,out,schedule,quiet=True,price=False)
            rows.append(dict(projections=entry['projections'],residue165=a,measure_file=out,restricted_lp=result['numerical_restricted_lp']))
        Path('data/core2310-branch165-summary.json').write_text(json.dumps(rows,indent=2));print('BRANCH165',len(rows),'seconds',round(time.time()-started,1),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['refine','branch165']);a=ap.parse_args();run() if a.mode=='refine' else branch165()
