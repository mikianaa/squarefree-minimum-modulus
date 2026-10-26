"""Finite-pool candidates for all 1800 small projections of a fixed core.

Every resulting candidate must be priced and audited by exact batch tools.
No finite-pool objective is used as an upper bound or proof.
"""
import argparse,itertools,json,time
from pathlib import Path
from density30_branch_lp import optimize


def run(limit=None):
    started=time.time();data=json.loads(Path('data/core2310-release66-1-9-measure.json').read_text())
    schedule=json.loads(Path('data/schedule-core2310-conditional-final.json').read_text())
    output=Path('data/core2310-all-projections-summary.json');rows=json.loads(output.read_text()) if output.exists() else []
    done={tuple(r['projections']) for r in rows};count=0
    for key in itertools.product(range(2),range(3),range(5),range(6),range(10)):
        if key in done:continue
        for d,a,colour in zip([22,33,55,66,110],key,[3,10,6,7,2]):
            data['family'][str(d)]=next(r for r in range(d) if r%(d//11)==a and r%11==colour)
        fixed={int(d):r for d,r in data['family'].items() if int(d) not in data['free_moduli']}
        support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
        data.update(support=support,base=sorted({x%210 for x in support}),fixed_family=fixed)
        name='-'.join(map(str,key));path=f'data/core2310-allproj-{name}-measure.json'
        result=optimize(data,path,schedule,quiet=True,price=False)
        rows.append(dict(projections=list(key),restricted_lp=result['numerical_restricted_lp'],measure_file=path));count+=1
        if count%30==0:
            output.write_text(json.dumps(rows,indent=2));print('CANDIDATES',len(rows),'seconds',round(time.time()-started,1),flush=True)
        if limit and count>=limit:break
    output.write_text(json.dumps(rows,indent=2));print('CANDIDATE SEARCH END',len(rows),'seconds',round(time.time()-started,1),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int);a=ap.parse_args();run(a.limit)
