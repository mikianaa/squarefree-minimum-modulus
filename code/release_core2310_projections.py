"""Progressively test all projections of 66, then other small stage-11 types."""
import argparse,json
from pathlib import Path
from density30_branch_lp import optimize


def run(worst_only=False):
    data=json.loads(Path('data/core2310-density30-5-measure.json').read_text())
    schedule=json.loads(Path('data/schedule-core2310-conditional-final.json').read_text())
    target=Path('data/core2310-release66-summary.json');summary=json.loads(target.read_text()) if target.exists() else []
    completed={(r['residue66'],r['residue110']) for r in summary}
    if summary:data=json.loads(Path(summary[-1]['measure_file']).read_text())
    for a in [0,2,3,4,5,1]:
        for b in ([5] if worst_only else [5,0,1,2,3,4,6,7,8,9]):
            if (a,b) in completed:continue
            data['family']['66']=next(r for r in range(66) if r%6==a and r%11==7)
            data['family']['110']=next(r for r in range(110) if r%10==b and r%11==2)
            fixed={int(d):r for d,r in data['family'].items() if int(d) not in data['free_moduli']}
            support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
            data.update(support=support,base=sorted({x%210 for x in support}),fixed_family=fixed)
            path=f'data/core2310-release66-{a}-{b}-measure.json';data=optimize(data,path,schedule,quiet=True)
            row=dict(residue66=a,residue110=b,value=data['numerical_power_bound'],measure_file=path);summary.append(row)
            target.write_text(json.dumps(summary,indent=2));print('RELEASE66',row,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--worst-only',action='store_true');a=ap.parse_args();run(a.worst_only)
