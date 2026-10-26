"""Probe all 70 base-projection pairs, using distinct fresh 11-colours.

The uniform-fibre domination argument covers arbitrary colours. Scores are
numerical until individual integer cores and a shared tail basis are audited.
"""
import json
from pathlib import Path
from optimize_conditional_lp import optimize


def run():
    data=json.loads(Path('data/core2310-split110-9-measure.json').read_text());summary=[]
    checkpoint=Path('data/core2310-split77-110-summary.json')
    if checkpoint.exists():
        summary=json.loads(checkpoint.read_text())
        if summary:data=json.loads(Path(summary[-1]['measure_file']).read_text())
    completed={(r['base_residue77'],r['base_residue110']) for r in summary}
    free=[154,165,231,330,385,462,770,1155,2310]
    for a in range(7):
        for b in range(10):
            if (a,b) in completed:continue
            data['family']['77']=next(r for r in range(77) if r%7==a and r%11==8)
            data['family']['110']=next(r for r in range(110) if r%10==b and r%11==2)
            fixed={int(d):r for d,r in data['family'].items() if int(d) not in free}
            support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
            data.update(support=support,base=sorted(set(x%210 for x in support)),fixed_family=fixed,free_moduli=free,support_mode='fixed_subfamily')
            source=f'data/core2310-split77-110-{a}-{b}-input.json';output=f'data/core2310-split77-110-{a}-{b}-measure.json'
            Path(source).write_text(json.dumps(data,separators=(',',':')))
            data=optimize(source,output,30,'data/schedule-core2310-conditional-final.json',free,False)
            row=dict(base_residue77=a,base_residue110=b,value=data['numerical_conditional_bound'],measure_file=output,residual=len(support))
            summary.append(row);print('SPLIT77-110',row,flush=True)
            Path('data/core2310-split77-110-summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':run()
