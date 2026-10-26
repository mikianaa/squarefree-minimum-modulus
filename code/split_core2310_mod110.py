"""Numerical ten-case fresh-colour split for modulus 110.

Uniform weights on each surviving 11-fibre enable a separate domination
argument for existing colours. Exhaustive integer audits are still needed.
"""
import json
from pathlib import Path
from optimize_conditional_lp import optimize


def run():
    data=json.loads(Path('data/core2310-expanded-9-measure.json').read_text());summary=[]
    free=[154,165,231,330,385,462,770,1155,2310]
    for a in range(10):
        residue=next(r for r in range(110) if r%10==a and r%11==2)
        data['family']['110']=residue
        fixed={int(d):r for d,r in data['family'].items() if int(d) not in free}
        support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
        data.update(support=support,base=sorted(set(x%210 for x in support)),fixed_family=fixed,free_moduli=free,support_mode='fixed_subfamily')
        source=f'data/core2310-split110-{a}-input.json';output=f'data/core2310-split110-{a}-measure.json'
        Path(source).write_text(json.dumps(data,separators=(',',':')))
        data=optimize(source,output,40,'data/schedule-core2310-conditional-final.json',free,False)
        row=dict(base_residue=a,mod110_residue=residue,value=data['numerical_conditional_bound'],measure_file=output,residual=len(support))
        summary.append(row);print('SPLIT110',row,flush=True)
        Path('data/core2310-split110-summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':run()
