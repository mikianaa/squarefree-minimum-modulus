"""Numerically release stage-11 conditions, reusing exhaustive pricing cuts.

Search output is not a covering certificate. The support is the complement
of the fixed subfamily; every released modulus is charged separately.
"""
import argparse,json
from pathlib import Path
from optimize_conditional_lp import optimize


def run(measure,schedule,prefix,start=8,end=15):
    data=json.loads(Path(measure).read_text())
    order=[2310,1155,770,462,385,330,231,165,154,110,77,66,55,33,22,11]
    summary=[]
    for count in range(start,end+1):
        free=sorted(order[:count]);fixed={int(d):r for d,r in data['family'].items() if int(d) not in free}
        support=[x for x in range(2310) if all(x%d!=r for d,r in fixed.items())]
        data.update(support=support,base=sorted(set(x%210 for x in support)),fixed_family=fixed,free_moduli=free,support_mode='fixed_subfamily')
        source=f'{prefix}-{count}-input.json';out=f'{prefix}-{count}-measure.json'
        Path(source).write_text(json.dumps(data,separators=(',',':')))
        data=optimize(source,out,40,schedule,free,True)
        row=dict(free_count=count,residual=len(support),value=data['numerical_conditional_bound'],measure_file=out)
        summary.append(row);print('BRANCH',row,flush=True)
        Path(f'{prefix}-summary.json').write_text(json.dumps(summary,indent=2))
        if row['value']>1.03:break


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('measure');ap.add_argument('schedule');ap.add_argument('prefix');ap.add_argument('--start',type=int,default=8);ap.add_argument('--end',type=int,default=15)
    a=ap.parse_args();run(a.measure,a.schedule,a.prefix,a.start,a.end)
