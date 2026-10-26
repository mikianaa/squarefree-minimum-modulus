"""Integer audit of all 70 fresh-colour projection pairs (77,110).

Case search can run concurrently with build; only completed summary entries
are consumed. A final audit requires exactly all 70 cases, each below one.
"""
import argparse,contextlib,hashlib,io,json,time
from fractions import Fraction
from pathlib import Path


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build():
    from core2310_conditional_envelope import run
    from audit_core2310_power_envelope import audit
    for a in range(7):
        for b in range(10):
            while True:
                try:summary=read('data/core2310-split77-110-summary.json')
                except (FileNotFoundError,json.JSONDecodeError):time.sleep(1);continue
                rows=[r for r in summary if r['base_residue77']==a and r['base_residue110']==b]
                if rows:break
                time.sleep(1)
            row=rows[0]
            if row['value']>=1:
                print('NOT YET BELOW ONE',a,b,row['value'],flush=True);continue
            core=f'data/core2310-split77-110-{a}-{b}-envelope.json';check=f'data/audit-core2310-split77-110-{a}-{b}-envelope.json'
            if Path(check).exists() and Path(core).exists() and read(check).get('core_sha256')==digest(core):
                print('ALREADY VERIFIED',a,b,flush=True);continue
            with contextlib.redirect_stdout(io.StringIO()):run(row['measure_file'],core);audit(core,check)
            print('INTEGER PAIR VERIFIED',a,b,flush=True)


def audit():
    primary=read('data/core2310-common-tail-basis.json');secondary=read('data/audit-core2310-common-tail-basis.json')
    assert primary['status']=='CERTIFIED' and secondary['status']=='VERIFIED'
    assert primary['loads']==secondary['loads']==[1,2,4,8,16,32]
    assert primary['schedule_sha256']==secondary['schedule_sha256']==digest(primary['schedule_file'])
    G=[max(Fraction(int(x),int(primary['scale'])),Fraction(int(y),int(secondary['scale']))) for x,y in zip(primary['upper_integers'],secondary['upper_integers'])]
    free=[154,165,231,330,385,462,770,1155,2310];basefamily=None;conditions=None;rows=[]
    for a in range(7):
        for b in range(10):
            corefile=f'data/core2310-split77-110-{a}-{b}-envelope.json';core=read(corefile);check=read(f'data/audit-core2310-split77-110-{a}-{b}-envelope.json');measure=read(core['measure_file'])
            assert check['status']=='VERIFIED' and check['core_sha256']==digest(corefile)
            assert core['initial_free_moduli']==measure['free_moduli']==free
            family={int(d):r for d,r in measure['family'].items()};assert set(family)=={d for d in range(7,2311) if 2310%d==0}
            assert family[77]%7==a and family[110]%10==b
            fixed={d:r for d,r in family.items() if d not in free}
            early={d:r for d,r in fixed.items() if d%11!=0};stage={d:r for d,r in fixed.items() if d%11==0}
            cond={d:(d//11,r%(d//11)) for d,r in stage.items() if d not in [11,77,110]}
            if basefamily is None:basefamily=early;conditions=cond
            else:assert early==basefamily and cond==conditions
            colours=[r%11 for r in stage.values()];assert len(set(colours))==len(colours)<11
            base=[x for x in range(210) if all(x%d!=r for d,r in early.items())]
            reference={x for x in range(2310) if all(x%d!=r for d,r in fixed.items())};assert reference==set(measure['support'])
            weights=core['weights_numerators'];den=core['weights_denominator'];one={}
            for x in base:
                lifts=[y for y in reference if y%210==x];active=sum(x%(d//11)==r%(d//11) for d,r in stage.items())
                assert len(lifts)==11-active>0
                values={weights[y] for y in lifts};assert len(values)==1;one[x]=values.pop()
            assert all(weights[y]==0 for y in range(2310) if y not in reference)
            caps=[]
            for d in free:
                direct=max(sum(weights[y] for y in range(r,2310,d)) for r in range(d))
                projection=max(sum(one[x] for x in base if x%(d//11)==r) for r in range(d//11));assert direct==projection;caps.append(direct)
            tail=sum(Fraction(w,core['power_seed_denominator'])*g for w,g in zip(core['power_seed_numerators'],G))
            total=tail+Fraction(sum(caps),den);assert total<1
            rows.append(dict(residue77=a,residue110=b,total_upper=str(total),total_upper_float=float(total),core_sha256=digest(corefile)))
    maximum=max(Fraction(r['total_upper']) for r in rows)
    result=dict(status='VERIFIED',claim='No distinct squarefree cover with all moduli >6, the fixed core-210 family, and these four stage-11 projection conditions',
                fixed_core210_family=basefamily,stage11_projection_conditions={d:dict(modulus=m,residue=r) for d,(m,r) in conditions.items()},
                fully_free_stage11_moduli=sorted(free+[11,77,110]),all_11_colours_free=True,omissions_allowed=True,
                case_count=len(rows),maximum_upper=str(maximum),maximum_upper_float=float(maximum),cases=rows,
                scope='Other core-210 families and other base projections at 22,33,55,66 remain unresolved; global upper bound remains 7')
    Path('data/audit-core2310-split77-110-branch.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','audit']);a=ap.parse_args();build() if a.mode=='build' else audit()
