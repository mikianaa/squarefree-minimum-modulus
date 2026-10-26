"""Ten certificates with a conditional permutation covering every 77 projection.

The reference 7-plane 6 is complete on the core-30 residual. Point density
depends only on x mod 30. Swapping 7-planes within each 30-fibre preserves
the two-level comparison and bounds all free-class capacities.
"""
import argparse,contextlib,hashlib,io,json,itertools
from fractions import Fraction
from pathlib import Path


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def search():
    from optimize_conditional_lp import optimize
    cuts=read('data/core2310-density30-probe.json')['pricing_cuts'];rows=[]
    for b in range(10):
        data=read(f'data/core2310-split110-{b}-measure.json');data['pricing_cuts']=cuts
        source=f'data/core2310-density30-{b}-input.json';out=f'data/core2310-density30-{b}-measure.json'
        Path(source).write_text(json.dumps(data,separators=(',',':')))
        result=optimize(source,out,25,'data/schedule-core2310-conditional-final.json',data['free_moduli'],True,30)
        cuts=result['pricing_cuts'];rows.append(dict(residue110=b,value=result['numerical_conditional_bound'],measure_file=out))
        Path('data/core2310-density30-summary.json').write_text(json.dumps(rows,indent=2));print('DENSITY30',rows[-1],flush=True)


def build():
    from core2310_conditional_envelope import run
    from audit_core2310_power_envelope import audit
    for b in range(10):
        with contextlib.redirect_stdout(io.StringIO()):
            core=f'data/core2310-density30-{b}-envelope.json'
            run(f'data/core2310-density30-{b}-measure.json',core)
            audit(core,f'data/audit-core2310-density30-{b}-envelope.json')
        print('DENSITY30 INTEGER VERIFIED',b,flush=True)


def audit(case_files=None,released=(77,110),output='data/audit-core2310-density30-transport.json',adaptive165=False):
    first=read('data/core2310-common-tail-basis.json');second=read('data/audit-core2310-common-tail-basis.json')
    assert first['status']=='CERTIFIED' and second['status']=='VERIFIED'
    assert first['loads']==second['loads']==[1,2,4,8,16,32]
    assert first['schedule_sha256']==second['schedule_sha256']==digest(first['schedule_file'])
    G=[max(Fraction(int(x),int(first['scale'])),Fraction(int(y),int(second['scale']))) for x,y in zip(first['upper_integers'],second['upper_integers'])]
    base_free=[154,165,231,330,385,462,770,1155,2310];early_reference=None;conditions_reference=None;rows=[]
    if case_files is None:case_files=[f'data/core2310-density30-{b}-envelope.json' for b in range(10)]
    covered=set();branches={};split_moduli=sorted(d for d in released if d!=77 and (not adaptive165 or d!=165))
    assert 77 in released and all(d in [22,33,55,66,77,110,165] for d in released)
    if adaptive165:assert set(released)=={22,33,55,66,77,110,165}
    for corefile in case_files:
        core=read(corefile);measure=read(core['measure_file']);cp=Path(corefile);check=read(cp.with_name('audit-'+cp.name))
        assert check['status']=='VERIFIED' and check['core_sha256']==digest(corefile)
        free=core['initial_free_moduli'];assert measure['free_moduli']==free
        assert free==base_free or (adaptive165 and free==[d for d in base_free if d!=165])
        family={int(d):r for d,r in measure['family'].items()};assert set(family)=={d for d in range(7,2311) if 2310%d==0}
        early={d:r for d,r in family.items() if d%11!=0};stage={d:r for d,r in family.items() if d%11==0 and d not in free}
        assert stage[77]%7==6;b=stage[110]%10
        key=tuple(stage[d]%(d//11) for d in split_moduli)
        conditioned165=stage[165]%15 if 165 in stage else None
        if adaptive165:
            branches.setdefault(key,set());assert conditioned165 not in branches[key];branches[key].add(conditioned165)
        else:assert key not in covered
        covered.add(key)
        conditions={d:dict(modulus=d//11,residue=r%(d//11)) for d,r in stage.items() if d not in [11,*released]}
        if early_reference is None:early_reference=early;conditions_reference=conditions
        else:assert early==early_reference and conditions==conditions_reference
        colours=[r%11 for r in stage.values()];assert len(set(colours))==len(colours)<11
        U={x for x in range(210) if all(x%d!=r for d,r in early.items())};R=sorted({x%30 for x in U})
        reference={x for x in range(2310) if x%210 in U and all(x%d!=r for d,r in stage.items())}
        assert reference==set(measure['support'])
        weights=core['weights_numerators'];den=core['weights_denominator'];density={};baseq={};sizes={}
        for x in R:
            indices=[y for y in reference if y%30==x];vals={weights[y] for y in indices};assert len(vals)==1;density[x]=vals.pop()
            assert any(y%7==6 for y in U if y%30==x)
        for y in U:
            sizes[y]=sum(z%210==y for z in reference)
            active=sum(y%(d//11)==r%(d//11) for d,r in stage.items());assert sizes[y]==11-active>0
            baseq[y]=sizes[y]*density[y%30]
        assert sum(baseq.values())==den and all(weights[y]==0 for y in range(2310) if y not in reference)
        # A complete reference plane exists. For any target plane a, swap
        # it with plane 6 when both survive, independently within each x30.
        # Sorted base masses, hence the second reduction's sums/maxima,
        # stay identical; each actual point density is at most density[x30].
        for a in range(7):
            transported={}
            for x in R:
                by7={y%7:y for y in U if y%30==x};assert 6 in by7
                for colour,y in by7.items():
                    source=y
                    if a in by7:
                        if colour==a:source=by7[6]
                        elif colour==6:source=by7[a]
                    transported[y]=baseq[source]
                assert sorted(transported[y] for y in by7.values())==sorted(baseq[y] for y in by7.values())
            for y,q in transported.items():
                active=sum(y%(d//11)==r%(d//11) for d,r in stage.items() if d!=77)+(y%7==a)
                actual_size=11-active
                assert density[y%30]<=q<=actual_size*density[y%30]
            assert sum(transported.values())==den
        caps=[]
        for d in free:
            m=d//11;direct=max(sum(weights[y] for y in range(r,2310,d)) for r in range(d))
            if m%7==0:
                projection=max(sum(density[x] for x in R if x%(m//7)==r) for r in range(m//7))
            else:
                projection=max(sum(density[y%30] for y in U if y%m==r) for r in range(m))
            assert direct==projection;caps.append(direct)
        tail=sum(Fraction(w,core['power_seed_denominator'])*g for w,g in zip(core['power_seed_numerators'],G))
        total=tail+Fraction(sum(caps),den);assert total<1
        rows.append(dict(projections=dict(zip(split_moduli,key)),conditioned165=conditioned165,residue110=b,total_upper=str(total),total_upper_float=float(total),core_sha256=digest(corefile)))
    assert covered==set(itertools.product(*(range(d//11) for d in split_moduli)))
    if adaptive165:assert all(values=={None} or values==set(range(15)) for values in branches.values())
    maximum=max(Fraction(r['total_upper']) for r in rows)
    result=dict(status='VERIFIED',claim='No distinct squarefree cover with all moduli >6 under the fixed core-210 family and the listed stage-11 projections',
                fixed_core210_family=early_reference,stage11_projection_conditions=conditions_reference,all_11_colours_free=True,
                fully_free_stage11_moduli=sorted(set(base_free+[11,*released])),omissions_allowed=True,case_count=len(rows),base_projection_cases=len(covered),transported_projection_pairs=7*len(rows),
                maximum_upper=str(maximum),maximum_upper_float=float(maximum),cases=rows,
                scope='Only this core-210 family and these projections; all other initial cases remain open, so the global upper bound is still 7')
    Path(output).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['search','build','audit']);args=ap.parse_args()
    {'search':search,'build':build,'audit':audit}[args.mode]()
