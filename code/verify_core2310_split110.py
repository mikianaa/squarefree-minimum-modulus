"""Build integer case certificates, then audit the complete modulus-110 split.

The ten fresh-colour cases dominate every colour by uniform-fibre spreading.
This eliminates the 110-residue condition but retains the other fixed classes.
"""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path


def read(p):return json.loads(Path(p).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def build():
    from core2310_conditional_envelope import run
    from audit_core2310_power_envelope import audit
    for a in range(10):
        source=f'data/core2310-split110-{a}-measure.json';core=f'data/core2310-split110-{a}-envelope.json';check=f'data/audit-core2310-split110-{a}-envelope.json'
        run(source,core);audit(core,check)
        print('INTEGER CASE DONE',a,flush=True)


def audit_split():
    primary=read('data/core2310-common-tail-basis.json');secondary=read('data/audit-core2310-common-tail-basis.json')
    assert primary['status']=='CERTIFIED' and secondary['status']=='VERIFIED'
    assert primary['loads']==secondary['loads']==[1,2,4,8,16,32]
    assert primary['schedule_sha256']==secondary['schedule_sha256']==digest(primary['schedule_file'])
    G=[max(Fraction(int(x),int(primary['scale'])),Fraction(int(y),int(secondary['scale']))) for x,y in zip(primary['upper_integers'],secondary['upper_integers'])]
    expected_free=[154,165,231,330,385,462,770,1155,2310];outer=None;rows=[]
    for a in range(10):
        corefile=f'data/core2310-split110-{a}-envelope.json';core=read(corefile);check=read(f'data/audit-core2310-split110-{a}-envelope.json');measure=read(core['measure_file'])
        assert check['status']=='VERIFIED' and check['core_sha256']==digest(corefile)
        assert measure['free_moduli']==core['initial_free_moduli']==expected_free
        family={int(d):r for d,r in measure['family'].items()};fixed={d:r for d,r in family.items() if d not in expected_free}
        assert set(family)=={d for d in range(7,2311) if 2310%d==0}
        assert family[110]%10==a and family[110]%11==2
        current={d:r for d,r in fixed.items() if d!=110}
        if outer is None:outer=current
        else:assert current==outer
        baseline={x for x in range(2310) if all(x%d!=r for d,r in outer.items())}
        reference={x for x in baseline if x%110!=family[110]}
        assert reference==set(measure['support'])
        weights=core['weights_numerators'];den=core['weights_denominator'];base=sorted({x%210 for x in baseline})
        one={};sizes={}
        for x in base:
            fibre=[y for y in reference if y%210==x];sizes[x]=len(fibre);assert fibre
            values={weights[y] for y in fibre};assert len(values)==1
            one[x]=values.pop()
            # All reference stage-11 colours are distinct. This is the
            # smallest possible fibre for these base projections, whatever
            # the actual colours (or any omitted classes) may be.
            active=sum(x%(d//11)==r%(d//11) for d,r in fixed.items() if d%11==0)
            assert sizes[x]==11-active
        colours=[r%11 for d,r in fixed.items() if d%11==0]
        assert len(set(colours))==len(colours)<11
        assert all(weights[y]==0 for y in range(2310) if y not in reference)
        # Every other 11-colour removes at most the fresh-colour removal.
        for residue in range(a,110,10):
            actual={x for x in baseline if x%110!=residue}
            assert all(sum(y%210==x for y in actual)>=sizes[x] for x in base)
        caps={}
        for d in expected_free:
            full=max(sum(weights[y] for y in range(r,2310,d)) for r in range(d))
            projected=max(sum(one[x] for x in base if x%(d//11)==r) for r in range(d//11))
            # Each actual point gets q_x/s_actual <= q_x/s_fresh = one[x]/den.
            assert full==projected
            caps[d]=full
        tail=sum(Fraction(w,core['power_seed_denominator'])*g for w,g in zip(core['power_seed_numerators'],G))
        total=tail+Fraction(sum(caps.values()),den);assert total<1
        rows.append(dict(base_residue=a,total_upper=str(total),total_upper_float=float(total),core_file=corefile,core_sha256=digest(corefile)))
    maximum=max(Fraction(r['total_upper']) for r in rows)
    result=dict(status='VERIFIED',claim='No distinct squarefree covering system with all moduli >6 can contain the listed fixed subfamily',
                fixed_family=outer,free_moduli=sorted(expected_free+[110]),case_count=10,covered_mod110_residues=110,omission_allowed=True,
                maximum_upper=str(maximum),maximum_upper_float=float(maximum),cases=rows,
                scope='One core-210 family and six fixed stage-11 classes; all other stage-11 residues and later primes are unrestricted. The global bound remains 7.')
    result['stronger_projection_claim']='No cover for the fixed core-210 family when each listed stage-11 modulus, if present, has the listed base projection; every 11-colour and omission is allowed'
    result['fixed_core210_family']={d:r for d,r in outer.items() if d%11!=0}
    result['stage11_projection_conditions']={d:dict(modulus=d//11,residue=r%(d//11)) for d,r in outer.items() if d%11==0 and d!=11}
    result['scope']='One core-210 family and five base-projection conditions; all 11-colours, all other stage-11 choices and all later primes are unrestricted. Global upper bound remains 7.'
    Path('data/audit-core2310-split110-branch.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','audit']);args=ap.parse_args()
    build() if args.mode=='build' else audit_split()
