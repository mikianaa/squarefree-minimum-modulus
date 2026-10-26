"""Combine exact unrestricted-class charges with independently audited tail.

Every future distortion preserves the core marginal. Consequently arbitrary
core-class residues cost at most their maximum Q-mass. Their sum plus the
later-prime loss strictly below one excludes the specified fixed subfamily.
"""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(corefile,schedule,coreaudit,primary,secondary,output):
    core=read(corefile);measure=read(core['measure_file']);a=read(primary);b=read(secondary);geometry=read(coreaudit)
    assert a['status']=='CERTIFIED' and b['status']=='VERIFIED' and geometry['status']=='VERIFIED'
    assert a['core_sha256']==b['core_sha256']==geometry['core_sha256']==digest(corefile)
    assert a['schedule_sha256']==b['schedule_sha256']==digest(schedule)
    free=measure['free_moduli'];family={int(d):r for d,r in measure['family'].items()};fixed={d:r for d,r in family.items() if d not in free}
    assert set(family)=={d for d in range(7,2311) if 2310%d==0}
    assert set(free)<=set(family) and len(set(free))==len(free)
    if 'initial_free_moduli' in core:assert core['initial_free_moduli']==free
    w=core['weights_numerators'];den=core['weights_denominator'];assert len(w)==2310 and sum(w)==den and den>0
    assert all(isinstance(v,int) and v>=0 for v in w)
    assert all(v==0 or all(x%d!=r for d,r in fixed.items()) for x,v in enumerate(w))
    caps={d:max(sum(w[x] for x in range(r,2310,d)) for r in range(d)) for d in free}
    extra=Fraction(sum(caps.values()),den)
    tail=max(Fraction(int(a['total_integer']),int(a['scale'])),Fraction(int(b['total_integer']),int(b['scale'])))
    total=tail+extra
    result=dict(status='VERIFIED' if total<1 else 'NOT_VERIFIED',claim='No distinct squarefree covering system with all moduli >6 can contain the specified fixed subfamily' if total<1 else 'This bound is insufficient',
                fixed_family=fixed,free_moduli=free,free_class_capacity_numerators=caps,capacity_denominator=den,
                free_charge=str(extra),tail_upper=str(tail),total_upper=str(total),total_upper_float=float(total),
                core_file=corefile,core_sha256=digest(corefile),schedule_file=schedule,schedule_sha256=digest(schedule),
                scope='All choices or omissions of the listed free moduli and all later primes; other initial configurations remain unresolved')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser()
    for name in ['corefile','schedule','coreaudit','primary','secondary','output']:ap.add_argument(name)
    audit(**vars(ap.parse_args()))
