"""Combine an audited tail with exact union charges for six unrestricted classes.

Later distortion steps preserve the core marginal. Therefore a class of
core modulus d costs at most max_a Q(x=a mod d), regardless of all future
choices. The reference measure may have artificial zeroes at other points;
only support outside the FIXED subfamily is necessary.
"""
import hashlib,json
from fractions import Fraction
from pathlib import Path


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit():
    corefile='data/core2310-six-free-envelope.json';schedule='data/schedule-core2310-six-free.json'
    core=read(corefile);measure=read(core['measure_file']);a=read('data/cert-core2310-six-free-tail.json');b=read('data/audit-core2310-six-free-tail.json')
    geometry=read('data/audit-core2310-six-free-envelope.json')
    assert a['status']=='CERTIFIED' and b['status']=='VERIFIED' and geometry['status']=='VERIFIED'
    assert a['core_sha256']==b['core_sha256']==geometry['core_sha256']==digest(corefile)
    assert a['schedule_sha256']==b['schedule_sha256']==digest(schedule)
    free=[330,385,462,770,1155,2310];assert measure['free_moduli']==free
    family={int(d):r for d,r in measure['family'].items()};fixed={d:r for d,r in family.items() if d not in free}
    assert set(fixed)|set(free)=={d for d in range(7,2311) if 2310%d==0}
    w=core['weights_numerators'];den=core['weights_denominator'];assert len(w)==2310 and sum(w)==den
    assert all(v==0 or all(x%d!=r for d,r in fixed.items()) for x,v in enumerate(w))
    caps={d:max(sum(w[x] for x in range(r,2310,d)) for r in range(d)) for d in free}
    extra=Fraction(sum(caps.values()),den)
    tail=max(Fraction(int(a['total_integer']),int(a['scale'])),Fraction(int(b['total_integer']),int(b['scale'])))
    total=tail+extra;assert total<1
    result=dict(status='VERIFIED',claim='No distinct squarefree covering system with all moduli >6 can contain the specified fixed subfamily',
                fixed_family=fixed,free_moduli=free,free_class_capacity_numerators=caps,capacity_denominator=den,
                free_charge=str(extra),tail_upper=str(tail),total_upper=str(total),total_upper_float=float(total),
                core_file=corefile,core_sha256=digest(corefile),schedule_file=schedule,schedule_sha256=digest(schedule),
                scope='All choices or omissions of the six listed moduli and all later primes; remaining early configurations are not excluded, so global upper bound stays 7')
    Path('data/audit-core2310-six-free-branch.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit()
