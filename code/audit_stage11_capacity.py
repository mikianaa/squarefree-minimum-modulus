"""Independent integer audit of the minimum stage-11 residual for ONE core.

Union capacity is a lower bound for every stage-11 extension of this fixed
core. An explicit extension attains it. No other core-210 family is excluded.
"""
import json
from pathlib import Path


def audit():
    old=json.loads(Path('data/core210-hard-prefix.json').read_text())
    extension=json.loads(Path('data/core2310-overlap-prefix.json').read_text())
    family={int(d):int(r) for d,r in old['family'].items()}
    base={x for x in range(210) if all(x%d!=r for d,r in family.items())}
    assert base==set(old['support']) and len(base)==104
    ext={int(d):int(r) for d,r in extension['extension'].items()}
    divisors=[d for d in range(1,211) if 210%d==0]
    assert set(ext)=={11*d for d in divisors}
    caps={d:max(sum(x%d==r for x in base) for r in range(d)) for d in divisors if d>1}
    lower=10*len(base)-sum(caps.values())
    # A translation by a multiple of 210 can always normalize the 11-class.
    domain={x for x in range(2310) if x%210 in base and x%11!=ext[11]}
    seen=set();sizes={}
    for d in divisors[1:]:
        covered={x for x in domain if x%(11*d)==ext[11*d]}
        assert len(covered)<=caps[d]
        assert not seen.intersection(covered)
        assert len(covered)==caps[d]
        seen.update(covered);sizes[d]=len(covered)
    residual=domain-seen
    assert residual==set(extension['support']) and len(residual)==lower==787
    result=dict(status='VERIFIED',claim='Minimum residual modulo 2310 after adding all stage-11 moduli to the specified fixed core-210 family is 787',
                base_family=family,base_uncovered=len(base),capacity_by_divisor=caps,
                total_capacity=sum(caps.values()),minimum_residual=lower,
                witness_file='data/core2310-overlap-prefix.json',
                scope='Conditional on one fixed core-210 family; neither all core-210 configurations nor the upper bound 6 is established')
    Path('data/audit-stage11-capacity.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit()
