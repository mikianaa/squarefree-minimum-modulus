"""Independently count an actual family by ordinary integer progressions."""
import hashlib
import json
import math
from pathlib import Path
from extremal_slice_capacity import extension_bound,exclusion_threshold


def main():
    path = Path('data/min-uncovered-9699690-witness.json')
    witness = json.loads(path.read_text())
    primes = [2,3,5,7,11,13,17,19]
    period = math.prod(primes)
    assert witness['N']==period
    family = [(r['modulus'],r['residue']) for r in witness['system']]
    assert len(set(d for d,a in family))==len(family)==254
    covered = bytearray(period)
    for d,a in family:
        assert d>=3 and period%d==0 and 0<=a<d
        for n in range(a,period,d):
            covered[n] = 1
    uncovered = covered.count(0)
    assert uncovered==witness['uncovered']
    first = extension_bound(primes,uncovered,23)['uncovered_lower']
    second = extension_bound(primes+[23],first,29)['uncovered_lower']
    out = dict(verified=True,period=period,classes=len(family),uncovered=uncovered,
               parity_uncovered=[covered[e::2].count(0) for e in (0,1)],
               sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
               cardinality_method_ceiling=[dict(new_prime=23,maximum_certifiable_lower=first,
                                               positive_threshold=exclusion_threshold(primes,23)),
                                           dict(new_prime=29,maximum_certifiable_lower=second)],
               scope='The witness is an upper bound. Method ceilings are NOT lower bounds for actual families. They limit only sequential use of the sharp unconstrained size-only cylinder bound, starting at this period.')
    Path('data/audit-extension-witness.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
