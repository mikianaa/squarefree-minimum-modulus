"""Full-period integer audit of saved six-, seven-, and eight-prime witnesses."""
import hashlib
import json
from pathlib import Path


def audit(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    n = data['N']
    remaining = n
    primes = []
    p = 2
    while p*p<=remaining:
        if remaining%p==0:
            primes.append(p)
            remaining //= p
            assert remaining%p!=0
        p += 1
    if remaining>1:
        primes.append(remaining)
    assert primes[0]==2
    system = [(r['modulus'],r['residue']) for r in data['system']]
    family = dict(system)
    assert len(family)==len(system)==2**len(primes)-2
    assert data['pattern']==[family[2*p]%2 for p in primes[1:]]
    assert all(family[p]==0 and family[2*p]%p==1 for p in primes[1:])
    covered = bytearray(n)
    for d,a in system:
        assert d>=3 and n%d==0 and 0<=a<d
        for x in range(a,n,d):
            covered[x] = 1
    left = [covered[e::2].count(0) for e in (0,1)]
    assert sum(left)==data['uncovered']
    return dict(file=path.as_posix(),period=n,classes=len(system),uncovered=sum(left),
                parity_uncovered=left,verified=True,sha256=hashlib.sha256(raw).hexdigest(),
                scope='upper witness only, not an optimality or covering certificate')


if __name__=='__main__':
    results = [audit(Path(f'data/min-uncovered-{n}-witness.json')) for n in (30030,510510,9699690)]
    Path('data/audit-partial-families.json').write_text(json.dumps(results,indent=2)+'\n')
    for row in results:
        print(json.dumps(row))
