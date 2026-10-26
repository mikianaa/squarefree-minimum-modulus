"""Integer witness: the free-class union bound is sharp for one fixed Q.

This is not a covering or a method barrier for other measures/schedules.
No solver status or numerical objective is trusted.
"""
import json,hashlib
from pathlib import Path


def audit():
    path=Path('data/core2310-hard-projection-fibre-envelope.json');core=json.loads(path.read_text());data=json.loads(Path(core['measure_file']).read_text())
    probe=json.loads(Path('data/core2310-capped-free-union-probe.json').read_text());classes={int(d):r for d,r in probe['actual_residue_witness'].items()}
    assert set(classes)==set(data['free_moduli']);w=core['weights_numerators'];den=core['weights_denominator'];covered=set();caps={}
    for d,r in classes.items():
        assert 0<=r<d
        points={x for x in range(2310) if x%d==r and w[x]>0};assert not covered&points
        cap=max(sum(w[x] for x in range(a,2310,d)) for a in range(d));assert sum(w[x] for x in points)==cap
        caps[d]=cap;covered|=points
    total=sum(w[x] for x in covered);assert total==sum(caps.values())
    result=dict(status='VERIFIED',claim='The union bound for these nine free moduli is attained by actual disjoint classes on the specified measure',
                class_witness=classes,union_numerator=total,denominator=den,union_mass=total/den,
                core_file=str(path),core_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='One fixed measure only. Does not exclude better weights, future overlap accounting, or a smaller global bound.')
    Path('data/audit-core2310-free-union-sharp.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit()
