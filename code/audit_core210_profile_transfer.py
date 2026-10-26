"""Transfer a complete density-30 branch proof to identical residual profiles.

The two conditional reductions and free-class projection bounds depend
only on the number of residual 7-lifts at each x mod 30, provided a complete
7-plane exists. CRT-coordinate permutations preserve every cylinder type.
Other profile types and non-equality prefixes remain outside this audit.
"""
import hashlib,itertools,json,time
from pathlib import Path


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def audit():
    started=time.time();proofpath='data/audit-core2310-complete-projection-tree.json';proof=read(proofpath)
    assert proof['status']=='VERIFIED' and proof['stage11_projection_conditions']=={} and proof['base_projection_cases']==1800
    family={int(d):r for d,r in proof['fixed_core210_family'].items()};U={x for x in range(210) if all(x%d!=r for d,r in family.items())}
    assert len(U)==104;R={x%30 for x in U};complete=[c for c in range(7) if sum(x%7==c for x in U)==len(R)];assert complete
    source_complete=complete[0];profile=[sum(x%30==z for x in U) for z in range(30)]
    crt30={(x%2,x%3,x%5):x for x in range(30)};crt210={(x%30,x%7):x for x in range(210)};profiles={}
    for p2 in itertools.permutations(range(2)):
        for p3 in itertools.permutations(range(3)):
            for p5 in itertools.permutations(range(5)):
                perm=[crt30[p2[z%2],p3[z%3],p5[z%5]] for z in range(30)];v=[0]*30
                for z in range(30):v[perm[z]]=profile[z]
                profiles.setdefault(tuple(v),perm)
    catpath='data/core210-maximum-prefixes.json';cat=read(catpath);prior=read('data/audit-core210-maximum-prefixes.json')
    assert prior['passed'] and prior['catalogue_sha256']==sha(catpath)
    fixed={int(d):r for d,r in cat['normalized_base'].items()};mods=cat['remaining_moduli'];full=(1<<210)-1
    masks={d:[sum(1<<x for x in range(r,210,d)) for r in range(d)] for d in [*fixed,*mods]}
    taken=0
    for d,r in fixed.items():taken|=masks[d][r]
    p30=[sum(1<<x for x in range(z,210,30)) for z in range(30)];p7=[sum(1<<x for x in range(c,210,7)) for c in range(7)]
    matches=[]
    for index,row in enumerate(cat['cases']):
        cover=taken
        for d,r in zip(mods,row):cover|=masks[d][r]
        remaining=full^cover;counts=tuple((remaining&m).bit_count() for m in p30)
        if counts not in profiles:continue
        n=sum(v>0 for v in counts);planes=[c for c in range(7) if (remaining&p7[c]).bit_count()==n]
        if not planes:continue
        assert remaining.bit_count()==104;perm=profiles[counts];plane=planes[0]
        # All possible reference fibre-size levels in the adaptive proof.
        # This independently checks the equal multisets used by reduction.
        for z in R:
            zz=perm[z]
            for t in range(4,11):
                source=sorted(t-(c==source_complete) if crt210[z,c] in U else 0 for c in range(7))
                target=sorted(t-(c==plane) if remaining>>crt210[zz,c]&1 else 0 for c in range(7))
                assert source==target
        matches.append(dict(catalogue_index=index,complete_7_plane=plane,crt30_permutation=perm))
    result=dict(status='VERIFIED',claim='These equality-case orbits inherit the complete later-stage exclusion from the proved residual profile',
                source_proof=proofpath,source_proof_sha256=sha(proofpath),catalogue_sha256=sha(catpath),
                excluded_orbits=len(matches),total_equality_orbits=cat['n_orbits'],matches=matches,seconds=time.time()-started,
                scope='Only the listed equality-case orbits; other profiles and larger residuals remain, so global upper bound is still 7')
    Path('data/audit-core210-profile-transfer.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='matches'},indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit()
