"""Independent completeness audit of the normalized equality-case catalogue.

Counts all normalized assignments without symmetry reduction, then checks
canonical representatives, stabilizers, validity and orbit-size totals.
"""
import itertools,json,sys,time,hashlib
from pathlib import Path
import numpy as np


def audit(catalogue,output):
    start=time.time();raw=Path(catalogue).read_bytes();data=json.loads(raw)
    base={7:0,10:0,14:1,15:6,21:16};mods=[30,35,42,70,105,210]
    b7=sum(1<<x for x in range(0,210,7))
    outside={d:[sum(1<<x for x in range(r,210,d) if x%7) for r in range(d)] for d in mods}
    used=sum(1<<x for x in range(210) if any(x%d==r for d,r in base.items()))
    expected_sizes=[6,6,5,3,2,1];counts=[0]*7
    def count(depth,taken):
        counts[depth]+=1
        if depth==5:
            # The last class is any one of the 105 remaining points.
            n=210-taken.bit_count();assert n==105;counts[6]+=n;return
        d=mods[depth]
        for mask in outside[d]:
            if mask.bit_count()==expected_sizes[depth] and not(mask&taken):count(depth+1,taken|mask)
    count(0,used)
    assert counts==data['normalized_assignments_by_depth']
    rows=np.array(data['cases'],dtype=np.int64);assert rows.shape==(data['n_orbits'],7)
    crt={(x%2,x%3,x%5,x%7):x for x in range(210)};permutations=[]
    for a in itertools.permutations(range(2,5)):
        for b in itertools.permutations(range(3,7)):
            pa=[0,1,*a];pb=[0,1,2,*b]
            permutations.append([crt[x%2,x%3,pa[x%5],pb[x%7]] for x in range(210)])
    maps={d:np.array([[g[r]%d for r in range(d)] for g in permutations],dtype=np.uint64) for d in mods}
    codes=set();orbit_total=0
    for lo in range(0,len(rows),1024):
        block=rows[lo:lo+1024];transformed=np.zeros((144,len(block)),dtype=np.uint64)
        original=np.zeros(len(block),dtype=np.uint64)
        for col,d in enumerate(mods):
            transformed=transformed*256+maps[d][:,block[:,col]]
            original=original*256+block[:,col].astype(np.uint64)
        assert np.array_equal(transformed.min(axis=0),original)
        stabilizers=np.sum(transformed==original,axis=0)
        assert np.array_equal(stabilizers,block[:,6])
        for row,code,stab in zip(block,original,stabilizers):
            assert int(code) not in codes;codes.add(int(code));orbit_total+=144//int(stab)
            taken=used
            for d,r,size in zip(mods,row,expected_sizes):
                mask=outside[d][int(r)];assert mask.bit_count()==size and not(mask&taken);taken|=mask
            assert taken.bit_count()==106
    assert orbit_total==counts[-1]
    result=dict(passed=True,n_orbits=len(rows),normalized_assignments=orbit_total,
                counts_without_symmetry=counts,canonical_and_stabilizer_checks=len(rows),
                minimum_uncovered=104,scope='Complete equality cases after the proved base normalization; not all prefixes',
                catalogue_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.time()-start)
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':audit(*sys.argv[1:])
