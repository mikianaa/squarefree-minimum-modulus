"""Enumerate equality cases of the elementary sharp 104-point residual bound.

After normalizing 7,10,14,15,21, the six remaining classes have disjoint
parts outside modulus 7. Residual coordinate permutations give a finite
canonical enumeration. No optimization solver is used.
"""
import itertools,json,sys,time
from pathlib import Path

BASE={7:0,10:0,14:1,15:6,21:16}
REST=[30,35,42,70,105,210]
EXPECTED=[6,6,5,3,2,1]


def enumerate_cases(output):
    start=time.time();points=list(range(210));bit=lambda xs:sum(1<<x for x in xs)
    b7=bit(x for x in points if x%7==0)
    initial=bit(x for x in points if any(x%d==r for d,r in BASE.items()))
    assert initial.bit_count()==83
    group=[]
    for perm5 in itertools.permutations([2,3,4]):
        for perm7 in itertools.permutations([3,4,5,6]):
            p5=[0,1,*perm5];p7=[0,1,2,*perm7]
            bycrt={(x%2,x%3,x%5,x%7):x for x in points}
            group.append([bycrt[x%2,x%3,p5[x%5],p7[x%7]] for x in points])
    images={d:[[g[r]%d for r in range(d)] for g in group] for d in REST}
    masks={d:[bit(x for x in points if x%d==r)&~b7 for r in range(d)] for d in REST}
    cases=[];depth_counts=[0]*7;weighted=[0]*7
    def visit(depth,family,covered,stabilizer):
        depth_counts[depth]+=1;weighted[depth]+=144//len(stabilizer)
        if depth==6:
            support=[x for x in points if not(covered>>x&1)]
            assert len(support)==104
            cases.append([family[d] for d in REST]+[len(stabilizer)])
            return
        d=REST[depth]
        for r,mask in enumerate(masks[d]):
            if mask.bit_count()!=EXPECTED[depth] or mask&covered:continue
            if r!=min(images[d][g][r] for g in stabilizer):continue
            sub=[g for g in stabilizer if images[d][g][r]==r]
            visit(depth+1,{**family,d:r},covered|mask,sub)
    visit(0,{},initial,list(range(144)))
    result=dict(minimum_uncovered=104,normalized_base=BASE,remaining_moduli=REST,
                stabilizer_group_size=144,orbit_counts_by_depth=depth_counts,
                normalized_assignments_by_depth=weighted,n_orbits=len(cases),case_columns=[*REST,'stabilizer_size'],cases=cases,
                scope='All equality cases for the 104-point core bound, not all possible prefixes',seconds=time.time()-start)
    Path(output).write_text(json.dumps(result,separators=(',',':')))
    print({k:v for k,v in result.items() if k!='cases'},flush=True)
    return result


def audit_bound(witness_file,output):
    witness=json.loads(Path(witness_file).read_text());family={int(d):r for d,r in witness['family'].items()}
    ds=[d for d in range(7,211) if 210%d==0]
    assert sorted(family)==ds
    covered={x for x in range(210) if any(x%d==r for d,r in family.items())}
    assert len(covered)==106
    assert [x for x in range(210) if x not in covered]==witness['support']
    # Check both unavoidable intersections outside a normalized 7-plane.
    pair_checks=0
    for d,e in [(10,21),(15,14)]:
        for a in range(d):
            for b in range(e):
                if b%7==0:continue
                intersection=[x for x in range(210) if x%d==a and x%e==b and x%7!=0]
                assert len(intersection)==1;pair_checks+=1
    sizes={d:max(sum(x%d==r and x%7!=0 for x in range(210)) for r in range(d)) for d in ds if d!=7}
    assert sum(sizes.values())==78
    assert sizes[14]==15 and sizes[21]==10
    # If either 14 or 21 lies inside the 7-plane, losing its entire
    # outside capacity is already stronger than subtracting two points.
    for d in (14,21):assert 180-(78-sizes[d])>=104
    bound=180-(78-2);assert bound==104
    result=dict(passed=True,minimum_uncovered=104,witness_uncovered=210-len(covered),
                pair_checks=pair_checks,outside_capacities=sizes,
                lower_bound_derivation='180-(78-1-1)=104; the two pairs use disjoint modulus sets',
                scope='Sharp finite-core residual theorem, not a bound on the full infinite-prime covering problem')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    if sys.argv[1]=='enumerate':enumerate_cases(sys.argv[2])
    else:audit_bound(sys.argv[2],sys.argv[3])
