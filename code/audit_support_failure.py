"""Integer check that particular weights cannot separate actual group choices.

One valid local F<=0 state per point upper-bounds the exact local minimum;
one actual assignment per group lower-bounds any unrestricted group maximum.
This checker needs no optimizer, exact local minimization, or numerical bound.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def audit(data):
    source=data['input'];R=source['core_residual'];q=source['high_primes'];pattern=source['high_pattern'];n=len(q)
    assert R and len(set(R))==len(R) and all(type(p) is int and 0<=p<210 for p in R)
    assert len(pattern)==n and all(b in (0,1) for b in pattern)
    assert q==sorted(set(q)) and all(type(p) is int and p>2 and all(p%d for d in range(2,math.isqrt(p)+1)) for p in q)
    subsets=[s for k in range(1,n+1) for s in itertools.combinations(range(n),k)]
    G=len(subsets);events=data['event_multipliers'];assignments=data['group_assignments']
    assert len(events)==len(assignments)==G
    shapes=[[d for d in range(3 if h<n else 1,211) if 210%d==0] for h in range(G)]
    for h,axis in enumerate(events):
        assert len(axis)==len(R)
        for row in axis:
            assert len(row)==len(shapes[h])+1 and row[0]==0 and all(type(v) is int and v<=0 for v in row)
    point_rows=data['point_states'];assert sorted(r['point_index'] for r in point_rows)==list(range(len(R)))
    local=0
    for row in point_rows:
        x=row['point_index'];k=row['counts'];D=[p-1-int(b==R[x]%2) for p,b in zip(q,pattern)]
        assert len(k)==G and all(type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16) for h,v in enumerate(k))
        f=[d-v for d,v in zip(D,k)]
        assert math.prod(f)<=sum(k[h]*math.prod(f[i] for i in range(n) if i not in S) for h,S in enumerate(subsets[n:],n))
        local-=sum(sum(events[h][x][:v+1]) for h,v in enumerate(k))
    group=0
    for h,residues in enumerate(assignments):
        assert len(residues)==len(shapes[h]) and all(type(a) is int and 0<=a<d for d,a in zip(shapes[h],residues))
        for x,p in enumerate(R):
            count=sum(p%d==a for d,a in zip(shapes[h],residues))
            group-=sum(events[h][x][:count+1])
    difference=local-group
    assert difference==data['claimed_score_difference']
    D=data['denominator'];assert type(D) is int and D>0
    return dict(verified=True,local_state_score=local,actual_group_score=group,score_difference=difference,
                denominator=D,these_weights_cannot_exclude_unrestricted_core=difference<=0,
                local_states_checked=len(point_rows),actual_group_assignments_checked=G,
                covering_constructed=False,global_problem_settled=False,
                scope='Only the specified fixed weights: exact local minima cannot exceed the displayed local-state score, and any unrestricted group bound must admit the displayed actual assignments.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();report=audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
