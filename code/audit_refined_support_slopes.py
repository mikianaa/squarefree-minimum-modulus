"""Audit exact rational exceptional-plane slopes over the entire plane domain."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import numpy as np
from audit_positive_state_dual import audit as audit_parent,minimum


def integer_minimum(target,volumes,penalties):
    """Forward exact-credit DP, independently of the generator's requirement DP."""
    assert type(target) is int and target>=0
    assert len(volumes)==len(penalties) and all(type(v) is int and v>=0 for v in volumes)
    assert all(row and row[0]==0 and all(type(v) is int and v>=0 for v in row) for row in penalties)
    infinity=1+sum(max(row) for row in penalties);assert 2*infinity<2**63
    dp=np.full(target+1,infinity,dtype=np.int64);dp[0]=0
    for volume,row in zip(volumes,penalties):
        nxt=np.full(target+1,infinity,dtype=np.int64)
        for k,penalty in enumerate(row):
            credit=volume*k
            if credit>=target:nxt[target]=min(int(nxt[target]),int(dp.min())+penalty)
            else:
                nxt[credit:target]=np.minimum(nxt[credit:target],dp[:target-credit]+penalty)
                nxt[target]=min(int(nxt[target]),int(dp[target-credit:].min())+penalty)
        dp=nxt
    return int(dp[target]) if dp[target]<infinity else None


def audit(data):
    parent=audit_parent(data)
    source,result=data['input'],data['result']
    base=result['pricing_history'][-1];old=result['positive_state_pricing']
    assert old['mode']=='nonpositive_support'
    integer_mode='integer_support_refinement' in result
    refined=result['integer_support_refinement' if integer_mode else 'refined_support_slopes'];assert refined['complete']
    rows=refined['points'];R=source['core_residual']
    assert [r['point_index'] for r in rows]==list(range(len(R)))
    minima=[];checks=0;exceptions=0;integer_checks=0
    for x,(point,row) in enumerate(zip(R,rows)):
        D=[q-1-(b==point%2) for q,b in zip(source['high_primes'],source['high_pattern'])];n=len(D)
        potentials=[[base['point_weight_denominator']*sum(g[x][:t+1]) for t in range(len(g[x]))] for g in base['event_multipliers']]
        envelope,free,groups=minimum(D,potentials,old['point_slope_numerators'][x],old['denominator'],True,True)
        max_credit=np.zeros(len(envelope),dtype=np.int64)
        for S in groups[n:]:
            volume=np.ones(len(envelope),dtype=np.int64)
            for i in range(n):
                if i not in S:volume*=free[i]
            max_credit+=16*volume
        feasible=np.prod(np.stack(free),axis=0)<=max_credit
        assert int(feasible.sum())==row['feasible_plane_count'] and int((~feasible).sum())==row['excluded_plane_count']
        bound=row['lower_numerator'];assert type(bound) is int and bound>=old['point_multipliers'][x]
        checked=np.zeros(len(envelope),dtype=bool)
        for special in row['refined_planes']:
            plane=special['plane'];assert len(plane)==n and all(type(v) is int and 0<=v<=min(14,d) for v,d in zip(plane,D))
            slope=Fraction(*special['slope']);assert slope>=0
            f=[d-v for d,v in zip(D,plane)]
            plane_penalty=-sum(potentials[h][v] for h,v in enumerate(plane))
            value=slope*math.prod(f)+plane_penalty
            volumes=[]
            for h,S in enumerate(groups[n:],n):
                volume=math.prod(f[i] for i in range(n) if i not in S)
                volumes.append(volume)
                value+=min(-potentials[h][k]-slope*volume*k for k in range(17))
            if integer_mode and 'integer_cost' in special:
                penalties=[[-v for v in line] for line in potentials[n:]]
                exact=integer_minimum(math.prod(f),volumes,penalties)
                assert exact is not None and plane_penalty+exact==special['integer_cost']>=bound
                tail=special['tail_counts'];assert len(tail)==len(volumes) and all(type(v) is int and 0<=v<=16 for v in tail)
                assert sum(v*k for v,k in zip(volumes,tail))>=math.prod(f)
                assert sum(line[k] for line,k in zip(penalties,tail))==exact
                integer_checks+=1
            elif integer_mode:assert -(-value.numerator//value.denominator)>=bound
            else:assert value>=bound
            index=np.ravel_multi_index(tuple(v-max(0,d-14) for v,d in zip(f,D)),tuple(min(14,d)+1 for d in D))
            checked[index]=True;exceptions+=1
        assert np.all(checked|(~feasible)|(envelope>=bound))
        if integer_mode:
            state=row['attaining_state'];assert len(state)==len(groups)
            assert all(type(v) is int and 0<=v<=(min(14,D[h]) if h<n else 16) for h,v in enumerate(state))
            f=[d-v for d,v in zip(D,state)]
            assert math.prod(f)<=sum(state[h]*math.prod(f[i] for i in range(n) if i not in S) for h,S in enumerate(groups[n:],n))
            assert -sum(potentials[h][v] for h,v in enumerate(state))==bound
        minima.append(bound);checks+=len(envelope)
    gain=sum(minima)-sum(old['point_multipliers'])
    pool=parent['pool_numerator']+gain if parent['pool_numerator'] is not None else None
    capacity=parent['unrestricted_capacity_numerator']+gain
    return dict(verified=True,point_multipliers=minima,denominator=old['denominator'],pool_numerator=pool,
                unrestricted_capacity_numerator=capacity,gain_over_parent_numerator=gain,
                plane_states_checked=checks,rational_slopes_checked=exceptions,
                integer_dp_checks=integer_checks,exact_integer_support_minima=integer_mode,
                excludes_stated_pool=pool is not None and pool>0,excludes_unrestricted_core_pattern=capacity>0,
                scope=('Exact integer event-penalty minimum at each point; independent attaining states do not form a covering'
                       if integer_mode else 'Certified nonpositive-support lower bound; fractional tail states are not integer covering assignments'))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();report=audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
