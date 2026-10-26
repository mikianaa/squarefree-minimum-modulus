"""Independently reconstruct selected case rewards and verify rational bounds."""
import argparse
import json
from pathlib import Path
from audit_weighted_coherent_dual import capacity


def audit(data):
    source=data['input'];R=source['core_residual'];n=len(source['high_primes']);G=2**n-1
    assert len(data['group_events'])==G
    improved=0;cuts=0;cliques=0;maximum=0
    for row in data['cases']:
        h=row['axis'];assert type(h) is int and 0<=h<G
        shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        fixed=row['fixed'];assert len(dict(fixed))==len(fixed)
        assert all(type(d) is int and d in shape and type(a) is int and 0<=a<d for d,a in fixed)
        remaining=[d for d in shape if d not in dict(fixed)]
        hits=[sum(p%d==a for d,a in fixed) for p in R];constant=0;events=[]
        for e in data['group_events'][h]:
            x,t,w=e['point'],e['threshold'],e['weight']
            assert type(x) is int and 0<=x<len(R) and type(t) is int and 1<=t<=len(shape) and type(w) is int and w>=0
            if t<=hits[x]:constant+=w
            elif t<=hits[x]+len(remaining):events.append(dict(point=x,threshold=t-hits[x],weight=w))
        assert constant==row['constant']
        old=constant+capacity(R,remaining,events,row['old_certificate'])
        new=constant+capacity(R,remaining,events,row['new_certificate'])
        assert old==row['old_bound'] and new==row['new_bound'] and new<=old
        improved+=new<old;maximum=max(maximum,old-new);cuts+=len(row['new_certificate'].get('pair_cuts',[]))
        cliques+=len(row['new_certificate'].get('clique_cuts',[]))
    return dict(verified=True,cases_checked=len(data['cases']),strictly_improved_cases=improved,
                positive_pair_multipliers_checked=cuts,maximum_bound_decrease=maximum,
                positive_clique_multipliers_checked=cliques,
                complete_group_coverage=False,covering_obstruction=False,
                scope='Every listed conditional bound is verified; this list is not a full group case partition')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();report=audit(json.loads(Path(args.source).read_text()))
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
