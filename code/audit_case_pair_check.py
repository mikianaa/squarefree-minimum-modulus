"""Exact verification of pair incompatibilities violated by relaxed case states."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from audit_case_relaxation_barrier import audit as audit_barrier


def score(points,shapes,scenario,x,y,t,u,conditional=False):
    assert 0<=x<y<len(points) and 1<=t<=len(shapes) and 1<=u<=len(shapes)
    fixed=scenario['fixed'];remaining=[d for d in shapes if d not in dict(fixed)]
    hits=[sum(p%d==a for d,a in fixed) for p in points];base=int(1 in remaining)
    capacity=(hits[x]+hits[y]+len(remaining)+sum(points[x]%d==points[y]%d for d in remaining)
              if conditional else len(shapes)+sum(points[x]%d==points[y]%d for d in shapes))
    assert t+u>capacity
    total=Fraction(0)
    for i,threshold in ((x,t),(y,u)):
        probabilities=[Fraction(*p) for p in scenario['count_distributions'][i]]
        assert len(probabilities)==len(remaining)-base+1 and min(probabilities)>=0 and sum(probabilities)==1
        total+=sum(p for j,p in enumerate(probabilities) if hits[i]+base+j>=threshold)
    return total,capacity


def audit(barrier,check):
    parent=audit_barrier(barrier);assert parent['verified']
    R=barrier['input']['core_residual'];n=len(barrier['input']['high_primes'])
    groups={row['axis']:row for row in barrier['group_mixtures']};counts=[];largest=Fraction(0)
    for key,conditional in (('individual_worst_global_pair_violations',False),('individual_worst_conditional_pair_violations',True)):
        for row in check[key]:
            h=row['axis'];shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
            scenario=groups[h]['scenarios'][row['scenario']]
            value,capacity=score(R,shapes,scenario,*row['points'],*row['thresholds'],conditional)
            gap=value-1;assert gap==Fraction(*row['excess']) and gap>0
        counts.append(len(check[key]))
    for row in check['mixture_global_pair_violations']:
        h=row['axis'];group=groups[h];shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        value=sum(Fraction(s['weight'],group['denominator'])*score(R,shapes,s,*row['points'],*row['thresholds'])[0]
                  for s in group['scenarios'])
        gap=value-1;assert gap==Fraction(*row['excess']) and gap>0;largest=max(largest,gap)
    return dict(verified=True,individual_global_violations=counts[0],individual_conditional_violations=counts[1],
                mixture_violations=len(check['mixture_global_pair_violations']),largest_mixture_excess=[largest.numerator,largest.denominator],
                excludes_this_fractional_barrier_witness=largest>0,covering_obstruction=False,
                scope='Checks claimed violations, not completeness of the violation list; the old witness fails stronger pair compatibility')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','check','output'):ap.add_argument('--'+key,required=True)
    args=ap.parse_args();report=audit(json.loads(Path(args.source).read_text()),json.loads(Path(args.check).read_text()))
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
