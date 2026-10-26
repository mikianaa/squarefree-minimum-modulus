"""Reconstruct a numerical curve-cone barrier with exact integer probabilities."""
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path


def quantize(values,denominator):
    values=[Fraction(max(0,float(v))) for v in values];total=sum(values);assert total>0
    scaled=[v*denominator/total for v in values]
    units=[v.numerator//v.denominator for v in scaled]
    missing=denominator-sum(units)
    for j in sorted(range(len(units)),key=lambda j:scaled[j]-units[j],reverse=True)[:missing]:units[j]+=1
    assert sum(units)==denominator and min(units)>=0
    return units


def tails(points,shapes,fixed,counts):
    remaining=[d for d in shapes if d not in dict(fixed)];base=int(1 in remaining)
    hits=[sum(p%d==a for d,a in fixed) for p in points]
    return [[sum(v for j,v in enumerate(row) if hits[x]+base+j>=t) for t in range(len(shapes)+1)]
            for x,row in enumerate(counts)]


def pair_maximum(points,shapes,fixed,counts):
    remaining=[d for d in shapes if d not in dict(fixed)];N=len(shapes)
    hit=[sum(p%d==a for d,a in fixed) for p in points]
    tail=tails(points,shapes,fixed,counts);largest=0
    for x,y in itertools.combinations(range(len(points)),2):
        cap=hit[x]+hit[y]+len(remaining)+sum(points[x]%d==points[y]%d for d in remaining)
        for t in range(max(1,cap+1-N),min(N,cap)+1):
            largest=max(largest,tail[x][t]+tail[y][cap+1-t])
    return largest


def rational_case(points,shapes,scenario,D):
    fixed=scenario['fixed'];resource_rows=[];supply=[0]*len(points)
    for row in scenario['resource_distributions']:
        units=quantize(row['probabilities'],D);table=dict(zip(row['residues'],units));d=row['resource']
        resource_rows.append(dict(resource=d,residues=row['residues'],probability_numerators=units))
        for x,p in enumerate(points):supply[x]+=table.get(p%d,0)
    counts=[];mean_repaired=0
    for x,row in enumerate(scenario['count_distributions']):
        units=quantize(row,D);excess=sum(j*v for j,v in enumerate(units))-supply[x]
        if excess>0:
            mean_repaired+=1
            for j in range(len(units)-1,0,-1):
                moved=min(units[j],max(0,(excess+j-1)//j))
                units[j]-=moved;units[0]+=moved;excess-=j*moved
            assert excess<=0
        counts.append(units)
    maximum=pair_maximum(points,shapes,fixed,counts)
    if maximum>D:
        counts=[[D-sum(v*D//maximum for v in row[1:])]+[v*D//maximum for v in row[1:]] for row in counts]
    assert pair_maximum(points,shapes,fixed,counts)<=D
    return dict(kind='fractional',fixed=fixed,resource_distributions=resource_rows,count_numerators=counts,
                repair=dict(mean_rows=mean_repaired,old_pair_maximum=maximum))


def certify(data,scenarios,basis,probe,D=10**8,mix_denominator=10**10):
    source=data['input'];assert source==scenarios['input']==probe['input']
    assert all(entry['input']==source for entry in basis)
    assert probe['numerical_margin']>0 and type(D) is int and D>0
    R=source['core_residual'];n=len(source['high_primes']);G=2**n-1
    events=[entry['result']['pricing_history'][-1]['event_multipliers'] for entry in [data]+basis]
    groups={row['axis']:row['scenarios'] for row in scenarios['group_scenarios']}
    point_mix=[];group_mix=[]
    for row in probe['point_mixtures']:
        weights=quantize([s['weight'] for s in row['states']],mix_denominator)
        point_mix.append(dict(point_index=row['point_index'],denominator=mix_denominator,
            states=[dict(weight=w,counts=s['counts']) for w,s in zip(weights,row['states']) if w]))
    for group in probe['group_mixtures']:
        h=group['axis'];shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        weights=quantize([s['weight'] for s in group['states']],mix_denominator);rows=[]
        for w,s in zip(weights,group['states']):
            if not w:continue
            if s['kind']=='actual':row=dict(kind='actual',assignment=data['result']['group_columns'][h][s['index']])
            else:row=rational_case(R,shapes,groups[h][s['index']],D)
            row['weight']=w;row['source_index']=s['index'];rows.append(row)
        group_mix.append(dict(axis=h,denominator=mix_denominator,states=rows))
    # Accumulate integer tail numerators. All count distributions have the same
    # probability denominator; no large least-common-multiple reconstruction.
    demand=[[[0]*(15 if h<n else 17) for p in R] for h in range(G)]
    point_den={row['point_index']:row['denominator'] for row in point_mix}
    for row in point_mix:
        x=row['point_index']
        for s in row['states']:
            for h,k in enumerate(s['counts']):
                for t in range(k+1):demand[h][x][t]+=s['weight']
    supply=[[[0]*(15 if h<n else 17) for p in R] for h in range(G)]
    group_den={row['axis']:row['denominator']*D for row in group_mix}
    for group in group_mix:
        h=group['axis'];shapes=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        for s in group['states']:
            if s['kind']=='actual':
                counts=[sum(p%d==a for d,a in zip(shapes,s['assignment'])) for p in R]
                tail=[[D*int(k>=t) for t in range(len(shapes)+1)] for k in counts]
            else:tail=tails(R,shapes,s['fixed'],s['count_numerators'])
            for x,line in enumerate(tail):
                for t,v in enumerate(line):supply[h][x][t]+=s['weight']*v
    surpluses=[];normalized=[]
    for b,axes in enumerate(events):
        for h,axis in enumerate(axes):
            for x,row in enumerate(axis):
                mass=-sum(row)
                if not mass or not -sum(events[0][h][x]):continue
                gap=Fraction(sum(-v*supply[h][x][t] for t,v in enumerate(row)),group_den[h])-Fraction(
                    sum(-v*demand[h][x][t] for t,v in enumerate(row)),point_den[x])
                assert gap>=0,'Rational reconstruction lost coordinate dominance'
                surpluses.append(dict(basis=b,axis=h,point_index=x,surplus=[gap.numerator,gap.denominator]))
                normalized.append(gap/mass)
    margin=min(normalized);assert margin>0
    return dict(input=source,basis_event_multipliers=events,probability_denominator=D,
                point_mixtures=point_mix,group_mixtures=group_mix,coordinate_surpluses=surpluses,
                minimum_normalized_surplus=[margin.numerator,margin.denominator],
                scope='All nonnegative combinations of the specified curves at positive-baseline-cap coordinates, when group bounds admit these explicit actual states and pair-compatible fractional case states. Not an actual covering, arbitrary-curve barrier, or global method optimality theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('source','scenarios','probe','output'):ap.add_argument('--'+key,required=True)
    ap.add_argument('--basis',action='append',default=[]);ap.add_argument('--probability-denominator',type=int,default=10**8)
    args=ap.parse_args();result=certify(json.loads(Path(args.source).read_text()),json.loads(Path(args.scenarios).read_text()),
        [json.loads(Path(p).read_text()) for p in args.basis],json.loads(Path(args.probe).read_text()),args.probability_denominator)
    result['sources']=dict(source=args.source,scenarios=args.scenarios,probe=args.probe,basis=args.basis)
    Path(args.output).write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps(dict(coordinates=len(result['coordinate_surpluses']),minimum_normalized_surplus=result['minimum_normalized_surplus'])))
