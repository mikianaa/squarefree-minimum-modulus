"""Certified case splitting of affine bounds on actual residue assignments."""
import argparse
import heapq
import json
from pathlib import Path
from event_score_majorant import optimized_upper_bound


def shifted(points,shapes,events,fixed):
    used=dict(fixed);remaining=[d for d in shapes if d not in used]
    hits=[sum(v%d==a for d,a in fixed) for v in points]
    constant=0;reduced=[]
    for e in events:
        threshold=e['threshold']-hits[e['point']]
        if threshold<=0:constant+=e['weight']
        elif threshold<=len(remaining):reduced.append(dict(point=e['point'],threshold=threshold,weight=e['weight']))
    return remaining,reduced,constant


def bound(points,shapes,events,depth=3,max_nodes=2000,pairs=False,cliques=False):
    assert points and len(set(shapes))==len(shapes)
    assert all(type(e['weight']) is int and e['weight']>=0 for e in events)
    nodes=[];heap=[]
    def add(fixed,parent_bound=None):
        remaining,reduced,constant=shifted(points,shapes,events,fixed)
        majorant=optimized_upper_bound
        if pairs:
            from pair_event_majorant import upper_bound as majorant
        if cliques:
            from clique_event_majorant import upper_bound as majorant
        certificate=majorant(points,remaining,[(e['point'],e['threshold']) for e in reduced],[e['weight'] for e in reduced])
        value=constant+certificate['bound']
        effective=min(value,parent_bound) if parent_bound is not None else value
        j=len(nodes);nodes.append(dict(fixed=fixed,constant=constant,affine=certificate,local_bound=value,effective_bound=effective))
        heapq.heappush(heap,(-effective,j));return j
    add([])
    while heap:
        negative,j=heapq.heappop(heap);node=nodes[j]
        used=dict(node['fixed']);remaining=[d for d in shapes if d not in used and d!=1]
        if len(used)>=depth or not remaining:break
        d=remaining[0];residues=sorted({x%d for x in points})
        if len(nodes)+len(residues)>max_nodes:break
        node['resource']=d
        node['children']=[add(node['fixed']+[[d,a]],node['effective_bound']) for a in residues]
    leaves=[row['effective_bound'] for row in nodes if 'children' not in row]
    return dict(method='branch_affine',bound=max(leaves),nodes=nodes,maximum_depth=depth,
                scope='All actual assignments; unoccupied residues are dominated because rewards are nonnegative')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--depth',type=int,default=3);ap.add_argument('--axis',type=int)
    ap.add_argument('--max-nodes',type=int,default=2000)
    ap.add_argument('--pairs',action='store_true')
    ap.add_argument('--cliques',action='store_true')
    args=ap.parse_args();data=json.loads(Path(args.source).read_text())
    source=data['input'];c=data['result']['pricing_history'][-1];n=len(source['high_primes']);rows=[]
    for h,row in enumerate(c['group_pricing']):
        if args.axis is not None and h!=args.axis:continue
        shape=[d for d in range(3 if h<n else 1,211) if 210%d==0]
        certificate=bound(source['core_residual'],shape,row['events'],args.depth,args.max_nodes,args.pairs,args.cliques)
        rows.append(dict(axis=h,certificate=certificate))
        output=dict(input=source,source=args.source,group_bounds=rows,complete=len(rows)==len(c['group_pricing']))
        Path(args.output).write_text(json.dumps(output,separators=(',',':'))+'\n')
        print(json.dumps(dict(axis=h,old_bound=row['capacity_certificate']['bound'],bound=certificate['bound'],nodes=len(certificate['nodes']))),flush=True)
