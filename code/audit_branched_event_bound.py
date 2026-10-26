"""Solver-free coverage and rational-bound audit of residue case splits."""
from audit_weighted_coherent_dual import capacity


def audit_bound(points,shapes,events,certificate):
    assert certificate['method']=='branch_affine' and points
    assert all(type(e['point']) is int and 0<=e['point']<len(points) and type(e['threshold']) is int
               and 1<=e['threshold']<=len(shapes) and type(e['weight']) is int and e['weight']>=0 for e in events)
    nodes=certificate['nodes'];seen=set()
    def visit(index,fixed,ancestor=None):
        assert type(index) is int and 0<=index<len(nodes) and index not in seen
        seen.add(index);node=nodes[index];assert node['fixed']==fixed
        used=dict(fixed);assert len(used)==len(fixed) and all(d in shapes and type(a) is int and 0<=a<d for d,a in fixed)
        remaining=[d for d in shapes if d not in used]
        counts=[sum(v%d==a for d,a in fixed) for v in points]
        constant=sum(e['weight'] for e in events if e['threshold']<=counts[e['point']])
        reduced=[dict(point=e['point'],threshold=e['threshold']-counts[e['point']],weight=e['weight']) for e in events
                 if counts[e['point']]<e['threshold']<=counts[e['point']]+len(remaining)]
        local=constant+capacity(points,remaining,reduced,node['affine'])
        assert constant==node['constant'] and local==node['local_bound']
        effective=min(local,ancestor) if ancestor is not None else local
        assert effective==node['effective_bound']
        if 'children' not in node:return effective
        d=node['resource'];assert d in remaining
        residues=sorted({x%d for x in points});assert len(node['children'])==len(residues)
        return max(visit(child,fixed+[[d,a]],effective) for child,a in zip(node['children'],residues))
    value=visit(0,[]);assert seen==set(range(len(nodes))) and value==certificate['bound']
    return value
