"""Reuse actual assignments across groups with identical core resources.

Maxima over this catalogue are finite-library scores, never unrestricted
upper bounds. Every returned candidate is an explicit valid assignment.
"""
import argparse
import json
from pathlib import Path
import numpy as np


class Catalogue:
    def __init__(self,points):
        self.points = list(points)
        self.rows,self.seen,self.arrays = {},{},{}

    def add(self,shape,residues):
        shape,residues = tuple(shape),tuple(residues)
        assert len(shape)==len(residues) and all(type(a) is int and 0<=a<d for a,d in zip(residues,shape))
        counts = tuple(sum(x%d==a for a,d in zip(residues,shape)) for x in self.points)
        if counts in self.seen.setdefault(shape,set()):return 0
        self.seen[shape].add(counts)
        self.rows.setdefault(shape,[]).append((counts,residues))
        self.arrays.pop(shape,None)
        return 1

    def add_result(self,source,result):
        assert source['core_residual']==self.points
        n = len(source['high_primes'])
        added = 0
        columns = result.get('group_columns')
        if columns is None:
            groups = sorted(result['coherent_primal_witness']['groups'],key=lambda row:row['axis'])
            assert [g['axis'] for g in groups]==list(range(2**n-1))
            columns = [[a['residues'] for a in g['assignments']] for g in groups]
        for h,rows in enumerate(columns):
            shape = [d for d in range(3 if h<n else 1,211) if 210%d==0]
            for residues in rows:added += self.add(shape,residues)
        return added

    def price(self,shape,events,weights,top=4):
        shape = tuple(shape)
        rows = self.rows.get(shape,[])
        if not rows:return []
        assert len(events)==len(weights) and all(type(w) is int and w>=0 for w in weights)
        rewards = np.zeros((len(self.points),len(shape)+1),dtype=np.int64)
        assert sum(weights)<2**63
        for (x,t),w in zip(events,weights):rewards[x,t:] += w
        if shape not in self.arrays:self.arrays[shape] = np.array([c for c,r in rows],dtype=np.int16)
        counts = self.arrays[shape]
        scores = sum(rewards[x,counts[:,x]] for x in range(len(self.points)))
        indices = sorted(range(len(rows)),key=lambda j:(-int(scores[j]),j))[:top]
        return [dict(residues=list(rows[j][1]),counts=list(rows[j][0]),objective=int(scores[j])) for j in indices]

    def export(self):
        return dict(families=[dict(shapes=list(shape),assignments=[list(r) for c,r in rows])
                              for shape,rows in sorted(self.rows.items())])

    @classmethod
    def load(cls,data):
        result = cls(data['input']['core_residual'])
        for family in data['result']['families']:
            for residues in family['assignments']:result.add(family['shapes'],residues)
        return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True);ap.add_argument('--inputs',nargs='+',required=True)
    ap.add_argument('--dual');ap.add_argument('--output',required=True)
    ap.add_argument('--pool-output',help='Export all reusable assignments to every identical resource group for independent audit')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    catalogue = Catalogue(data['input']['core_residual'])
    additions = []
    for path in args.inputs:
        old = json.loads(Path(path).read_text())
        additions.append(dict(source=path,new_profiles=catalogue.add_result(old['input'],old['result'])))
    result = catalogue.export();result['sources'] = additions
    if args.dual:
        saved = json.loads(Path(args.dual).read_text());assert saved['input']==data['input']
        c = saved['result']['pricing_history'][-1]
        scored = []
        for h,row in enumerate(c['group_pricing']):
            shape = [d for d in range(3 if h<len(data['input']['high_primes']) else 1,211) if 210%d==0]
            events,weights = [(e['point'],e['threshold']) for e in row['events']],[e['weight'] for e in row['events']]
            scored.append(catalogue.price(shape,events,weights))
        numerator = sum(c['point_multipliers'])-c['point_weight_denominator']*sum(rows[0]['objective'] for rows in scored)
        result['dual_diagnostic'] = dict(source=args.dual,numerator=numerator,denominator=c['denominator'],
                                         candidates=scored,scope='lower bound evaluated on this larger finite catalogue only')
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    if args.pool_output:
        n = len(data['input']['high_primes'])
        columns = []
        for h in range(2**n-1):
            shape = tuple(d for d in range(3 if h<n else 1,211) if 210%d==0)
            columns.append([list(r) for c,r in catalogue.rows[shape]])
        Path(args.pool_output).write_text(json.dumps(dict(input=data['input'],
            result=dict(group_columns=columns),provenance=dict(catalogue=args.output)),separators=(',',':'))+'\n')
    print(json.dumps(dict(families=[dict(resources=len(f['shapes']),profiles=len(f['assignments'])) for f in result['families']],
                          additions=additions,dual_numerator=result.get('dual_diagnostic',{}).get('numerator'))))
