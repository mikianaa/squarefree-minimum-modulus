"""Incompatible high-hit events at two specified core points.

For k projection shapes, z(x)+z(y)<=k+c where c shapes identify x and y.
Thus [z(x)>=t]+[z(y)>=u]<=1 whenever t+u>k+c. It suffices to retain
the boundary t+u=k+c+1; larger threshold pairs are weaker.
"""
import itertools
import numpy as np


def catalogue(R,dimensions):
    core = [d for d in range(3,211) if 210%d==0]
    result = []
    for h in range((1<<dimensions)-1):
        shapes = core if h<dimensions else [1,2]+core
        limit = len(shapes)
        for x,y in itertools.combinations(range(len(R)),2):
            cap = limit+sum(R[x]%d==R[y]%d for d in shapes)
            for t in range(1,limit+1):
                u = cap+1-t
                if 1<=u<=limit:
                    result.append((h,x,y,t,u))
    return result


def violations(witness,R,dimensions,rows):
    hist = np.zeros((len(R),(1<<dimensions)-1,17),dtype=np.int64)
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[row['point_index'],h,v] += row['weight']
    tails = np.cumsum(hist[:,:,::-1],axis=2)[:,:,::-1]
    result = []
    for h,x,y,t,u in rows:
        gap = int(tails[x,h,t])+int(tails[y,h,u])-witness['denominator']
        if gap>0:
            result.append((gap,(h,x,y,t,u)))
    return sorted(result,reverse=True)
