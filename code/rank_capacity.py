"""Integer superlevel-set bounds implied by subset capacity H.

If m_e points have group hit count at least t_e, then t0*m0+t1*m1<=H(m0,m1).
Maximizing a nonnegative linear weight over these feasible integer pairs
therefore bounds the corresponding weighted number of threshold events.
"""
import itertools
import numpy as np

WEIGHTS = [(1,0),(0,1),(1,1),(2,1),(1,2)]


def catalogue(H,dimensions,weights=WEIGHTS):
    H = np.array(H,dtype=np.int64)
    a,b = np.indices(H.shape)
    templates = []
    for kind,limit in ((0,14),(1,16)):
        table = H+(a+b+np.maximum(a,b) if kind else 0)
        for t0,t1 in itertools.product(range(1,limit+1),repeat=2):
            allowed = t0*a+t1*b<=table
            for w0,w1 in weights:
                bound = int((w0*a+w1*b)[allowed].max())
                templates.append((kind,t0,t1,w0,w1,bound))
    return [(h,t0,t1,w0,w1,bound) for h in range((1<<dimensions)-1)
            for kind,t0,t1,w0,w1,bound in templates if kind==int(h>=dimensions)]


def violations(witness,dimensions,rows):
    hist = np.zeros(((1<<dimensions)-1,2,17),dtype=np.int64)
    for row in witness['support']:
        for h,v in enumerate(row['counts']):
            hist[h,row['parity'],v] += row['weight']
    tails = np.cumsum(hist[:,:,::-1],axis=2)[:,:,::-1]
    scale = witness['denominator']
    result = []
    for h,t0,t1,w0,w1,bound in rows:
        used = w0*int(tails[h,0,t0])+w1*int(tails[h,1,t1])
        if used>scale*bound:
            result.append((used-scale*bound,(h,t0,t1,w0,w1,bound)))
    return sorted(result,reverse=True)
