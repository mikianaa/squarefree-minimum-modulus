"""Independent integer enumeration of the no-modulus-10 core envelope.

Normalize the modulus-15 plane to x3=x5=0, and the nonredundant
modulus-30 point to x2=0, (x3,x5)=(0,1),(1,0), or (1,1).
Enumerate all 810000 choices of one cylinder per nontrivial divisor of 30.
No orbit catalogue, fibre rule, or primary count constructor is imported.
"""
import json,itertools,hashlib,sys
from pathlib import Path
import numpy as np


def audit(output):
    points=list(itertools.product(range(2),range(3),range(5)))
    subsets=[(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]
    sizes=[2,3,5,6,10,15,30]
    tables=[]
    for subset in subsets:
        residues=list(itertools.product(*[range((2,3,5)[i]) for i in subset]))
        tables.append(np.array([[all(x[i]==a for i,a in zip(subset,r)) for x in points] for r in residues],dtype=np.int64))
    masks=[]
    for a,b in [(0,1),(1,0),(1,1)]:
        masks.append(np.array([not(x[1]==0 and x[2]==0 or x==(0,a,b)) for x in points],dtype=np.int64))
    W=np.array(masks).T
    assert np.all(W.sum(axis=0)==27)
    maxima=np.zeros((3,9),dtype=np.int64);args=np.zeros((3,9),dtype=np.int64)
    total=int(np.prod(sizes))
    for start in range(0,total,4096):
        ix=np.arange(start,min(total,start+4096));digits=ix.copy()
        counts=np.ones((len(ix),30),dtype=np.int64)
        for size,tab in zip(sizes,tables):
            counts+=tab[digits%size];digits//=size
        assert np.all(digits==0)
        for j in range(8):
            values=np.maximum(counts-j,0)@W
            which=np.argmax(values,axis=0);tops=values[which,np.arange(3)]
            better=tops>maxima[:,j];maxima[better,j]=tops[better];args[better,j]=ix[which[better]]
    expected=[68,41,20,12,4,3,2,1,0]
    assert maxima.tolist()==[expected]*3
    law=[expected[j-1]-2*expected[j]+expected[j+1] for j in range(1,8)]+[expected[7]]
    assert law==[6,13,0,7,0,0,0,1] and sum(law)==27
    # Re-evaluate each maximizing catalogue using ordinary Python point tests.
    for case,mask in enumerate(masks):
        for j in range(8):
            index=int(args[case,j]); choices=[]
            for size in sizes:choices.append(index%size);index//=size
            value=0
            for x,keep in zip(points,mask):
                if not keep:continue
                count=1
                for subset,choice in zip(subsets,choices):
                    residues=list(itertools.product(*[range((2,3,5)[i]) for i in subset]))
                    count+=all(x[i]==a for i,a in zip(subset,residues[choice]))
                value+=max(0,count-j)
            assert value==expected[j]
    result=dict(passed=True,n_configurations=total,n_core_patterns=3,
                uncovered_points=27,stoploss_numerators=expected,denominator=27,
                law_numerators=law,count_seed_numerators=[6,13,7,1],
                argmax_indices=args.tolist(),scope='complete integer core comparison; no floating solver')
    Path(output).write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':audit(sys.argv[1])
