"""Independent integer audit of conditional reduction and power-of-two seed.

Rebuilds residual geometry, groups with explicit modular indexing, and core-30
catalogues by mixed-radix decoding; checks all seven necessary thresholds.
"""
import json,sys,math,time,hashlib
from fractions import Fraction
from pathlib import Path
import numpy as np


def audit(source,output):
    t=time.time();raw=Path(source).read_bytes();data=json.loads(raw)
    measure=json.loads(Path(data['measure_file']).read_text());family={int(d):r for d,r in measure['family'].items()}
    assert set(family)=={d for d in range(7,2311) if 2310%d==0}
    fixed=family
    if measure.get('support_mode')=='fixed_subfamily':
        free=set(measure['free_moduli']);assert free==set(data['initial_free_moduli']) and free<=set(family)
        fixed={d:r for d,r in family.items() if d not in free}
        assert fixed=={int(d):r for d,r in measure['fixed_family'].items()}
    residual={x for x in range(2310) if all(x%d!=r for d,r in fixed.items())}
    assert residual==set(measure['support'])
    w=data['weights_numerators'];den=data['weights_denominator']
    assert len(w)==2310 and sum(w)==den and den>0 and all(isinstance(v,int) and v>=0 for v in w)
    assert all(v==0 or x in residual for x,v in enumerate(w))
    groups=[{1:[4*v for v in w]}];period=2310
    for p in data['elimination_order']:
        nxt=[];n=period//p
        for group in groups:
            left={};right={}
            for k,v in group.items():
                left.setdefault(k,[0]*n);left.setdefault(2*k,[0]*n);right.setdefault(2*k,[0]*n)
                for x in range(n):
                    vals=[v[y] for y in range(x,period,n)];cap=max(vals);assert cap%2==0
                    left[k][x]+=sum(vals)-cap;left[2*k][x]+=cap//2;right[2*k][x]+=cap//2
            nxt.extend([left,right])
        groups=nxt;period=n
    assert period==30
    knots=[0,1,2,4,8,16,32];maxima=np.zeros((4,len(knots)),dtype=np.int64);divs=[2,3,5,6,10,15,30]
    assert 32*4*den<2**63
    for lo in range(0,math.prod(divs),4096):
        indices=np.arange(lo,min(lo+4096,math.prod(divs)));code=indices.copy();counts=np.ones((len(indices),30),dtype=np.int64)
        for d in reversed(divs):
            residues=code%d;code//=d;counts+=(np.arange(30)[None,:]%d==residues[:,None])
        assert np.max(code)==0
        for g,group in enumerate(groups):
            for j,threshold in enumerate(knots):
                values=sum(np.maximum(k*counts-threshold,0)@np.array(v,dtype=np.int64) for k,v in group.items())
                maxima[g,j]=max(maxima[g,j],int(values.max()))
    H=[int(x) for x in maxima.sum(axis=0)]
    assert H==[data['envelope_numerators'][j] for j in knots]
    assert H[0]-H[1]==4*den and H[-1]==0
    slopes=[Fraction(H[i+1]-H[i],knots[i+1]-knots[i]) for i in range(6)]
    masses=[slopes[i+1]-slopes[i] for i in range(5)]+[-slopes[-1]]
    assert all(m>=0 for m in masses) and sum(masses)==4*den
    assert [m/(4*den) for m in masses]==[Fraction(n,data['power_seed_denominator']) for n in data['power_seed_numerators']]
    result=dict(status='VERIFIED',claim='The specified probability measure on the fixed residual core has the claimed convex power-of-two comparison',
                core_file=source,core_sha256=hashlib.sha256(raw).hexdigest(),catalogues=math.prod(divs),thresholds=knots,
                residual=len(residual),seconds=time.time()-t,scope='One fixed subfamily; later continuation and any free initial classes are charged separately')
    Path(output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    audit(*sys.argv[1:])
