"""Independent int64 table audit of every catalogue in a core batch.

Uses mixed-radix catalogue generation, explicit modular-index reduction,
and direct pointwise lookup/addition instead of the generator's dot products.
"""
import argparse,hashlib,json,math,time
from fractions import Fraction
from pathlib import Path
import numpy as np


def audit(manifest_file):
    started=time.time();manifest=json.loads(Path(manifest_file).read_text());tables=[];cores=[]
    for entry in manifest:
        raw=Path(entry['core_file']).read_bytes();core=json.loads(raw);data=json.loads(Path(core['measure_file']).read_text())
        family={int(d):r for d,r in data['family'].items()};free=set(core['initial_free_moduli']);assert free==set(data['free_moduli'])
        assert set(family)=={d for d in range(7,2311) if 2310%d==0} and free<=set(family)
        fixed={d:r for d,r in family.items() if d not in free};assert fixed=={int(d):r for d,r in data['fixed_family'].items()}
        support={x for x in range(2310) if all(x%d!=r for d,r in fixed.items())};assert support==set(data['support'])
        weights=core['weights_numerators'];den=core['weights_denominator'];assert len(weights)==2310 and sum(weights)==den>0
        assert all(isinstance(v,int) and v>=0 and (v==0 or x in support) for x,v in enumerate(weights))
        groups=[{1:[4*v for v in weights]}];period=2310
        for p in [11,7]:
            out=[];n=period//p
            for group in groups:
                left={};right={}
                for k,v in group.items():
                    left.setdefault(k,[0]*n);left.setdefault(2*k,[0]*n);right.setdefault(2*k,[0]*n)
                    for x in range(n):
                        fibre=[v[y] for y in range(x,period,n)];cap=max(fibre);assert cap%2==0
                        left[k][x]+=sum(fibre)-cap;left[2*k][x]+=cap//2;right[2*k][x]+=cap//2
                out.extend([left,right])
            groups=out;period=n
        assert 32*4*den<2**63
        for group in groups:
            for j in [1,2,4,8,16]:
                table=np.array([[sum(v[x]*max(k*m-j,0) for k,v in group.items()) for m in range(9)] for x in range(30)],dtype=np.int64)
                tables.append(table)
        cores.append((core,hashlib.sha256(raw).hexdigest()))
    T=np.stack(tables,axis=-1)
    # Each summand is nonnegative and increasing in the catalogue count.
    # This bound certifies every partial sum, permitting narrower exact
    # arithmetic and keeping the working set in cache.
    bound=T[:,8,:].sum(axis=0);dtype=np.uint32 if int(max(bound))<2**32 else np.int64
    assert int(max(bound))<2**63;T=T.astype(dtype)
    maximum=np.zeros(len(tables),dtype=dtype);divs=[2,3,5,6,10,15,30];chunk=96
    for lo in range(0,math.prod(divs),chunk):
        indices=np.arange(lo,min(lo+chunk,math.prod(divs)));code=indices.copy();counts=np.ones((len(indices),30),dtype=np.int8)
        for d in reversed(divs):residues=code%d;code//=d;counts+=(np.arange(30)[None,:]%d==residues[:,None])
        assert max(code)==0
        values=np.zeros((len(indices),len(tables)),dtype=dtype)
        for x in range(30):values+=T[x,counts[:,x],:]
        maximum=np.maximum(maximum,values.max(axis=0))
    for i,(core,sha) in enumerate(cores):
        maxima=maximum[20*i:20*(i+1)].reshape(4,5);assert maxima.tolist()==core['group_power_envelope_numerators']
        h=[int(v) for v in maxima.sum(axis=0)];H=[h[0]+4*core['weights_denominator']]+h+[0];assert H==core['power_envelope_numerators']
        knots=[0,1,2,4,8,16,32];slopes=[Fraction(H[j+1]-H[j],knots[j+1]-knots[j]) for j in range(6)]
        masses=[slopes[j+1]-slopes[j] for j in range(5)]+[-slopes[-1]];den=4*core['weights_denominator']
        assert all(v>=0 for v in masses) and sum(masses)==den
        assert [v/den for v in masses]==[Fraction(v,core['power_seed_denominator']) for v in core['power_seed_numerators']]
        source=Path(manifest[i]['core_file']);out=source.with_name('audit-'+source.name)
        result=dict(status='VERIFIED',core_file=str(source),core_sha256=sha,catalogues=810000,thresholds=knots,
                    arithmetic=str(np.dtype(dtype)),partial_sum_bound=int(max(bound)),
                    scope='Independent overflow-checked integer core audit; initial free-class charges and infinite continuation are separate')
        out.write_text(json.dumps(result,indent=2))
    print('BATCH INDEPENDENT INTEGER AUDITS',len(cores),'seconds',time.time()-started,flush=True)


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('manifest');a=ap.parse_args();audit(a.manifest)
