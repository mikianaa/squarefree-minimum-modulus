"""Independent arithmetic audit of the seven-prime plane certificates.

Uses clipped-capacity hinge functions and polynomial convolution, rather than
balanced load vectors and the three-index loss tensor in the producer. Tail
capacities are reconstructed by enumerating ordinary divisors. Allocation
knapsacks use a Python dictionary, independently of the producer's array DP.
"""
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np

DS = (1,3,5,7,15,21,35,105)
HIGH = (11,13,17)
DIVISORS = sorted(math.prod(p for i,p in enumerate((2,3,5,7,11,13,17)) if mask>>i&1) for mask in range(128))


def arrays(profile,pattern,overrides=None):
    result = []
    for e in (0,1):
        target = profile[e]
        if target==0:
            result.append(np.zeros((128,)*3,dtype=np.int64))
            continue
        axis_hist = []
        axis_hinge = []
        for axis,p in enumerate(HIGH):
            size = p-1-int(pattern[axis]==e)
            histograms = []
            hinges = []
            for mask in range(128):
                odd = [profile[2+2*j+e] for j in range(7)]
                if overrides is not None:
                    odd[0] = overrides[axis][e]
                weights = odd+[profile[2+2*j+e] for j in range(7) if ((mask>>(6-j))&1)==e]
                h = [max(0,size*t-sum(min(t,w) for w in weights)) for t in range(target+2)]
                hist = [h[1]]+[h[t+1]-2*h[t]+h[t-1] for t in range(1,target)]
                assert min(hist)>=0 and sum(hist)<=size
                histograms.append(hist)
                hinges.append(h)
            axis_hist.append(np.array(histograms,dtype=np.int64))
            axis_hinge.append(np.array(hinges,dtype=np.int64))
        # Truncated polynomial product of the first two integer distributions.
        convolution = np.zeros((128,128,target),dtype=np.int64)
        for i in range(target):
            for j in range(target-i):
                convolution[:,:,i+j] += np.outer(axis_hist[0][:,i],axis_hist[1][:,j])
        third = axis_hinge[2][:,target-np.arange(target)]
        value = convolution.reshape(128*128,target) @ third.T
        result.append(value.reshape((128,)*3))
    return result


def tail(profile,pattern):
    odd = [0,0]
    even = []
    for d in DIVISORS:
        if sum(d%p==0 for p in HIGH)<2:
            continue
        core = math.gcd(d,210)
        if core%2==0:
            core //= 2
        j = DS.index(core)
        pair = []
        for e in (0,1):
            cap = profile[e] if j==0 else profile[2*j+e]
            pair.append(cap*math.prod(p-1-int(pattern[i]==e) for i,p in enumerate(HIGH) if d%p))
        if d%2:
            odd = [x+y for x,y in zip(odd,pair)]
        else:
            even.append(tuple(pair))
    assert len(even)==32
    return odd,even


def check(profile,pattern,certificate,overrides=None,refined=None):
    values = arrays(profile,pattern,overrides)
    lower = int((values[0]+values[1]).min())
    odd,even = tail(profile,pattern)
    upper = sum(odd)+sum(max(c) for c in even)
    assert lower==certificate['plane_deficit_lower']
    assert upper==certificate['line_point_capacity_upper']
    assert lower-upper==certificate['gap']
    assert (lower>upper)==certificate['excluded']
    if refined is not None:
        best = {}
        for x,y in zip(values[0].flat,values[1].flat):
            x,y = int(x),int(y)
            if y<best.get(x,10**12):
                best[x] = y
        frontier = []
        previous = 10**12
        for x,y in sorted(best.items()):
            if y<previous:
                frontier.append([x,y])
                previous = y
        assert frontier==refined['frontier']
        states = {0:0}
        for gain,cost in even:
            updated = dict(states)
            for used,value in states.items():
                updated[used+cost] = max(updated.get(used+cost,-1),value+gain)
            states = updated
        survivors = []
        for x,y in frontier:
            need = max(0,x-odd[0])
            budget = sum(c[1] for c in even)-max(0,y-odd[1])
            if budget>=0 and any(used<=budget and value>=need for used,value in states.items()):
                survivors.append([x,y])
        assert survivors==refined['surviving_lower_pairs']
        assert (not survivors)==refined['excluded']
    return True


def main():
    start = time.monotonic()
    profiles = json.loads(Path('data/core-geometry-profile-audit.json').read_text())
    base = json.loads(Path('data/plane-capacity-510510.json').read_text())
    refined = json.loads(Path('data/refined-plane-capacity-510510.json').read_text())
    lookup = {(tuple(x['core_pattern']),tuple(x['high_pattern']),x['profile_index']):x for x in refined}
    checked = 0
    summaries = []
    assert len(base)==32
    for row in base:
        source = next(x for x in profiles if x['pattern']==row['core_pattern'])
        assert len(source['easiest_profiles'])==len(row['tests'])
        for i,(p,c) in enumerate(zip(source['easiest_profiles'],row['tests'])):
            check(p,row['high_pattern'],c,refined=lookup.get((tuple(row['core_pattern']),tuple(row['high_pattern']),i)))
            checked += 1
        summaries.append(dict(core_pattern=row['core_pattern'],high_pattern=row['high_pattern'],profiles=len(row['tests']),verified=True))
        print(json.dumps(summaries[-1]),flush=True)
    extra = 0
    for suffix in ('010','011'):
        data = json.loads(Path(f'data/core-triple-refinement-{suffix}.json').read_text())
        for t in data['tests']:
            check(t['profile'],[1,1,1],t['simple'],refined=t.get('refined'))
            extra += 1
    coupled = json.loads(Path('data/coupled-core3-510510.json').read_text())
    for t in coupled['tests']:
        p = t['profile']
        u,v = t['layer1_mod3_counts']
        for c in t['cases']:
            overrides = [(0,u) if choice==1 else (p[0],v) for choice in c['core_mod3_choices']]
            check(p,[1,1,1],c,overrides=overrides)
            assert c['excluded']
            extra += 1
    result = dict(base_profiles_checked=checked,additional_profiles_checked=extra,
                  cases=summaries,verified=True,seconds=round(time.monotonic()-start,3),
                  scope='independent arithmetic audit; core-catalogue completeness is audited separately')
    Path('data/audit-plane-capacity-510510.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)


if __name__=='__main__':
    main()
