"""Choose additional nonnegative Lagrange coefficients from prior certificates.

Past slopes are heuristic choices only. New lower bounds must be recomputed
with the current penalties and independently audited with the saved coefficients.
"""
from collections import Counter
from fractions import Fraction


def coefficients(default,planes,old_denominator,new_denominator,limit=12,maximum=None,prior=(),reward_ratio=None):
    assert type(limit) is int and limit>=0
    assert old_denominator>0 and new_denominator>0
    assert default and all(type(v) is int and v>=0 for v in default)
    base=set(default)
    if maximum is not None:assert max(base)<=maximum
    if not limit:return sorted(base)
    scale=Fraction(new_denominator,old_denominator) if reward_ratio is None else Fraction(*reward_ratio)
    assert scale>=0
    frequencies=Counter()
    for row in planes:
        num,den=row['slope'];assert type(num) is int and num>=0 and type(den) is int and den>0
        value=Fraction(num,den)*scale
        # Nearby integer coefficients give valid bounds independently of whether
        # either is an optimal slope for the new penalties.
        for v in {value.numerator//value.denominator,-(-value.numerator//value.denominator)}:
            if v not in base and (maximum is None or v<=maximum):frequencies[v]+=1
    for v in prior:
        assert type(v) is int and v>=0
        value=v*scale
        rounded=-(-value.numerator//value.denominator)
        if rounded not in base and (maximum is None or rounded<=maximum):frequencies[rounded]+=1
    if len(frequencies)<=limit:return sorted(base|set(frequencies))
    # Keep frequent critical slopes and a spread across the remaining range.
    chosen=set(v for v,count in sorted(frequencies.items(),key=lambda item:(-item[1],item[0]))[:(limit+1)//2])
    others=sorted(set(frequencies)-chosen);needed=limit-len(chosen)
    for j in range(needed):chosen.add(others[j*(len(others)-1)//max(1,needed-1)])
    return sorted(base|chosen)


def safe_maximum(model,parity,potentials):
    """Match the conservative int64 bound used by positive_state_dual."""
    fixed=sum(max(map(abs,row)) for row in potentials)
    factor=int(model.costs[parity].max())+16*sum(int(v.max()) for v in model.volumes[parity])
    assert factor>0 and fixed<2**63
    return (2**63-1-fixed)//factor
