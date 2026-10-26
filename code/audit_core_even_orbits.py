"""Check canonical-even and full-even enumeration give the same core orbits.

This audits the first-use label reduction, not the separate inclusion
reduction or the conditional-tail inequalities. NumPy is required.
"""
import itertools
import json
from pathlib import Path
import time
from core_orbit_reduction import CoreOrbits
from core_residual_geometry import pair_residuals


def main():
    start = time.monotonic()
    results = []
    for suffix in itertools.product((0,1),repeat=2):
        p = (0,)+suffix
        model = CoreOrbits(p)
        _,first,_ = pair_residuals(p)
        expected = set(model.canonical(first))
        _,full,_ = pair_residuals(p,False)
        observed = set(model.canonical(full))
        assert expected==observed
        row = dict(pattern=list(p),canonical_source_masks=len(first),full_even_source_masks=len(full),
                   orbits=len(expected),same_orbits=True)
        results.append(row)
        print(json.dumps(row|{'seconds':round(time.monotonic()-start,3)}),flush=True)
    result = dict(verified=True,cases=results,seconds=round(time.monotonic()-start,3),
                  scope='full-even versus canonical-even pair residuals have identical CRT-label orbits')
    Path('data/audit-core-even-orbits.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
