"""Conditional-tail bounds for the entire low-anchor pattern 000 at 9699690."""
import json
from pathlib import Path
from portable_tail_certificates import CoreCapacityUpper,PortableDuals


def main():
    read = lambda name: json.loads(Path('data',name).read_text())
    base = next(r for r in read('three-axis-210-certificate.json')['cases'] if r['pattern']==[0,0,0])
    model = CoreCapacityUpper([0,0,0])
    duals = PortableDuals(read('conditional-tail-core-pattern-certificates.json'))
    rows = [dict(core_counts=[a,b],bounds=duals.bounds(model.compression[:a+1,:b+1]).tolist()) for a,b in base['residual_pair_frontier']]
    result = dict(period=9699690,core_pattern=[0,0,0],high_patterns=duals.patterns,
                  frontier=base['residual_pair_frontier'],cases=rows,
                  uniform_lower=min(min(r['bounds']) for r in rows),
                  scope='all normalized 210 cores with anchor parity pattern 000, and all 16 higher patterns')
    Path('data/eight-prime-anchor000-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}))


if __name__=='__main__':
    main()
