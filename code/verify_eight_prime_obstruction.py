"""Re-run all four branch audits and summarize the eight-prime obstruction.

The summary is written only when all four independently checked branches
have positive bounds. --record-existing only summarizes already saved audit
reports; it does not itself reverify their mathematical assertions.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def record_summary():
    names = ['audit-eight-prime-anchor000.json','audit-compressed-tail-branch-001.json',
             'audit-paired-tail-branch-010.json','audit-paired-tail-branch-011.json']
    audits = [json.loads(Path('data',name).read_text()) for name in names]
    patterns = [[0,*p] for p in itertools.product((0,1),repeat=2)]
    bounds = []
    for pattern,report in zip(patterns,audits):
        assert report['verified'] and report['period']==9699690 and report['core_pattern']==pattern
        lower = report.get('uniform_uncovered_lower',report.get('uniform_lower'))
        assert type(lower) is int and lower>0
        bounds.append(lower)
    files = names+['eight-prime-anchor000-certificate.json','conditional-tail-core-pattern-certificates.json',
                   'core-hit-concentration-7-15.json','core-fibre-eight-plane.json',
                   'compressed-tail-branch-001.json']
    for tag in ('010','011'):
        files.extend([f'paired-tail-refinement-{tag}.json',f'joint-tail-profiles-{tag}.json',f'joint-tail-envelopes-{tag}.json'])
    hashes = {name:hashlib.sha256(Path('data',name).read_bytes()).hexdigest() for name in files}
    result = dict(excluded_primorial=9699690,excluded_prime_support=[2,3,5,7,11,13,17,19],
                  minimum_allowed_modulus=3,core_patterns=patterns,branch_lower_bounds=bounds,
                  uniform_uncovered_lower=min(bounds),counterexample_prime_factors_at_least=9,
                  next_primorial=223092870,solver_dependent=False,proof_assistant_verified=False,
                  global_problem_settled=False,
                  extension_to_arbitrary_support='BBMST even-modulus theorem plus odd-prime compression',
                  certificate_sha256=hashes,
                  scope='summary of the four saved independent branch audits; not a proof of the global minimum-2 conjecture')
    Path('data/eight-prime-obstruction.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--record-existing',action='store_true')
    args = ap.parse_args()
    if not args.record_existing:
        from audit_eight_prime_anchor_branch import main as audit_000
        from audit_compressed_tail_branch import audit as audit_001
        from audit_paired_tail_branch import audit as audit_paired
        audit_000()
        result = audit_001(json.loads(Path('data/compressed-tail-branch-001.json').read_text()))
        Path('data/audit-compressed-tail-branch-001.json').write_text(json.dumps(result,indent=2)+'\n')
        for pattern in ([0,1,0],[0,1,1]):
            result = audit_paired(pattern)
            tag = ''.join(map(str,pattern))
            Path(f'data/audit-paired-tail-branch-{tag}.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps(record_summary(),indent=2))


if __name__=='__main__':
    main()
