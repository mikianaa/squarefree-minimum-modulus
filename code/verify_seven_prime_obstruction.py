"""Re-run both independent audits and record the finite obstruction summary."""
import hashlib
import json
from pathlib import Path
import audit_core_completeness
import audit_plane_capacity


def record_summary():
    core = json.loads(Path('data/audit-core-completeness-510510.json').read_text())
    plane = json.loads(Path('data/audit-plane-capacity-510510.json').read_text())
    assert core['verified'] and plane['verified']
    assert len(core['cases'])==4 and len(plane['cases'])==32
    assert plane['base_profiles_checked']==2936 and plane['additional_profiles_checked']==116
    files = ['core-geometry-profile-audit.json','plane-capacity-510510.json',
             'refined-plane-capacity-510510.json','core-triple-refinement-010.json',
             'core-triple-refinement-011.json','coupled-core3-510510.json',
             'audit-core-completeness-510510.json','audit-plane-capacity-510510.json']
    hashes = {f:hashlib.sha256((Path('data')/f).read_bytes()).hexdigest() for f in files}
    result = dict(excluded_primorial=510510,excluded_prime_support=[2,3,5,7,11,13,17],
                  minimum_allowed_modulus=3,anchor_patterns_excluded=32,
                  counterexample_prime_factors_at_least=8,next_primorial=9699690,
                  arithmetic='exact integer enumeration and two independent arithmetic implementations',
                  solver_dependent=False,proof_assistant_verified=False,global_problem_settled=False,
                  extension_to_arbitrary_support='BBMST even-modulus theorem plus prime compression',
                  certificate_sha256=hashes)
    Path('data/seven-prime-obstruction.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    audit_core_completeness.main()
    audit_plane_capacity.main()
    print(json.dumps(record_summary(),indent=2))
