"""SAT encodings and RUP proof audits for weighted core event bounds.

Ask whether one residue per resource can make the weighted threshold score
exceed a proposed bound. Threshold flags imply sufficient actual hits. Fresh
copies of flags encode integer weights in an ordinary cardinality constraint.
An UNSAT answer is accepted only after its proof passes separate propagation.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import threading
import time
from pysat.card import CardEnc,EncType
from pysat.solvers import Solver
from sat_uncovered_audit import check_rup


def build(R,shapes,events,bound):
    clauses,variables,top = [],{},0
    for d in shapes:
        literals = []
        for a in sorted({v%d for v in R}):
            top += 1
            variables[d,a] = top
            literals.append(top)
        clauses.append(literals)
        enc = CardEnc.atmost(literals,bound=1,top_id=top,encoding=EncType.seqcounter)
        clauses.extend(enc.clauses)
        top = max(top,enc.nv)
    weighted = []
    for row in events:
        x,t,w = (row[k] for k in ('point','threshold','weight'))
        assert type(x) is int and 0<=x<len(R)
        assert type(t) is int and 1<=t<=len(shapes)
        assert type(w) is int and w>0
        top += 1
        flag = top
        hits = [variables[d,R[x]%d] for d in shapes]
        enc = CardEnc.atleast(hits,bound=t,top_id=top,encoding=EncType.seqcounter)
        clauses.extend([-flag]+clause for clause in enc.clauses)
        top = max(top,enc.nv)
        for i in range(w):
            top += 1
            clauses.extend(([-top,flag],[top,-flag]))
            weighted.append(top)
    enc = CardEnc.atleast(weighted,bound=bound+1,top_id=top,encoding=EncType.seqcounter)
    clauses.extend(enc.clauses)
    return clauses,max(top,enc.nv)


def verify_bound(R,shapes,cut,proof_path,seconds=60,existing=False):
    start = time.monotonic()
    clauses,variables = build(R,shapes,cut['events'],cut['bound'])
    cnf_hash = hashlib.sha256(json.dumps(clauses,separators=(',',':')).encode()).hexdigest()
    if existing:
        with gzip.open(proof_path,'rt') as f:
            proof = f.read().splitlines()
        answer = False
    else:
        with Solver(name='glucose4',bootstrap_with=clauses,with_proof=True) as solver:
            timer = threading.Timer(seconds,solver.interrupt)
            timer.daemon = True
            timer.start()
            try:
                answer = solver.solve_limited(expect_interrupt=True)
            finally:
                timer.cancel()
            proof = solver.get_proof() if answer is False else None
    result = dict(axis=cut['axis'],bound=cut['bound'],variables=variables,clauses=len(clauses),cnf_sha256=cnf_hash,
                  status='UNSAT' if answer is False else 'SAT' if answer is True else 'UNKNOWN')
    if proof is not None:
        result['rup_audit'] = check_rup(clauses,proof)
        raw = ('\n'.join(proof)+'\n').encode()
        if not existing:
            with gzip.open(proof_path,'wb') as f:
                f.write(raw)
        result.update(proof_path=str(proof_path),proof_sha256=hashlib.sha256(raw).hexdigest())
    result['verified'] = result.get('rup_audit',{}).get('verified',False)
    result['seconds'] = time.monotonic()-start
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=60)
    ap.add_argument('--existing',action='store_true')
    args = ap.parse_args()
    data = json.loads(Path(args.source).read_text())
    base = json.loads(Path(data['source']).read_text())
    R,n = base['input']['core_residual'],len(base['input']['high_primes'])
    reports = []
    for cut in data['cuts']:
        shapes = [d for d in range(1 if cut['axis']>=n else 3,211) if 210%d==0]
        proof = Path(args.output).with_name(Path(args.output).stem+f"-axis{cut['axis']}.rup.gz")
        report = verify_bound(R,shapes,cut,proof,args.time,args.existing)
        reports.append(report)
        Path(args.output).write_text(json.dumps(dict(source=args.source,reports=reports,
                                                   all_verified=len(reports)==len(data['cuts']) and all(r['verified'] for r in reports)),indent=2)+'\n')
        print(json.dumps(report),flush=True)
