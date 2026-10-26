"""Independent integer-residue SAT encoding of a fixed-pattern deficit bound.

Ask whether <= limit residues can remain uncovered. Uses literal residues modulo
N, cardinality networks, and optional first-use label ordering. A Glucose RUP
proof, when produced, is checked clause by clause with Minisat's unit propagation.
This checks the propositional proof, not the mathematical normalization theorem.
"""
import argparse
import gzip
import hashlib
import json
import math
import threading
import time
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


def build(primes,pattern,limit,symmetry=True):
    n = 2*math.prod(primes)
    anchors = [(p,0) for p in primes]
    anchors += [(2*p,1 if b else p+1) for p,b in zip(primes,pattern)]
    remaining = [a for a in range(n) if not any(a%d == r for d,r in anchors)]
    anchor_moduli = {d for d,r in anchors}
    ds = [d for d in range(3,n+1) if n%d == 0 and d not in anchor_moduli]
    variables = {}
    clauses = []
    top = 0
    previous = {p:{a:[] for a in range(1,p)} for p in primes}
    for d in ds:
        residues = sorted({a%d for a in remaining})
        lits = []
        for r in residues:
            top += 1
            variables[d,r] = top
            lits.append(top)
            if symmetry:
                for p in primes:
                    if d%p == 0 and r%p >= 3:
                        clauses.append([-top]+previous[p][r%p-1])
        clauses.append(lits)
        enc = CardEnc.atmost(lits,bound=1,top_id=top,encoding=EncType.seqcounter)
        clauses.extend(enc.clauses)
        top = max(top,enc.nv)
        for r in residues:
            for p in primes:
                if d%p == 0:
                    previous[p][r%p].append(variables[d,r])
    uncovered = []
    for a in remaining:
        top += 1
        uncovered.append(top)
        clauses.append([variables[d,a%d] for d in ds]+[top])
    enc = CardEnc.atmost(uncovered,bound=limit,top_id=top,encoding=EncType.cardnetwrk)
    clauses.extend(enc.clauses)
    return clauses,max(top,enc.nv)


def check_rup(clauses,proof):
    checked = 0
    with Solver(name='minisat22',bootstrap_with=clauses) as checker:
        for line in proof:
            if line.startswith('d '):
                continue  # Retaining already proved consequences is sound.
            clause = [int(v) for v in line.split() if int(v)]
            consistent,_ = checker.propagate(assumptions=[-v for v in clause])
            if consistent:
                return dict(verified=False,checked_additions=checked,reason='not RUP; may require RAT verification')
            checker.add_clause(clause)
            checked += 1
            if checked % 10000 == 0:
                print(json.dumps(dict(rup_additions_checked=checked)),flush=True)
        consistent,_ = checker.propagate()
        return dict(verified=not consistent,checked_additions=checked,
                    reason='root contradiction after verified RUP additions' if not consistent else 'no final contradiction')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--primes',type=int,nargs='+',default=[3,5,7,11])
    ap.add_argument('--pattern',required=True)
    ap.add_argument('--limit',type=int,required=True)
    ap.add_argument('--time',type=float,default=120)
    ap.add_argument('--output',required=True)
    ap.add_argument('--no-symmetry',action='store_true')
    ap.add_argument('--solver',default='glucose4')
    ap.add_argument('--skip-rup-check',action='store_true')
    a = ap.parse_args()
    start = time.monotonic()
    pattern = list(map(int,a.pattern))
    clauses,nvars = build(a.primes,pattern,a.limit,not a.no_symmetry)
    with Solver(name=a.solver,bootstrap_with=clauses,with_proof=True) as solver:
        if a.solver.startswith('cadical'):
            deadline = time.monotonic()+a.time
            answer = None
            while answer is None and time.monotonic() < deadline:
                solver.conf_budget(2000)
                answer = solver.solve_limited()
        else:
            timer = threading.Timer(a.time,solver.interrupt)
            timer.daemon = True
            timer.start()
            try:
                answer = solver.solve_limited(expect_interrupt=True)
            finally:
                timer.cancel()
        proof = solver.get_proof() if answer is False else None
        stats = solver.accum_stats()
    result = dict(N=2*math.prod(a.primes),pattern=pattern,uncovered_at_most=a.limit,
                  primes=a.primes,label_symmetry=not a.no_symmetry,solver=a.solver,
                  status='UNSAT' if answer is False else 'SAT' if answer is True else 'UNKNOWN',
                  variables=nvars,clauses=len(clauses),seconds=round(time.monotonic()-start,3),stats=stats)
    print(json.dumps(result),flush=True)
    if proof is not None:
        result['rup_audit'] = (dict(verified=False,reason='not requested') if a.skip_rup_check else check_rup(clauses,proof))
        raw = ('\n'.join(proof)+'\n').encode()
        path = Path(a.output).with_suffix('.rup.gz')
        with gzip.open(path,'wb') as f:
            f.write(raw)
        result['proof_path'] = str(path)
        result['proof_sha256'] = hashlib.sha256(raw).hexdigest()
        print(json.dumps(result['rup_audit']),flush=True)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
