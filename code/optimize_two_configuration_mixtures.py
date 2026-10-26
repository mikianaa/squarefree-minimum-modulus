"""Optimize independent group probabilities using point-LP supporting cuts.

Reported dual bounds are exploratory floating-point Benders bounds. Saved
primal configurations and their rational costs are independently audited.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import time
from ortools.linear_solver import pywraplp
from two_configuration_mixtures import structure,counts_for,point_coupling,evaluate,assignments_from_integer


def optimize(source,A,B,seconds=60,max_iterations=500,verbose=False):
    start = time.monotonic()
    groups,shapes = structure(source)
    n,G,R = len(source['high_primes']),len(groups),source['core_residual']
    ca,cb = counts_for(source,A),counts_for(source,B)
    solver = pywraplp.Solver.CreateSolver('GLOP')
    theta = [solver.NumVar(0,int(any(a[h]!=b[h] for a,b in zip(ca,cb))),f't{h}') for h in range(G)]
    costs = []
    for x in R:
        D = [q-1-(b==x%2) for q,b in zip(source['high_primes'],source['high_pattern'])]
        lower = -16*sum(math.prod(D[i] for i in range(n) if i not in S) for S in groups[n:])
        costs.append(solver.NumVar(lower,solver.infinity(),f'c{x}'))
    solver.Minimize(sum(costs))
    current = [0.0]*G
    best_value,best_theta = math.inf,None
    lower = -math.inf
    history = []
    for iteration in range(max_iterations):
        solutions = [point_coupling(source,x,a,b,current) for x,a,b in zip(R,ca,cb)]
        value = sum(r['objective'] for r in solutions)
        if value<best_value:
            best_value,best_theta = value,current[:]
        for z,row in zip(costs,solutions):
            solver.Add(z>=row['supporting_intercept']+sum(g*t for g,t in zip(row['gradient'],theta)))
        solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
        if solver.Solve()!=pywraplp.Solver.OPTIMAL:
            break
        lower = solver.Objective().Value()
        history.append(dict(iteration=iteration+1,upper=best_value,lower=lower,seconds=time.monotonic()-start))
        if verbose and (iteration<5 or (iteration+1)%10==0):
            print(json.dumps(history[-1]),flush=True)
        if best_value-lower<1e-4 or time.monotonic()-start>=seconds:
            break
        current = [max(0,min(1,t.solution_value())) for t in theta]
    rounded = [Fraction(round(p*10**9),10**9) for p in best_theta]
    result = evaluate(source,A,B,rounded,save_witness=True)
    result.update(assignments_A=A,assignments_B=B,probabilities=[[p.numerator,p.denominator] for p in rounded],
                  exploratory_lower=lower,numerical_gap=best_value-lower,history=history,
                  iterations=len(history),seconds=time.monotonic()-start,
                  scope='optimized probabilities for one fixed pair per group; no global optimality or obstruction claim')
    return result


def choose_alternatives(source,pool):
    R = source['core_residual']
    groups,shapes = structure(source)
    prices = pool.get('certificate',{}).get('event_multipliers')
    reports = {r['axis']:r for r in pool.get('group_reports',[])}
    selected = []
    for h,rows in enumerate(pool['group_columns']):
        events = reports.get(h,{}).get('candidate_cut',{}).get('events')
        def score(residues):
            counts = [sum(x%d==a for d,a in zip(shapes[h],residues)) for x in R]
            if events:
                return sum(e['weight'] for e in events if counts[e['point']]>=e['threshold'])
            if prices:
                return -sum(sum(prices[h][x][1:v+1]) for x,v in enumerate(counts))
            return sum(counts)
        selected.append(max(rows,key=score))
    return selected


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--integer',required=True)
    ap.add_argument('--pool',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=60)
    args = ap.parse_args()
    data,integer,pool = [json.loads(Path(p).read_text()) for p in (args.source,args.integer,args.pool)]
    for key in ('core_residual','high_primes','high_pattern'):
        assert data['input'][key]==integer['input'][key]==pool['input'][key]
    A = assignments_from_integer(data['input'],integer['result'])
    B = choose_alternatives(data['input'],pool['result'])
    result = optimize(data['input'],A,B,args.time,verbose=True)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(seconds=result['seconds'],iterations=result['iterations'],gap=result['numerical_gap'],audit=result['audit'])))
