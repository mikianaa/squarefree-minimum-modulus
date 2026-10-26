"""Benders decomposition of the explicit coherent-group relaxation.

Point distributions use integer-priced columns; a small-row master mixes
actual group configurations. Finite-pool lower estimates are not global
obstructions. Rational primal witnesses are checked separately.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from ortools.linear_solver import pywraplp
from coherent_transport_repair import rational_witness,drop_cost_bounds
from audit_coherent_groups import audit
from multitail_lp import LocalStates
from group_event_coherence import maximize
from two_configuration_mixtures import structure,assignments_from_integer
from coherent_point_transport import solve_point as compact_point
from event_score_majorant import upper_bound as event_majorant
from coherent_group_catalogue import Catalogue


class FixedValue:
    def __init__(self,value):self.value=value
    def solution_value(self):return self.value


class PointOracle:
    def __init__(self,source,model,index,seed,bounded_duals=True):
        self.index,self.parity = index,source['core_residual'][index]%2
        self.model = model
        self.shapes = structure(source)[1]
        self.solver = pywraplp.Solver.CreateSolver('GLOP')
        self.solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
        self.mass = self.solver.Constraint(1,1)
        self.links = [[None]+[self.solver.Constraint(-self.solver.infinity(),0) for d in shape]
                      for shape in self.shapes]
        self.objective = self.solver.Objective()
        self.objective.SetMinimization()
        self.cost_scale = max(1,model.cost(self.parity,(0,)*len(self.shapes)))
        self.numerical_retries = 0
        self.slacks = []
        self.bounds = drop_cost_bounds(model.high,model.pattern,self.parity)
        if bounded_duals:
            # Paying B_h for one missing threshold unit is exact: decreasing
            # that count by one raises F by at most B_h. Quantile repair
            # realizes this operation for arbitrary marginal distributions.
            for h,line in enumerate(self.links):
                for t,c in enumerate(line[1:],1):
                    slack = self.solver.NumVar(0,self.solver.infinity(),f'slack{h}_{t}')
                    c.SetCoefficient(slack,-1)
                    self.objective.SetCoefficient(slack,self.bounds[h]/self.cost_scale)
                    self.slacks.append((h,t,slack))
        self.columns,self.seen = [],set()
        self.add((0,)*len(self.shapes))
        for row in seed:
            if row['point_index']==index:
                self.add(row['counts'])

    def add(self,state):
        state = tuple(state)
        if state in self.seen:
            return 0
        self.seen.add(state)
        var = self.solver.NumVar(0,self.solver.infinity(),str(len(self.columns)))
        self.mass.SetCoefficient(var,1)
        self.objective.SetCoefficient(var,self.model.cost(self.parity,state)/self.cost_scale)
        for h,v in enumerate(state):
            for t in range(1,v+1):
                self.links[h][t].SetCoefficient(var,1)
        self.columns.append((self.index,state,var))
        return 1

    def solve(self,supply,passes=2):
        for h,line in enumerate(supply):
            for t,value in enumerate(line[1:],1):
                self.links[h][t].SetUb(value)
        scale = 10**6
        answer = None
        for step in range(passes):
            status = self.solver.Solve()
            if status!=pywraplp.Solver.OPTIMAL:
                self.numerical_retries += 1
                self.solver.SetSolverSpecificParametersAsString('use_dual_simplex: false\nuse_preprocessing: false')
                status = self.solver.Solve()
                self.solver.SetSolverSpecificParametersAsString('use_dual_simplex: true\nuse_preprocessing: true')
            assert status==pywraplp.Solver.OPTIMAL,(self.index,status,len(self.columns))
            multipliers = [[0]+[min(0,round(c.dual_value()*self.cost_scale*scale)) for c in line[1:]] for line in self.links]
            common = [np.cumsum(np.array(line,dtype=np.int64)) for line in multipliers]
            priced = self.model.price(self.parity,common,scale,batch=12,diverse=True)
            mass_dual = self.mass.dual_value()*self.cost_scale
            constant = priced[0][0]
            bound = (constant+sum(v*p for row,prob in zip(multipliers,supply) for v,p in zip(row,prob)))/scale
            answer = dict(constant=constant,multipliers=multipliers,denominator=scale,
                          primal=self.objective.Value()*self.cost_scale,lower=bound,
                          slack_cost=sum(self.bounds[h]*var.solution_value() for h,t,var in self.slacks))
            if answer['primal']-bound<1e-5:
                break
            # Do not modify the LP after the last solve: its solution is
            # needed immediately to reconstruct a feasible primal witness.
            if step+1<passes:
                added = sum(self.add(state) for value,state in priced if value/scale<mass_dual-1e-5)
                if not added:
                    break
        self.pending = [state for value,state in priced if value/scale<mass_dual-1e-5]
        return answer


def group_prices(source,model,shapes,cut_rows,z_lower,seconds,pointwise=False,pool=None,catalogue=None):
    """Snapshot and round a master dual, then search all group assignments.

    CP upper bounds here are exploratory until independently proof-checked.
    All returned assignments are usable columns regardless of proof status.
    """
    scale,event_scale = 10**6,1000
    by_point = [[] for x in source['core_residual']]
    for x,cut,constraint in cut_rows:
        by_point[x].append((cut,max(0,constraint.dual_value())))
    weights = [[[0]*(len(shape)+1) for x in source['core_residual']] for shape in shapes]
    constant = 0
    point_weights = []
    for x,rows in enumerate(by_point):
        total = max(1,sum(value for cut,value in rows))
        masses = [math.floor(value/total*scale) for cut,value in rows]
        assert sum(masses)<=scale
        point_weights.append(sum(masses) if pointwise else scale)
        constant += (point_weights[-1]-sum(masses))*scale*z_lower[x]
        for (cut,value),w in zip(rows,masses):
            assert cut['denominator']==scale
            constant += w*cut['constant']
            for h,line in enumerate(cut['multipliers']):
                for t,v in enumerate(line):weights[h][x][t] -= w*v
    denominator = scale*scale
    quantum = denominator//event_scale
    priced = []
    event_multipliers = [[[0]*len(line) for line in group] for group in weights]
    for h,shape in enumerate(shapes):
        events = [(x,t) for x in range(len(by_point)) for t in range(1,len(shape)+1) if weights[h][x][t]]
        integers = [(weights[h][x][t]+quantum-1)//quantum for x,t in events]
        for (x,t),w in zip(events,integers):event_multipliers[h][x][t] = -w
        known = catalogue.price(shape,events,integers) if catalogue is not None else []
        row = maximize(source['core_residual'],shape,events,integers,seconds,collect=True,
                       hint=known[0]['residues'] if known else None)
        row.setdefault('candidate_assignments',[]).extend(known)
        row['catalogue_best_score'] = known[0]['objective'] if known else None
        if catalogue is not None:
            for candidate in row['candidate_assignments']:catalogue.add(shape,candidate['residues'])
        capacity = event_majorant(source['core_residual'],shape,events,integers)
        row['capacity_certificate'] = capacity
        row['cp_upper_bound'] = row.get('upper_bound')
        if 'upper_bound' not in row or capacity['bound']<=row['upper_bound']:
            row['upper_bound'] = capacity['bound']
            row['bound_source'] = 'integer_capacity'
        else:row['bound_source'] = 'cp_sat_needs_independent_proof'
        row['events'] = [dict(point=x,threshold=t,weight=w) for (x,t),w in zip(events,integers)]
        priced.append(row)
    master_numerator = constant-sum(row.get('upper_bound',0)*quantum for row in priced)
    complete = all('upper_bound' in row for row in priced)
    # Reprice the aggregated dual on the entire point-state domain. This is
    # stronger than the cut combination and makes its certificate independent
    # of the numerical master solution and the finite configuration pool.
    point_prices = []
    for x,v in enumerate(source['core_residual']):
        common = [scale*np.cumsum(np.array(group[x],dtype=np.int64)) for group in event_multipliers]
        point_prices.append(model.price(v%2,common,event_scale*point_weights[x])[0][0])
    numerator = sum(point_prices)-scale*sum(row.get('upper_bound',0) for row in priced)
    capacity_numerator = sum(point_prices)-scale*sum(row['capacity_certificate']['bound'] for row in priced)
    pool_upper = None if pool is None else [max(sum(e['weight'] for e in row['events'] if counts[e['point']]>=e['threshold'])
                         for residues,counts,var in configs) for row,configs in zip(priced,pool)]
    return dict(numerator=numerator if complete else None,denominator=event_scale*scale,
                lower=numerator/(event_scale*scale) if complete else None,
                master_lower=master_numerator/denominator if complete else None,
                point_multipliers=point_prices,event_multipliers=event_multipliers,
                point_weights=point_weights,point_weight_denominator=scale,
                capacity_numerator=capacity_numerator,capacity_lower=capacity_numerator/(event_scale*scale),
                pool_upper_bounds=pool_upper,pool_sizes=[len(rows) for rows in pool] if pool is not None else None,
                pool_numerator=sum(point_prices)-scale*sum(pool_upper) if pool_upper is not None else None,
                group_pricing=priced,event_scale=event_scale,
                scope='rounded valid master dual; CP-SAT group upper bounds need independent proof checks')


def solve(source,pool,integer=None,seed=None,seconds=120,passes=2,warm_start=None,bounded_duals=True,
          pricing_seconds=0,pricing_every=5,pointwise=False,level_fraction=0,compact_states=0,catalogue=None):
    start = time.monotonic()
    model = LocalStates(source['high_primes'],source['high_pattern'])
    R,n,G = source['core_residual'],model.n,len(model.groups)
    groups,shapes = structure(source)
    if catalogue is not None:assert catalogue.points==R
    solver = pywraplp.Solver.CreateSolver('GLOP')
    solver.SetSolverSpecificParametersAsString('use_dual_simplex: true')
    objective = solver.Objective(); objective.SetMinimization()
    mass = [solver.Constraint(1,1) for h in groups]
    configs,seen = [[] for h in groups],[set() for h in groups]
    cut_rows = []
    assert 0<=level_fraction<1
    z,z_lower,objective_vars = [],[],[]
    for x in R:
        D = [q-1-(b==x%2) for q,b in zip(source['high_primes'],source['high_pattern'])]
        lower = -16*sum(math.prod(D[i] for i in range(n) if i not in S) for S in groups[n:])
        z_lower.append(lower)
        var = solver.NumVar(lower,solver.infinity(),f'cost{x}')
        z.append(var)
        if pointwise:
            positive = solver.NumVar(0,solver.infinity(),f'positive{x}')
            solver.Add(positive>=var)
            objective.SetCoefficient(positive,1)
            objective_vars.append(positive)
        else:
            objective.SetCoefficient(var,1)
            objective_vars.append(var)
    level_constraint = solver.Constraint(-solver.infinity(),solver.infinity())
    for var in objective_vars:level_constraint.SetCoefficient(var,1)
    overlaps,overlap_seen = [],set()
    def coefficient(cut,h,count):
        return sum(cut['multipliers'][h][1:count+1])/cut['denominator']
    def add_config(h,residues):
        if catalogue is not None:catalogue.add(shapes[h],residues)
        counts = tuple(sum(x%d==a for a,d in zip(residues,shapes[h])) for x in R)
        if counts in seen[h]:
            return 0
        seen[h].add(counts)
        var = solver.NumVar(0,solver.infinity(),f'g{h}_{len(configs[h])}')
        mass[h].SetCoefficient(var,1)
        for x,cut,c in cut_rows:
            c.SetCoefficient(var,-coefficient(cut,h,counts[x]))
        configs[h].append((tuple(residues),counts,var))
        return 1
    A = assignments_from_integer(source,integer) if integer else None
    warm_best = (warm_start or {}).get('coherent_primal_witness')
    warm_groups = {r['axis']:r['assignments'] for r in warm_best['groups']} if warm_best else {}
    for h,rows in enumerate(pool['group_columns']):
        if A:
            add_config(h,A[h])
        for row in rows:
            add_config(h,row)
        for row in warm_groups.get(h,[]):add_config(h,row['residues'])
        assert configs[h]
    seed_rows = list((seed or {}).get('support',[]))
    if warm_start:
        seed_rows += warm_start.get('point_columns',[])
    if warm_best:seed_rows += warm_best['support']
    if A:
        for x,v in enumerate(R):
            state = [counts[x] for rows in configs for residues,counts,var in rows[:1]]
            for h in range(n):
                state[h] = min(state[h],source['high_primes'][h]-1-(source['high_pattern'][h]==v%2))
            # With deterministic supplies F is affine in each plane count,
            # and decreasing in each tail count. All box vertices suffice.
            for plane in itertools.product(*[(0,state[h]) for h in range(n)]):
                seed_rows.append(dict(point_index=x,counts=list(plane)+state[n:]))
    oracles = [PointOracle(source,model,x,seed_rows,bounded_duals) for x in range(len(R))]
    def add_cut(x,cut):
        c = solver.Constraint(cut['constant']/cut['denominator'],solver.infinity())
        c.SetCoefficient(z[x],1)
        for h,rows in enumerate(configs):
            for residues,counts,var in rows:
                c.SetCoefficient(var,-coefficient(cut,h,counts[x]))
        cut_rows.append((x,cut,c))
    if warm_start:
        for row in warm_start.get('point_cuts',[]):
            add_cut(row['point_index'],row['cut'])
    # Start from a known coherent choice, independently of the master LP.
    probabilities = [[1.0]+[0.0]*(len(rows)-1) for rows in configs]
    best,best_value,history = None,math.inf,[]
    if warm_best:
        report = audit(dict(input=source,result=dict(coherent_primal_witness=warm_best)))
        best = warm_best
        key = 'positive_point_means_numerator' if pointwise else 'objective_numerator'
        best_value = report[key]/report['denominator']
        for h,rows in enumerate(configs):
            values = {}
            for assignment in warm_groups[h]:
                counts = tuple(sum(v%d==a for d,a in zip(shapes[h],assignment['residues'])) for v in R)
                values[counts] = values.get(counts,0)+assignment['weight']/warm_best['denominator']
            probabilities[h] = [values.get(counts,0) for residues,counts,var in rows]
    pricing_history = []
    for iteration in itertools.count(1):
        supplies = [[[0.0]*(len(shape)+1) for h,shape in enumerate(shapes)] for x in R]
        for h,rows in enumerate(configs):
            for (residues,counts,var),p in zip(rows,probabilities[h]):
                if p:
                    for x,v in enumerate(counts):
                        for t in range(1,v+1):supplies[x][h][t] += p
        answers,compact_solved,current_point_columns = [],0,[]
        for x,(oracle,supply) in enumerate(zip(oracles,supplies)):
            compact = compact_point(source,R[x],supply,compact_states) if compact_states else None
            if compact is not None:
                compact_solved += 1
                for state,p in compact['support']:oracle.add(state)
                multipliers = compact['event_multipliers']
                scale = compact['denominator']
                common = [np.cumsum(np.array(line,dtype=np.int64)) for line in multipliers]
                constant = model.price(R[x]%2,common,scale)[0][0]
                lower = (constant+sum(v*p for row,prob in zip(multipliers,supply) for v,p in zip(row,prob)))/scale
                answers.append(dict(constant=constant,multipliers=multipliers,denominator=scale,
                                    primal=compact['objective'],lower=lower,slack_cost=0))
                current_point_columns.extend((x,tuple(state),FixedValue(p)) for state,p in compact['support'])
                oracle.pending = []
            else:
                answers.append(oracle.solve(supply,passes))
                current_point_columns.extend(oracle.columns)
        upper = sum(max(0,row['primal']) if pointwise else row['primal'] for row in answers)
        # rational_witness expects objects with solution_value. The current
        # group probabilities precede the next master solve and must be kept.
        current_configs = [[(residues,counts,FixedValue(p)) for (residues,counts,var),p in zip(rows,probs)]
                           for rows,probs in zip(configs,probabilities)]
        candidate = rational_witness(source,model,current_point_columns,current_configs,shapes)
        key = 'positive_point_means_numerator' if pointwise else 'objective_numerator'
        value = candidate[key]/candidate['denominator']
        if value<best_value:
            best,best_value = candidate,value
            audit(dict(input=source,result=dict(coherent_primal_witness=best)))
        for x,row in enumerate(answers):add_cut(x,row)
        solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
        if solver.Solve()!=pywraplp.Solver.OPTIMAL:break
        lower = objective.Value()
        probabilities = [[max(0,var.solution_value()) for residues,counts,var in rows] for rows in configs]
        probabilities = [[p/sum(row) for p in row] for row in probabilities]
        converged = best_value-lower<1e-3
        new_groups = 0
        if (pricing_seconds and (not pointwise or lower>1e-6)
                and (iteration%pricing_every==0 or converged) and time.monotonic()-start<seconds):
            priced = group_prices(source,model,shapes,cut_rows,z_lower,pricing_seconds,pointwise,configs,catalogue)
            pricing_history.append(dict(iteration=iteration,**priced))
            for h,row in enumerate(priced['group_pricing']):
                candidates = row.get('candidate_assignments',[])+([row] if 'residues' in row else [])
                for candidate in candidates:new_groups += add_config(h,candidate['residues'])
                probabilities[h] += [0.0]*(len(configs[h])-len(probabilities[h]))
        stabilized = False
        if (level_fraction and not converged and not new_groups and best is not None
                and time.monotonic()-start<seconds):
            # Select a nearby point of an intermediate objective level set.
            # The earlier primary optimum remains the only pool lower bound.
            # Sum min(p_j,center_j) maximizes overlap (minimizes L1 distance).
            for h,rows in enumerate(configs):
                for residues,counts,var in rows:
                    if (h,counts) not in overlap_seen:
                        overlap_seen.add((h,counts))
                        overlap = solver.NumVar(0,1,f'overlap{h}_{len(overlaps)}')
                        solver.Add(overlap<=var)
                        overlaps.append((h,counts,overlap))
            centers = [dict() for h in groups]
            for group in best['groups']:
                h = group['axis']
                for assignment in group['assignments']:
                    counts = tuple(sum(v%d==a for d,a in zip(shapes[h],assignment['residues'])) for v in R)
                    centers[h][counts] = centers[h].get(counts,0)+assignment['weight']/best['denominator']
            for h,counts,var in overlaps:
                var.SetUb(centers[h].get(counts,0))
                objective.SetCoefficient(var,-1)
            for var in objective_vars:objective.SetCoefficient(var,0)
            level_constraint.SetUb(lower+level_fraction*(best_value-lower))
            solver.SetTimeLimit(max(1,int(1000*(seconds-time.monotonic()+start))))
            if solver.Solve()==pywraplp.Solver.OPTIMAL:
                probabilities = [[max(0,var.solution_value()) for residues,counts,var in rows] for rows in configs]
                probabilities = [[p/sum(row) for p in row] for row in probabilities]
                stabilized = True
            # Restoring the primary problem invalidates the solver solution;
            # all probabilities have been snapshotted before this mutation.
            level_constraint.SetUb(solver.infinity())
            for h,counts,var in overlaps:objective.SetCoefficient(var,0)
            for var in objective_vars:objective.SetCoefficient(var,1)
        row = dict(iteration=iteration,upper=best_value,pool_lower=lower if not new_groups else None,
                   pool_lower_before_group_additions=lower,
                   level_stabilized=stabilized,
                   compact_points_solved=compact_solved,
                   current_repaired_value=value,penalized_primal=upper,
                   oracle_gap=sum(a['primal']-a['lower'] for a in answers),
                   slack_cost=sum(a['slack_cost'] for a in answers),
                   numerical_retries=sum(o.numerical_retries for o in oracles),
                   added_group_columns=new_groups,
                   exploratory_global_lower=pricing_history[-1]['lower'] if pricing_history else None,
                   point_columns=sum(len(o.columns) for o in oracles),group_columns=sum(map(len,configs)),
                   seconds=time.monotonic()-start)
        history.append(row);print(json.dumps(row),flush=True)
        if (best_value==0 if pointwise else best_value<0) or (converged and not new_groups) or time.monotonic()-start>=seconds:break
        for oracle in oracles:
            for state in oracle.pending:oracle.add(state)
    result = dict(coherent_primal_witness=best,history=history,seconds=time.monotonic()-start,
                  group_columns=[[list(residues) for residues,counts,var in rows] for rows in configs],
                  point_columns=[dict(point_index=x,counts=list(state)) for oracle in oracles for x,state,var in oracle.columns],
                  point_cuts=[dict(point_index=x,cut=cut) for x,cut,c in cut_rows],
                  bounded_duals=bounded_duals,
                  pricing_history=pricing_history,
                  objective_kind='positive_point_means' if pointwise else 'raw_sum',
                  objective_value=best_value,
                  level_fraction=level_fraction,
                  compact_states=compact_states,
                  catalogue=catalogue.export() if catalogue is not None else None,
                  scope='explicit finite configuration pool; numerical lower bounds do not exclude omitted assignments')
    result['audit'] = audit(dict(input=source,result=result))
    return result


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True)
    ap.add_argument('--pool',required=True)
    ap.add_argument('--integer',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--time',type=float,default=120)
    ap.add_argument('--passes',type=int,default=2)
    ap.add_argument('--warm-start')
    ap.add_argument('--unbounded-duals',action='store_true')
    ap.add_argument('--pricing-seconds',type=float,default=0)
    ap.add_argument('--pricing-every',type=int,default=5)
    ap.add_argument('--pointwise',action='store_true')
    ap.add_argument('--level-fraction',type=float,default=0)
    ap.add_argument('--compact-states',type=int,default=0)
    ap.add_argument('--catalogue')
    args = ap.parse_args()
    data,pool,integer = [json.loads(Path(p).read_text()) for p in (args.source,args.pool,args.integer)]
    for key in ('core_residual','high_primes','high_pattern'):
        assert data['input'][key]==pool['input'][key]==integer['input'][key]
    warm = json.loads(Path(args.warm_start).read_text()) if args.warm_start else None
    if warm:assert warm['input']==data['input']
    catalogue = Catalogue.load(json.loads(Path(args.catalogue).read_text())) if args.catalogue else None
    result = solve(data['input'],pool['result'],integer['result'],data['result'].get('point_primal_witness'),args.time,args.passes,
                   warm_start=warm['result'] if warm else None,bounded_duals=not args.unbounded_duals,
                   pricing_seconds=args.pricing_seconds,pricing_every=args.pricing_every,pointwise=args.pointwise,
                   level_fraction=args.level_fraction,compact_states=args.compact_states,catalogue=catalogue)
    Path(args.output).write_text(json.dumps(dict(input=data['input'],result=result),separators=(',',':'))+'\n')
    print(json.dumps(dict(seconds=result['seconds'],iterations=len(result['history']),audit=result['audit'])))
