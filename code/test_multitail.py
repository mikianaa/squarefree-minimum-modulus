import copy
import itertools
import json
from pathlib import Path
import unittest
import numpy as np
from multitail_lp import LocalStates
from audit_multitail import audit,local_minima
from tail_primal_witness import exact_witness
from rank_capacity import catalogue,violations
from audit_fibre_concentration import gray_capacities
from audit_multitail_identity import audit as audit_identity
from point_pair_capacity import catalogue as pair_catalogue,violations as pair_violations
from collections import Counter
from core_subset_thresholds import separate
from separate_subset_events import separate as separate_any_subset
from joint_triple_events import feasible as triple_feasible


class MultitailTests(unittest.TestCase):
    def test_joint_triple_feasibility_matches_resource_enumeration(self):
        shapes = [2,3,5,6]
        for points in ([1,2,4],[1,7,11],[0,6,12]):
            assignments = itertools.product(*({x%d for x in points} for d in shapes))
            hits = [tuple(sum(x%d==a for d,a in zip(shapes,choices)) for x in points) for choices in assignments]
            for thresholds in itertools.product(range(len(shapes)+1),repeat=3):
                brute = any(all(v>=t for v,t in zip(row,thresholds)) for row in hits)
                self.assertEqual(triple_feasible(points,shapes,thresholds),brute)

    def test_triple_audit_rejects_pair_feasible_witness(self):
        data = json.loads(Path('data/nine-prime-fixed-core-pair-rank-lp.json').read_text())
        data['result']['triple_cuts_enabled'] = True
        with self.assertRaises(AssertionError):
            audit_identity(data)

    def test_subset_separator_agrees_with_all_four_point_checks(self):
        data = json.loads(Path('data/nine-prime-fixed-core-triple-rank-lp.json').read_text())
        exhaustive = separate(data,4)
        general = separate_any_subset(data)
        def keys(rows):
            return {(r['axis'],tuple(r['points']),tuple(r['thresholds']),r['violation_numerator']) for r in rows}
        self.assertEqual(keys(exhaustive['violations']),keys(general['violations']))
        self.assertEqual(len(exhaustive['violations']),3)

    def test_subset_separator_two_point_result_matches_bruteforce(self):
        shapes = [d for d in range(3,211) if 210%d==0]
        R = [13,19]
        cap = len(shapes)+sum(R[0]%d==R[1]%d for d in shapes)
        for a,b in ((14,14),(1,1),(10,10),(8,9)):
            data = dict(input=dict(core_residual=R,high_primes=[11]),result=dict(point_primal_witness=dict(
                denominator=10,support=[dict(point_index=x,counts=[v],weight=10) for x,v in enumerate((a,b))])))
            result = separate_any_subset(data)
            self.assertEqual(bool(result['violations']),a+b>cap)
            self.assertTrue(all(r['status'] in ('OPTIMAL','INFEASIBLE') for r in result['axis_reports']))

    def test_stored_subset_audit_rejects_invalid_capacity(self):
        data = json.loads(Path('data/nine-prime-fixed-core-subset-rank-lp.json').read_text())
        self.assertTrue(audit_identity(data)['verified'])
        data['result']['subset_constraints'][0][2] = [1,1,1]
        with self.assertRaises(AssertionError):
            audit_identity(data)

    def test_three_point_cut_adds_information_beyond_pair_capacities(self):
        data = json.loads(Path('data/nine-prime-fixed-core-pair-rank-lp.json').read_text())
        report = separate(data)
        self.assertEqual(len(report['violations']),5)
        row = report['violations'][0]
        R = [data['input']['core_residual'][i] for i in row['points']]
        shapes = [d for d in range(1,211) if 210%d==0]
        cap = sum(max(Counter(x%d for x in R).values()) for d in shapes)
        self.assertEqual(cap,21)
        self.assertEqual(sum(row['thresholds']),22)
        for i,j in itertools.combinations(range(3),2):
            pair_cap = len(shapes)+sum(R[i]%d==R[j]%d for d in shapes)
            self.assertLessEqual(row['thresholds'][i]+row['thresholds'][j],pair_cap)

    def test_pair_event_cuts_reject_previously_feasible_rank_witness(self):
        data = json.loads(Path('data/nine-prime-fixed-core-identity-rank-lp.json').read_text())
        R = data['input']['core_residual']
        witness = data['result']['point_primal_witness']
        self.assertTrue(pair_violations(witness,R,5,pair_catalogue(R,5)))
        data['result']['pair_cuts_enabled'] = True
        with self.assertRaises(AssertionError):
            audit_identity(data)

    def test_pair_event_caps_for_all_two_point_resource_choices(self):
        R = [13,19]
        shapes = [d for d in range(3,211) if 210%d==0]
        options = [[(int(R[0]%d==a),int(R[1]%d==a)) for a in {R[0]%d,R[1]%d}] for d in shapes]
        rows = pair_catalogue(R,1)
        for assignment in itertools.product(*options):
            a,b = map(sum,zip(*assignment))
            witness = dict(denominator=1,support=[dict(point_index=0,counts=[a],weight=1),dict(point_index=1,counts=[b],weight=1)])
            self.assertFalse(pair_violations(witness,R,1,rows))

    def test_rank_inequalities_for_every_small_projection_assignment(self):
        R,moduli = [1,2,4,5],[2,3,6]
        H = gray_capacities(R,moduli)
        rows = catalogue(H,1)
        for choices in itertools.product(*(range(d) for d in moduli)):
            witness = dict(denominator=1,support=[dict(parity=x%2,counts=[sum(x%d==a for d,a in zip(moduli,choices))],weight=1) for x in R])
            self.assertFalse(violations(witness,1,rows))

    def test_rank_cuts_reject_the_audited_identity_relaxation_witness(self):
        data = json.loads(Path('data/nine-prime-fixed-core-identity-lp.json').read_text())
        source,witness = data['input'],data['result']['point_primal_witness']
        converted = dict(denominator=witness['denominator'],support=[dict(parity=source['core_residual'][r['point_index']]%2,counts=r['counts'],weight=r['weight']) for r in witness['support']])
        self.assertTrue(violations(converted,5,catalogue(source['H_upper'],5)))

    def test_identity_primal_rejects_modified_projection_mass(self):
        data = json.loads(Path('data/nine-prime-fixed-core-identity-lp.json').read_text())
        self.assertTrue(audit_identity(data)['verified'])
        data['result']['point_primal_witness']['projections'][0]['assignments'][0]['weight'] += 1
        with self.assertRaises(AssertionError):
            audit_identity(data)

    def test_pricing_against_all_two_prime_local_states(self):
        model = LocalStates([3,5],[0,1])
        common = [np.array([-max(0,v-i-1)*(i+2) for v in range(15 if i<2 else 17)],dtype=np.int64) for i in range(3)]
        for e in (0,1):
            D = [q-1-(b==e) for q,b in zip(model.high,model.pattern)]
            values = []
            for state in itertools.product(*(range(d+1) for d in D),range(17)):
                values.append(7*model.cost(e,state)-sum(int(c[v]) for c,v in zip(common,state)))
            priced = model.price(e,common,7,batch=5)
            self.assertEqual(priced[0][0],min(values))
            for value,state in priced:
                self.assertEqual(value,7*model.cost(e,state)-sum(int(c[v]) for c,v in zip(common,state)))
            source = dict(high_primes=model.high,high_pattern=model.pattern)
            # The independent checker takes threshold multipliers; test its
            # minima separately with a single known threshold in each group.
            certificate = dict(denominator=7,tail_multipliers=[dict(axis=i,thresholds=[i+1,i+1],multiplier=-(i+2)) for i in range(3)])
            independent,_ = local_minima(source,certificate)
            self.assertEqual(independent[e],min(values))

    def test_four_prime_certificate_reproduced_by_general_audit(self):
        data = json.loads(Path('data/critical-010-paired-bound.json').read_text())
        data['input']['high_primes'] = [11,13,17,19]
        self.assertEqual(audit(data,primal=False)['integer_lower'],1324)

    def test_general_rational_repair_includes_correct_zero_state(self):
        class Variable:
            def solution_value(self):
                return 0.99
        H = [[0],[14]]
        rows = []
        for h,limit in enumerate((14,14,16)):
            for a,b in itertools.product(range(limit+1),repeat=2):
                bound = max(v+(x+y+max(x,y) if h==2 else 0)-a*x-b*y for x,line in enumerate(H) for y,v in enumerate(line))
                rows.append((h,a,b,bound,None))
        witness = exact_witness([(0,(9,11,16),Variable())],[0,1],[1,0],rows,[11,13],[1,1])
        self.assertTrue(all(len(r['counts'])==3 for r in witness['support']))
        self.assertLess(witness['objective_numerator'],0)
        data = dict(input=dict(high_primes=[11,13],high_pattern=[1,1],counts=[1,0],H_upper=H),
                    result=dict(primal_relaxation_witness=witness))
        self.assertTrue(audit(data)['verified'])
        bad = copy.deepcopy(data)
        bad['result']['primal_relaxation_witness']['support'][0]['weight'] += 1
        with self.assertRaises(AssertionError):
            audit(bad)


if __name__=='__main__':
    unittest.main()
