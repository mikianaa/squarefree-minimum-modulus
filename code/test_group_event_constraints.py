import copy
import itertools
import json
from pathlib import Path
import unittest
from pysat.solvers import Solver
from audit_group_event_cut import build
from audit_multitail_identity import audit,audit_group_cut
from group_event_coherence import maximize,simplify_candidate
from multitail_identity import solve
from joint_subset_events import attainable,partitions,feasible_integer
from sat_uncovered_audit import check_rup


class GroupEventTests(unittest.TestCase):
    def test_resume_accepts_valid_nonboundary_pair_cut_and_rejects_invalid_one(self):
        data = json.loads(Path('data/nine-prime-fixed-core-group-events-lp.json').read_text())
        R,n = data['input']['core_residual'],len(data['input']['high_primes'])
        rows = data['result']['pair_constraints']
        def cap(row):
            h,x,y,t,u = row
            shapes = [d for d in range(1 if h>=n else 3,211) if 210%d==0]
            return len(shapes)+sum(R[x]%d==R[y]%d for d in shapes)
        self.assertTrue(any(r[3]+r[4]>cap(r)+1 for r in rows))
        result = solve(data['input'],seconds=0,pair_cuts=True,warm_start=data['result'])
        self.assertEqual(result['pair_constraints'],rows)
        bad = copy.deepcopy(data['result'])
        bad['pair_constraints'][0][3:] = [1,1]
        with self.assertRaises(AssertionError):
            solve(data['input'],seconds=0,pair_cuts=True,warm_start=bad)

    def test_small_integer_separation_retains_exact_violation(self):
        raw = json.loads(Path('data/nine-prime-group-coherence-separation.json').read_text())
        data = json.loads(Path(raw['source']).read_text())
        candidate = next(r['cut'] for r in raw['results'] if r['axis']==0)
        cut = simplify_candidate(data,candidate,scales=(16,))
        self.assertIsNotNone(cut)
        self.assertEqual(cut['bound'],6)
        used = sum(r['weight']*e['weight'] for r in data['result']['point_primal_witness']['support']
                   for e in cut['events'] if r['point_index']==e['point'] and r['counts'][0]>=e['threshold'])
        self.assertEqual(used-cut['bound']*cut['denominator'],cut['violation_numerator'])
        self.assertGreater(cut['violation_numerator'],0)

    def test_four_point_downward_closure_matches_direct_assignments(self):
        points,shapes = [1,2,4,7],[2,3,5]
        vectors = [tuple(sum(x%d==a for d,a in zip(shapes,choices)) for x in points)
                   for choices in itertools.product(*({x%d for x in points} for d in shapes))]
        possible = attainable(partitions(points,shapes),4)
        for thresholds in itertools.product(range(4),repeat=4):
            brute = any(all(v>=t for v,t in zip(row,thresholds)) for row in vectors)
            self.assertEqual(bool(possible[thresholds]),brute)
            self.assertEqual(feasible_integer(points,shapes,thresholds),brute)

    def test_weighted_sat_encoding_matches_all_resource_choices(self):
        R,shapes = [1,2,4],[2,3,6]
        events = [dict(point=0,threshold=2,weight=2),dict(point=0,threshold=3,weight=1),
                  dict(point=1,threshold=2,weight=3),dict(point=2,threshold=1,weight=1)]
        best = 0
        for choices in itertools.product(*({x%d for x in R} for d in shapes)):
            counts = [sum(x%d==a for d,a in zip(shapes,choices)) for x in R]
            best = max(best,sum(r['weight'] for r in events if counts[r['point']]>=r['threshold']))
        priced = maximize(R,shapes,[(r['point'],r['threshold']) for r in events],[r['weight'] for r in events],collect=True)
        self.assertTrue(priced['candidate_assignments'])
        for row in priced['candidate_assignments']:
            counts = [sum(x%d==a for d,a in zip(shapes,row['residues'])) for x in R]
            self.assertEqual(row['counts'],counts)
            self.assertEqual(row['objective'],sum(r['weight'] for r in events if counts[r['point']]>=r['threshold']))
        self.assertEqual(priced['objective'],best)
        self.assertEqual(priced['upper_bound'],best)
        for bound in range(sum(r['weight'] for r in events)):
            clauses,_ = build(R,shapes,events,bound)
            # This Windows Glucose build crashes on shutdown for some SAT
            # formulas when proof logging is enabled. Check both outcomes
            # with Minisat, and request proofs only for the UNSAT instances.
            with Solver(name='minisat22',bootstrap_with=clauses) as solver:
                answer = solver.solve()
                self.assertEqual(answer,best>bound)
            if not answer:
                with Solver(name='glucose4',bootstrap_with=clauses,with_proof=True) as solver:
                    self.assertFalse(solver.solve())
                    proof = solver.get_proof()
                self.assertTrue(check_rup(clauses,proof)['verified'])

    def test_new_group_cut_rejects_old_fractional_witness(self):
        data = json.loads(Path('data/nine-prime-fixed-core-subset-rank-lp.json').read_text())
        cuts = json.loads(Path('data/nine-prime-group-event-certified-cuts.json').read_text())['cuts']
        data['result']['group_event_constraints'] = [next(c for c in cuts if c['axis']==0)]
        with self.assertRaises(AssertionError):
            audit(data)

    def test_bound_tampering_changes_proof_encoding(self):
        data = json.loads(Path('data/nine-prime-fixed-core-subset-rank-lp.json').read_text())
        cut = json.loads(Path('data/nine-prime-group-event-certified-cuts.json').read_text())['cuts'][0]
        cut = copy.deepcopy(cut)
        cut['bound'] -= 1
        with self.assertRaises(AssertionError):
            audit_group_cut(data['input'],cut)


if __name__=='__main__':
    unittest.main()
