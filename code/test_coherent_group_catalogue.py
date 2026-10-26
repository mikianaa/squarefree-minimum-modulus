import itertools
import unittest
from contextlib import redirect_stdout
import io
from unittest.mock import patch
from coherent_group_catalogue import Catalogue
from group_event_coherence import maximize


class CatalogueTests(unittest.TestCase):
    def test_scores_match_exhaustive_assignments_and_roundtrip(self):
        points,shape = [0,1,4,7],[2,3]
        catalogue = Catalogue(points)
        assignments = list(itertools.product(range(2),range(3)))
        for residues in assignments:catalogue.add(shape,residues)
        events,weights = [(0,1),(1,2),(2,1),(3,2)],[2,13,5,19]
        scores = [sum(w for (x,t),w in zip(events,weights)
                      if sum(points[x]%d==a for d,a in zip(shape,residues))>=t)
                  for residues in assignments]
        result = catalogue.price(shape,events,weights,top=100)
        self.assertEqual([r['objective'] for r in result],sorted(scores,reverse=True))
        restored = Catalogue.load(dict(input=dict(core_residual=points),result=catalogue.export()))
        self.assertEqual(restored.price(shape,events,weights,top=100),result)
        self.assertEqual(catalogue.price([2,5],events,weights),[])

    def test_reuse_witness_and_refresh_cached_scores(self):
        source = dict(core_residual=[1,2],high_primes=[11])
        shape = [d for d in range(3,211) if 210%d==0]
        catalogue = Catalogue(source['core_residual'])
        witness = dict(coherent_primal_witness=dict(groups=[dict(axis=0,assignments=[dict(residues=[1]*14)])]))
        self.assertEqual(catalogue.add_result(source,witness),1)
        self.assertEqual(catalogue.price(shape,[(1,14)],[9])[0]['objective'],0)
        self.assertEqual(catalogue.add_result(source,dict(group_columns=[[[2]*14]])),1)
        self.assertEqual(catalogue.price(shape,[(1,14)],[9])[0]['objective'],9)
        self.assertEqual(catalogue.add(shape,[2]*14),0)

    def test_hint_survives_timeout_and_only_returns_valid_candidates(self):
        points,shape,events,weights = [1,2],[3,5],[(0,1),(1,2)],[4,9]
        # Both hinted residues miss the point set; normalization must gain hits.
        result = maximize(points,shape,events,weights,seconds=0,collect=True,hint=[0,0])
        self.assertNotIn('upper_bound',result)
        self.assertNotIn('objective',result)
        self.assertEqual(len(result['candidate_assignments']),1)
        candidate = result['candidate_assignments'][0]
        self.assertEqual(candidate['counts'],[2,0])
        self.assertEqual(candidate['objective'],4)
        solved = maximize(points,shape,events,weights,seconds=1,hint=[0,0])
        self.assertEqual(solved['objective'],9)
        self.assertGreaterEqual(solved['upper_bound'],9)

    def test_master_adds_catalogue_assignment_when_solver_times_out(self):
        from decomposed_coherent_groups import solve
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        shape = [d for d in range(3,211) if 210%d==0]
        catalogue = Catalogue(source['core_residual'])
        balanced = [1]*7+[2]*7
        catalogue.add(shape,balanced)
        with patch('decomposed_coherent_groups.maximize',side_effect=lambda *a,**kw:dict(status='UNKNOWN')):
            with redirect_stdout(io.StringIO()):
                result = solve(source,dict(group_columns=[[[1]*14,[2]*14]]),seconds=5,
                               pricing_seconds=.1,pricing_every=1,catalogue=catalogue)
        self.assertIn(balanced,result['group_columns'][0])
        self.assertAlmostEqual(result['objective_value'],5,places=6)
        self.assertTrue(result['audit']['verified'])
        self.assertTrue(all(row['group_pricing'][0]['cp_upper_bound'] is None for row in result['pricing_history']))


if __name__=='__main__':unittest.main()
