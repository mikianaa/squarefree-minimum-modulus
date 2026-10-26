"""Validation of the geometry and new comparison calculations."""
import hashlib,json,unittest
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from core210_joint_probe import DIVS,price
from core210_continuous_hull import hull

ROOT=Path(__file__).resolve().parents[1]

class Core210Tests(unittest.TestCase):
    def test_fractional_and_capped_pricing_witnesses(self):
        d=json.loads((ROOT/'data/core210-hard2-dyadic-probe.json').read_text());p=np.array(d['support']);w=np.array(d['weights'])
        for t,cap in [(1.25,16),(2.5,16),(3,11)]:
            value,rows,residues=price(p,w,t,np.random.default_rng(710),random_starts=5,sweeps=3,ceiling=cap)
            scores=[]
            for row,a in zip(rows,residues):
                counts=np.array([sum(x%m==r for m,r in zip(DIVS,a)) for x in p])
                expected=np.maximum(np.minimum(counts,cap)-t,0)
                np.testing.assert_array_equal(row,expected);scores.append(expected@w)
            self.assertAlmostEqual(value,max(scores),places=11)

    def test_uniform_full_core_hull(self):
        atoms,weights=hull(dict(support=list(range(210)),weights=[1/210]*210,pool=[]))
        original_values=np.array([1,2,4,8,16]);original_weights=np.array([48,92,56,13,1])/210
        for t in np.arange(0,16.01,.25):
            self.assertAlmostEqual(weights@np.maximum(atoms-t,0),original_weights@np.maximum(original_values-t,0),places=10)

    def test_dyadic_basis_matches_full_core_comparison(self):
        b=json.loads((ROOT/'data/core210-hard2-dyadic-basis.json').read_text())
        g=json.loads((ROOT/b['schedule_file']).read_text())['G']
        f=np.array([1,2,4,8,16]);w=np.array([48,92,56,13,1])/210
        value=b['constant']+b['mean_coefficient']*(w@f-1)+b['top_coefficient']/210
        value+=sum(c*(w@np.maximum(f-float(F(t)),0)) for t,c in zip(b['thresholds'],b['coefficients']))
        self.assertAlmostEqual(value,sum(v*g[int(x)] for x,v in zip(f,w)),places=10)

    def test_catalogue_audit_and_hard_case(self):
        raw=(ROOT/'data/core210-maximum-prefixes.json').read_bytes();d=json.loads(raw)
        a=json.loads((ROOT/'data/audit-core210-maximum-prefixes.json').read_text())
        self.assertEqual(hashlib.sha256(raw).hexdigest(),a['catalogue_sha256']);self.assertTrue(a['passed'])
        self.assertEqual(a['normalized_assignments'],35490420);self.assertEqual(a['n_orbits'],288846)
        h=json.loads((ROOT/'data/core210-hard-prefix-classification.json').read_text())
        self.assertEqual(d['cases'][h['catalogue_index_zero_based']][:6],h['canonical_residues'])

    def test_audited_fixed_schedule_bound(self):
        a=json.loads((ROOT/'data/audit-core210-fixed-schedule.json').read_text())
        self.assertTrue(a['passed']);self.assertGreater(F(a['lower_bound_exact']),F(1001192376,10**9))
        self.assertEqual(a['point_checks'],104);self.assertEqual(a['positive_dual_terms'],123)

if __name__=='__main__':unittest.main()
