import copy
import itertools
import json
from pathlib import Path
import unittest
import numpy as np
from paired_slice_capacity import paired_capacity
from audit_paired_slice_capacity import audit as audit_pairs
from audit_tail_primal import audit as audit_primal
from audit_abstract_tail_bound import audit as audit_dual
from audit_fibre_concentration import gray_capacities
from audit_joint_tail_envelopes import histogram_profile, from_profile
from joint_tail_profiles import joint_profile
from portable_tail_certificates import CoreCapacityUpper
from portable_tail_certificates import concentration_bounds
from audit_paired_tail_branch import direct_concentration


class PairedCapacityTests(unittest.TestCase):
    def test_independent_direct_concentration_transform(self):
        H = np.array([[0,5,8],[5,9,12],[8,12,15]],dtype=np.int64)
        self.assertTrue(np.array_equal(direct_concentration(H),np.stack([concentration_bounds(H,s) for s in (14,16)])))

    def test_pair_table_equals_exhaustive_subset_capacity(self):
        R = [1,2,4,5,7,8,11,14]
        for moduli in ([2,3],[3,5],[3,6]):
            certificate = paired_capacity(R,moduli)
            self.assertEqual(certificate['H_upper'],gray_capacities(R,moduli))
            self.assertTrue(audit_pairs(certificate)['verified'])

    def test_matching_upper_sandwiches_exact_capacity(self):
        for R in ([],[1],[2],[1,2,4,5,7,8,11,14]):
            moduli = [2,3,5,6]
            certificate = paired_capacity(R,moduli)
            H = np.array(certificate['H_upper'])
            self.assertTrue(np.all(H>=np.array(gray_capacities(R,moduli))))
            layers = [[n for n in R if n%2==e] for e in (0,1)]
            for a,b in itertools.product(range(len(layers[0])+1),range(len(layers[1])+1)):
                independent = sum(max(min(a,sum(n%d==r for n in layers[0]))+
                                      min(b,sum(n%d==r for n in layers[1])) for r in range(d)) for d in moduli)
                self.assertLessEqual(H[a,b],independent)
            self.assertTrue(audit_pairs(certificate)['verified'])

    def test_rejects_repeated_resource_in_matching(self):
        data = paired_capacity([1,2,3,4],[2,3,5,6])
        data['matching_witnesses'][0][0] = [[0,1],[0,3]]
        with self.assertRaises(AssertionError):
            audit_pairs(data)

    def test_rejects_understated_capacity(self):
        data = paired_capacity([1,2,3,4],[2,3,5,6])
        data['H_upper'][1][1] -= 1
        with self.assertRaises(AssertionError):
            audit_pairs(data)

    def test_profiles_agree_with_independent_ordinary_histograms(self):
        for pattern in itertools.product((0,1),repeat=2):
            model = CoreCapacityUpper((0,)+pattern)
            mask = sum(1<<j for j in range(0,len(model.numbers),4))
            R = [n for j,n in enumerate(model.numbers) if mask>>j&1]
            profile = histogram_profile(R)
            self.assertEqual(joint_profile(model,mask),profile)
            self.assertTrue(np.array_equal(model.upper(mask),from_profile(profile,model.compression)))

    def test_rational_primal_and_mass_tampering(self):
        data = json.loads(Path('data/joint-tail-010-primal-probe.json').read_text())
        self.assertTrue(audit_primal(data)['verified'])
        data['result']['primal_relaxation_witness']['support'][0]['weight'] += 1
        with self.assertRaises(AssertionError):
            audit_primal(data)

    def test_dual_rejects_mismatched_capacity(self):
        data = json.loads(Path('data/critical-010-paired-bound.json').read_text())
        capacity = json.loads(Path('data/critical-010-paired-capacity.json').read_text())
        self.assertEqual(audit_dual(data,capacity)['integer_lower'],1324)
        bad = copy.deepcopy(capacity)
        bad['H_upper'][1][1] += 1
        with self.assertRaises(AssertionError):
            audit_dual(data,bad)


if __name__=='__main__':
    unittest.main()
