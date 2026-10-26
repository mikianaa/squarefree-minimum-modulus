import copy
import itertools
import json
from pathlib import Path
import unittest

from audit_event_capacity_bound import audit as audit_bound
from audit_multitail_identity import audit_group_cut
from weighted_event_capacity import capacities


class WeightedEventCapacityTests(unittest.TestCase):
    def test_capacity_is_exact_for_small_weighted_resource_assignments(self):
        points,weights,moduli = [1,2,4,7],[1,3,2,1],[2,3,5]
        best = max(sum(w*sum(x%d==a for d,a in zip(moduli,choices))
                       for x,w in zip(points,weights))
                   for choices in itertools.product(*(range(d) for d in moduli)))
        self.assertEqual(best,sum(capacities(points,weights,moduli)))

    def fixture(self):
        source = dict(core_residual=[67,89,97,139,164],high_primes=[11])
        weights = [2,2,1,1,2]
        moduli = [d for d in range(3,211) if 210%d==0]
        cut = dict(axis=0,bound=4,
                   events=[dict(point=x,threshold=t,weight=1) for x,t in enumerate([4,8,4,4,4])],
                   capacity_certificate=dict(point_weights=weights,
                       resource_capacities=capacities(source['core_residual'],weights,moduli)),
                   proof_audit=dict(kind='weighted_capacity',verified=False))
        return source,cut

    def test_five_point_cut_uses_arithmetic_not_verified_flag(self):
        source,cut = self.fixture()
        self.assertEqual(audit_group_cut(source,cut),0)
        self.assertEqual(sum(cut['capacity_certificate']['resource_capacities']),39)

    def test_weighted_cut_rejects_changed_threshold_bound_or_resource_capacity(self):
        source,cut = self.fixture()
        variants = []
        bad = copy.deepcopy(cut); bad['events'][1]['threshold'] = 7; variants.append(bad)
        bad = copy.deepcopy(cut); bad['bound'] = 3; variants.append(bad)
        bad = copy.deepcopy(cut); bad['capacity_certificate']['resource_capacities'][0] -= 1; variants.append(bad)
        bad = copy.deepcopy(cut); bad['events'][1]['point'] = 0; variants.append(bad)
        for bad in variants:
            with self.assertRaises(AssertionError):
                audit_group_cut(source,bad)

    def test_compact_bound_checks_every_event_profile_and_is_sharp(self):
        certificate = json.loads(Path('data/small-event-bound-capacity-certificate.json').read_text())
        result = audit_bound(certificate)
        self.assertEqual(result['threshold_profiles_checked'],3072)
        self.assertEqual(result['bound'],6)
        self.assertEqual(result['attained_score'],6)
        certificate['bound'] = 5
        with self.assertRaises(AssertionError):
            audit_bound(certificate)


if __name__=='__main__':
    unittest.main()
