import copy
import json
from pathlib import Path
import unittest
from audit_coherent_groups import audit
from coherent_group_model import solve
from multitail_lp import LocalStates
from event_dual_seed import convert


class CoherentGroupTests(unittest.TestCase):
    def test_pointwise_audit_rejects_false_no_positive_fibre_claim(self):
        source = dict(core_residual=[1,2],high_primes=[3,5],high_pattern=[0,1])
        shapes = [[d for d in range(3 if h<2 else 1,211) if 210%d==0] for h in range(3)]
        witness = dict(denominator=1,objective_numerator=-13,positive_point_means_numerator=3,
                       groups=[dict(axis=h,assignments=[dict(residues=[1%d for d in shape],weight=1)]) for h,shape in enumerate(shapes)],
                       support=[dict(point_index=0,counts=[2,3,16],weight=1),
                                dict(point_index=1,counts=[0,0,1],weight=1)])
        data = dict(input=source,result=dict(coherent_primal_witness=witness))
        report = audit(data)
        self.assertTrue(report['negative_cost'])
        self.assertFalse(report['all_point_means_nonpositive'])
        self.assertEqual(report['point_cost_numerators'],[-16,3])
        witness['positive_point_means_numerator'] = 0
        with self.assertRaises(AssertionError):audit(data)

    def test_event_seed_preserves_every_old_penalty_value(self):
        data = json.loads(Path('data/nine-prime-group-refinement-round3-model.json').read_text())
        source,c = data['input'],data['result']['certificate']
        self.assertEqual(c['denominator'],10**6)
        arrays = convert(source,c)
        for h,group in enumerate(arrays):
            for x,line in enumerate(group):
                parity = source['core_residual'][x]%2
                for v in range(len(line)):
                    expected = c['link_multipliers'][x][h]*v
                    expected += sum(r['multiplier']*max(0,v-r['thresholds'][parity]) for r in c['tail_multipliers'] if r['axis']==h)
                    expected += sum(r['multiplier']*r['weights'][parity] for r in c.get('rank_multipliers',[]) if r['axis']==h and v>=r['thresholds'][parity])
                    for kind in ('pair_multipliers','subset_multipliers'):
                        expected += sum(r['multiplier'] for r in c.get(kind,[]) if r['axis']==h
                                        for point,t in zip(r['points'],r['thresholds']) if point==x and v>=t)
                    expected += sum(r['multiplier']*e['weight'] for r in c.get('group_event_multipliers',[]) if r['axis']==h
                                    for e in r['events'] if e['point']==x and v>=e['threshold'])
                    self.assertEqual(sum(line[:v+1]),expected)

    def test_stabilized_pricing_retains_small_exact_bound(self):
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        result = solve(source,seconds=10,target_lower=5,stabilization=.8)
        self.assertEqual(result['certificate']['integer_lower'],5)

    def integer_fixture(self):
        data = json.loads(Path('data/nine-prime-fixed-core-integer-projections.json').read_text())
        source,result = data['input'],data['result']
        model = LocalStates(source['high_primes'],source['high_pattern'])
        choices = {(r['axis'],r['modulus']):r['residue'] for r in result['choices']}
        shapes = [[d for d in range(1 if h>=model.n else 3,211) if 210%d==0] for h in range(len(model.groups))]
        groups = [dict(axis=h,assignments=[dict(residues=[choices[h,d] for d in shape],weight=1)]) for h,shape in enumerate(shapes)]
        support = []
        for x,v in enumerate(source['core_residual']):
            counts = [sum(v%d==choices[h,d] for d in shape) for h,shape in enumerate(shapes)]
            for h in range(model.n):
                counts[h] = min(counts[h],source['high_primes'][h]-1-(source['high_pattern'][h]==v%2))
            support.append(dict(point_index=x,counts=counts,weight=1))
        witness = dict(denominator=1,objective_numerator=result['integer_expression'],groups=groups,support=support)
        return dict(input=source,result=dict(coherent_primal_witness=witness))

    def test_known_integer_configuration_is_feasible_for_full_group_model(self):
        report = audit(self.integer_fixture())
        self.assertTrue(report['verified'])
        self.assertFalse(report['negative_cost'])
        self.assertEqual(report['objective_numerator'],678850)
        self.assertEqual(report['threshold_links_checked'],10692)

    def test_group_distribution_mass_tampering_is_rejected(self):
        data = self.integer_fixture()
        data['result']['coherent_primal_witness']['groups'][0]['assignments'][0]['weight'] += 1
        with self.assertRaises(AssertionError):
            audit(data)

    def test_point_event_not_supported_by_group_configuration_is_rejected(self):
        data = self.integer_fixture()
        witness = data['result']['coherent_primal_witness']
        witness['support'][0]['counts'][-1] = 16
        with self.assertRaises(AssertionError):
            audit(data)

    def test_one_group_pricing_recovers_exact_small_optimum(self):
        # Every one of 14 resources hits at most one of these points. Their
        # free-coordinate counts are 10 and 9, so exactly 5 remain at best.
        source = dict(core_residual=[1,2],high_primes=[11],high_pattern=[0])
        result = solve(source,seconds=10,target_lower=5)
        self.assertEqual(result['certificate']['integer_lower'],5)


if __name__=='__main__':
    unittest.main()
