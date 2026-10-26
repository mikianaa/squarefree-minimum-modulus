import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from audit_group_event_cut import verify_bound
from audit_multitail_identity import audit_group_cut
from transfer_group_event_bounds import transfer


class TransferGroupEventTests(unittest.TestCase):
    def test_all_transferred_formulas_have_identical_resource_domains(self):
        data = json.loads(Path('data/nine-prime-group-refinement-round3-model.json').read_text())
        source,cuts = data['input'],data['result']['group_event_constraints']
        expanded = transfer(source,cuts)
        self.assertEqual(len(expanded),327)
        self.assertEqual(len(transfer(source,expanded)),327)
        for cut in expanded:
            self.assertEqual(cut['axis']<5,cut['transferred_from_axis']<5)
            self.assertEqual(cut['proof_audit']['axis'],cut['axis'])

    def test_actual_proof_reuse_requires_identical_inputs_and_resource_set(self):
        source = dict(core_residual=[0,1],high_primes=[11,13,17])
        shapes = [d for d in range(3,211) if 210%d==0]
        cut = dict(axis=0,bound=1,events=[dict(point=x,threshold=14,weight=1) for x in (0,1)])
        temporary_root = Path('tmp').resolve()
        temporary_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=temporary_root) as directory:
            assert Path(directory).resolve().is_relative_to(temporary_root)
            cut['proof_audit'] = verify_bound(source['core_residual'],shapes,cut,Path(directory)/'proof.gz')
            self.assertTrue(cut['proof_audit']['verified'])
            cache = set()
            audit_group_cut(source,cut,cache)
            self.assertEqual(len(cache),1)
            moved = transfer(source,[cut])[1]
            with patch('sat_uncovered_audit.check_rup',side_effect=AssertionError('unexpected replay')):
                self.assertEqual(audit_group_cut(source,moved,cache),0)
            bad = copy.deepcopy(moved)
            bad['events'][0]['threshold'] = 13
            with self.assertRaises(AssertionError):
                audit_group_cut(source,bad,cache)
            bad = copy.deepcopy(moved)
            bad['axis'] = bad['proof_audit']['axis'] = 3
            with self.assertRaises(AssertionError):
                audit_group_cut(source,bad,cache)
            self.assertEqual(len(cache),1)


if __name__=='__main__':
    unittest.main()
