import json
from pathlib import Path
import unittest
import numpy as np
from audit_eight_prime_anchor_branch import explicit_capacity
from audit_compressed_tail_branch import local_minima
from portable_tail_certificates import CoreCapacityUpper


class CompressedTailBranchTests(unittest.TestCase):
    def test_independent_point_weights_match_all_four_anchor_boxes(self):
        for p in ([0,0,0],[0,0,1],[0,1,0],[0,1,1]):
            np.testing.assert_array_equal(CoreCapacityUpper(p).compression,explicit_capacity(p))

    def test_new_local_auditor_reproduces_existing_integer_certificate(self):
        data = json.loads(Path('data/aggregate-conditional-tail-eight-lp.json').read_text())
        actual,checked = local_minima(data['certificate'],data['high_pattern'])
        self.assertEqual(actual,data['certificate']['point_multipliers'])
        self.assertEqual(checked,59175)
        bad = json.loads(json.dumps(data['certificate']))
        bad['tail_multipliers'][0]['multiplier'] = 1
        with self.assertRaises(AssertionError):
            local_minima(bad,data['high_pattern'])


if __name__=='__main__':
    unittest.main()
