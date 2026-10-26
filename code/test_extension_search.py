"""Check actual congruence reconstruction from the compressed search grid."""
import json
from pathlib import Path
import unittest
from search_core_extensions import search,crt_extend


class ExtensionSearchTests(unittest.TestCase):
    def test_crt_reconstruction(self):
        a,d = crt_extend(13,30,[4,7],[11,13])
        self.assertEqual(d,4290)
        self.assertEqual((a%30,a%11,a%13),(13,4,7))
        self.assertTrue(0<=a<d)

    def test_one_greedy_pass_preserves_actual_family_bound(self):
        root = Path(__file__).resolve().parents[1]
        source = json.loads((root/'data/min-uncovered-30030-witness.json').read_text())
        out = search(source,[11,13],[1,1],seconds=0,seed=3)
        self.assertLessEqual(out['uncovered'],source['uncovered'])
        family = [(r['modulus'],r['residue']) for r in out['system']]
        count = sum(not any(n%d==a for d,a in family) for n in range(30030))
        self.assertEqual(count,out['uncovered'])
        self.assertEqual(out['pattern'],[0,1,1,1,1])


if __name__=='__main__':
    unittest.main()
