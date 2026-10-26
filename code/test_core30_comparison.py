"""Checks of core geometry and the stored arithmetic certificates."""
import hashlib,itertools,json,unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

class CoreComparisonTests(unittest.TestCase):
    def test_core_geometries_and_capacity(self):
        points=list(itertools.product(range(2),range(3),range(5)))
        for a,b in [(0,1),(1,0),(1,1)]:
            u=[x for x in points if not(x[1]==0 and x[2]==0 or x==(0,a,b))]
            self.assertEqual(len(u),27)
            self.assertEqual(sorted(sum(x[0]==i for x in u) for i in range(2)),[13,14])
            total=0
            for mask in range(8):
                axes=[i for i in range(3) if mask>>i&1];counts={}
                for x in u:
                    key=tuple(x[i] for i in axes);counts[key]=counts.get(key,0)+1
                total+=max(counts.values())
            self.assertEqual(total,68)

    def test_envelope_and_interpolation(self):
        h=[68,41,20,12,4,3,2,1,0];masses=[6,13,7,1];values=[1,2,4,8]
        for j in range(8):
            for part in range(5):
                t=F(j)+F(part,5)
                actual=sum(w*max(F(0),f-t) for f,w in zip(values,masses))
                upper=h[j]*(1-F(part,5))+h[j+1]*F(part,5)
                self.assertEqual(actual,upper)
        self.assertEqual(sum(masses),27)
        self.assertEqual(sum(w*3**k for k,w in enumerate(masses)),135)

    def test_exhaustive_audit_matches_primary_rows(self):
        a=json.loads((ROOT/'data/audit-core30-no10-envelope.json').read_text())
        b=json.loads((ROOT/'data/core30-comparison-d5-0.json').read_text())
        self.assertTrue(a['passed']);self.assertEqual(a['n_configurations'],810000)
        for case in [21,26,27]:
            self.assertEqual(b['rows'][case]['stoploss_numerators'],a['stoploss_numerators'])
            self.assertEqual(b['rows'][case]['law_numerators'],a['law_numerators'])

    def test_certificates_and_schedule_identity(self):
        raw=(ROOT/'data/schedule-core30-no10.json').read_bytes();digest=hashlib.sha256(raw).hexdigest()
        a=json.loads((ROOT/'data/cert-core30-no10.json').read_text())
        b=json.loads((ROOT/'data/audit-core30-no10-continuation.json').read_text())
        self.assertEqual(a['schedule_sha256'],digest);self.assertEqual(b['schedule_sha256'],digest)
        x=F(a['total_integer'])/F(a['scale']);y=F(b['normalized_upper'])
        self.assertLess(x,F(989676,1000000));self.assertLess(y,F(989676,1000000))
        self.assertLess(abs(x-y),F(1,10**12))
        self.assertEqual(a['parameters']['count_cap'],120);self.assertEqual(b['parameters']['count_cap'],100)

if __name__=='__main__':unittest.main()
