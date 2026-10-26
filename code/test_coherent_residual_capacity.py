import math
import unittest

from coherent_residual_capacity import cap_over_coherent, divisors_gt1


class TestCoherentResidualCapacity(unittest.TestCase):
    def test_brute_force_small_L(self):
        L = 2310
        divs = [d for d in range(2, L + 1) if L % d == 0 and d > 1]
        for m in (2, 3, 6, 30, 105, 210, 2310):
            R = set(range(0, L, m))
            brute = sum(max(sum(1 for x in R if x % d == a) for a in range(d))
                        for d in divs)
            cap, Rsize = cap_over_coherent(m, L, divs)
            self.assertEqual(Rsize, len(R))
            self.assertEqual(cap, brute, msg=f"m={m}")

    def test_known_nine_prime_values(self):
        cap, Rsize = cap_over_coherent(210)
        self.assertEqual(Rsize, 1062347)
        self.assertEqual(cap, 22161973)
        cap, Rsize = cap_over_coherent(223092870)
        self.assertEqual(Rsize, 1)
        self.assertEqual(cap, 511)
        self.assertEqual(len(divisors_gt1()), 511)


if __name__ == "__main__":
    unittest.main()
