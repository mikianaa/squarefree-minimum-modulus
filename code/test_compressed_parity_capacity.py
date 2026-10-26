"""Exhaustive checks of compressed parity bounds, independent of their weights."""
from collections import Counter
import itertools
import unittest
from compressed_parity_capacity import capacity,extend_total


def subsets(points):
    return [tuple(p for i,p in enumerate(points) if mask>>i&1) for mask in range(1<<len(points))]


def direct_capacity(layers,k):
    total = 0
    for s in range(1,1<<k):
        cs = [Counter(tuple(p[i] for i in range(k) if s>>i&1) for p in layer) for layer in layers]
        total += max((cs[0]+cs[1]).values(),default=0)
        total += max((max(c.values(),default=0) for c in cs),default=0)
    return total


class CompressedParityTests(unittest.TestCase):
    def test_exact_for_all_small_layer_sizes(self):
        for dims in [((1,2),(2,1)),((2,2),(2,2)),((1,3),(2,2))]:
            options = [subsets(list(itertools.product(*(range(d) for d in ds)))) for ds in dims]
            best = {}
            for layers in itertools.product(*options):
                key = tuple(map(len,layers))
                best[key] = max(best.get(key,0),direct_capacity(layers,2))
            for counts,v in best.items():
                self.assertEqual(capacity(dims,counts),v,(dims,counts))

    def test_extension_bound_on_all_small_actual_sets(self):
        dims = ((1,3),(2,2))
        options = [subsets(list(itertools.product(*(range(d) for d in ds)))) for ds in dims]
        # A total-capacity deficit can be weaker than separately clipping the
        # layers. Enumerate every even-shape allocation directly for the latter.
        actual = {b:[None]*8 for b in (0,1)}
        for layers in itertools.product(*options):
            counts = list(map(len,layers))
            caps = []
            for s in range(1,4):
                caps.append([max(Counter(tuple(p[i] for i in range(2) if s>>i&1) for p in layer).values(),default=0)
                             for layer in layers])
            for b in (0,1):
                best = min(sum(max(0,(4-(b==e))*counts[e]-sum(c[e]*(1+(a==e)) for c,a in zip(caps,allocation)))
                               for e in (0,1)) for allocation in itertools.product((0,1),repeat=3))
                r = sum(counts)
                old = actual[b][r]
                actual[b][r] = best if old is None else min(old,best)
        for b in (0,1):
            for lower in range(8):
                self.assertEqual(extend_total(dims,lower,5,b)['uncovered_lower'],min(actual[b][lower:]))


if __name__=='__main__':
    unittest.main()
