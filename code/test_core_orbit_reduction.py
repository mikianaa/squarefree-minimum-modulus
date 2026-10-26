import itertools
import random
import unittest
from core_orbit_reduction import CoreOrbits,minimal_orbits,SubsetIndex


def ordinary_images(points,mask):
    index = {point:j for j,point in enumerate(points)}
    active = [point for j,point in enumerate(points) if mask>>j&1]
    images = set()
    for p5 in itertools.permutations((2,3,4)):
        labels5 = {1:1,**dict(zip((2,3,4),p5))}
        for p7 in itertools.permutations((2,3,4,5,6)):
            labels7 = {1:1,**dict(zip((2,3,4,5,6),p7))}
            image = sum(1<<index[e,(a[0],labels5[a[1]],labels7[a[2]])] for e,a in active)
            images.add(image)
    return images


class CoreOrbitTests(unittest.TestCase):
    def test_subset_index_against_direct_search(self):
        rng = random.Random(12)
        index = SubsetIndex(10)
        saved = []
        for batch in range(4):
            added = [rng.randrange(1<<10) for _ in range(30)]
            saved.extend(added)
            index.extend(added)
            for mask in range(1<<10):
                self.assertEqual(index.contains_subset(mask),any(s&mask==s for s in saved))

    def test_orbits_against_independent_coordinate_relabelling(self):
        rng = random.Random(210)
        for pattern in ((0,0,0),(0,0,1),(0,1,0),(0,1,1)):
            model = CoreOrbits(pattern)
            masks = [0,(1<<len(model.points))-1]+[rng.getrandbits(len(model.points)) for _ in range(3)]
            for mask in masks:
                expected = ordinary_images(model.points,mask)
                self.assertEqual(model.orbit(mask),expected)
                canonical = model.canonical([mask])[0]
                self.assertIn(canonical,expected)
                self.assertEqual(set(model.canonical(expected)),{canonical})

    def test_inclusion_reduction_covers_every_input_orbit(self):
        rng = random.Random(511)
        model = CoreOrbits((0,1,1))
        masks = [rng.getrandbits(len(model.points)) for _ in range(20)]
        masks += [m|rng.getrandbits(len(model.points)) for m in masks]
        canonical = model.canonical(masks)
        kept = minimal_orbits(model,canonical)
        images = set().union(*(ordinary_images(model.points,m) for m in kept))
        self.assertTrue(all(any(x&m==x for x in images) for m in canonical))
        self.assertTrue(set(kept)<=set(canonical))


if __name__=='__main__':
    unittest.main()
