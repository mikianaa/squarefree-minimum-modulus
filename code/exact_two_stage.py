"""Exact worst case of loss_5 + loss_7 over ALL families at cutoff C0 = 8 (squarefree, distinct moduli).

Stage 5 (prime 5): allowed moduli 10, 15, 30 (modulus 5 <= C0 excluded).  Measure P_2 is uniform on
Z_2 x Z_3 (delta_2 = delta_3 = 0, no loss at p = 2, 3 since 2, 3, 6 <= C0).
Stage 7: allowed moduli 7m, m in {2,3,5,6,10,15,30} (7 itself excluded).  Measure P_5 = BBMST fibre
construction from the ACTUAL stage-5 family with distortion delta_5 = D5/L:
    fibre (x2,x3), covered set C (alpha = |C|/5):
      alpha <= delta: uncovered weight 1/(5(1-alpha)), covered 0
      alpha >  delta: uncovered 1/(5(1-delta)),  covered (alpha-delta)/(alpha(1-delta))/5     (times P_2(x2,x3) = 1/6)
Stage-7 loss = sum_x P_5(x) (alpha_7(x) - d7)_+/(1-d7), alpha_7(x) = |{classes mod 7 of active moduli}|/7.
Monotonicity: for a FIXED P_5 adding a stage-7 modulus or separating two colliding classes never
decreases alpha_7 pointwise, so the stage-7 maximum is attained with all 7 moduli present and 7
distinct classes mod 7 (possible: 7 classes), i.e. alpha_7(x) = n(x)/7, n(x) = #{m : x = a_m mod m}.
That leaves prod m = 810000 residue configurations, enumerated exhaustively.
Stage 5 is enumerated over all 11*16*31 = 5456 families (absent / residue mod m / class mod 5), reduced
to orbits of the covered-set table under Aff(Z_2) x Aff(Z_3) x Aff(Z_5) (order 240), which is a symmetry
of the whole problem (progressions with squarefree moduli are coordinate cylinders).
All arithmetic is exact: integers (numpy int64, magnitudes checked) and Fractions.

Usage: python exact_two_stage.py D5 D7 [out.json]       (deltas D/10^6)
"""
import itertools, json, sys, time
from fractions import Fraction
import numpy as np

L = 10 ** 6
X30 = [(x % 2, x % 3, x % 5) for x in range(30)]           # index x in Z_30 <-> coordinates (CRT)


def stage5_tables():
    """All stage-5 families -> covered-set table cov[(x2,x3)] = bitmask over Z_5."""
    opts2 = [None] + [(a, c) for a in range(2) for c in range(5)]        # modulus 10: residue mod 2, class mod 5
    opts3 = [None] + [(a, c) for a in range(3) for c in range(5)]        # modulus 15
    opts6 = [None] + [(a, c) for a in range(6) for c in range(5)]        # modulus 30: residue mod 6
    fams = []
    for o2, o3, o6 in itertools.product(opts2, opts3, opts6):
        cov = {}
        for x2 in range(2):
            for x3 in range(3):
                mask = 0
                if o2 and x2 == o2[0]: mask |= 1 << o2[1]
                if o3 and x3 == o3[0]: mask |= 1 << o3[1]
                if o6 and x2 == o6[0] % 2 and x3 == o6[0] % 3: mask |= 1 << o6[1]
                cov[(x2, x3)] = mask
        fams.append((o2, o3, o6, cov))
    return fams


def group_elements():
    """Aff(Z_2) x Aff(Z_3) x Aff(Z_5): maps on coordinates."""
    aff = lambda p: [(a, b) for a in range(1, p) for b in range(p)]
    return list(itertools.product(aff(2), aff(3), aff(5)))


def canon(cov, G):
    best = None
    for (g2, g3, g5) in G:
        # image table: new coordinates (g2(x2), g3(x3)), new classes g5(c)
        tab = [0] * 6
        for (x2, x3), mask in cov.items():
            y2 = (g2[0] * x2 + g2[1]) % 2; y3 = (g3[0] * x3 + g3[1]) % 3
            nm = 0
            for c in range(5):
                if mask >> c & 1: nm |= 1 << ((g5[0] * c + g5[1]) % 5)
            tab[y2 * 3 + y3] = nm
        t = tuple(tab)
        if best is None or t < best: best = t
    return best


def p5_weights(tab, D5):
    """Exact P_5 on Z_30 for covered-set table tab (index x2*3+x3 -> mask), as Fractions."""
    d = Fraction(D5, L)
    w = []
    for x in range(30):
        x2, x3, x5 = X30[x]
        mask = tab[x2 * 3 + x3]; j = bin(mask).count("1"); a = Fraction(j, 5)
        covered = mask >> x5 & 1
        if a <= d:
            v = Fraction(0) if covered else 1 / (5 * (1 - a))
        else:
            v = (a - d) / (a * (1 - d) * 5) if covered else 1 / (5 * (1 - d))
        w.append(v / 6)
    assert sum(w) == 1
    return w


def loss5(tab, D5):
    d = Fraction(D5, L); tot = Fraction(0)
    for i in range(6):
        a = Fraction(bin(tab[i]).count("1"), 5)
        if a > d: tot += (a - d) / (1 - d) / 6
    return tot


def stage7_counts():
    """n(x) for every residue configuration of the 7 stage-7 moduli: array (810000, 30) int8."""
    N = np.zeros((1, 30), dtype=np.int8)
    for m in (2, 3, 5, 6, 10, 15, 30):
        I = np.array([[1 if x % m == a else 0 for x in range(30)] for a in range(m)], dtype=np.int8)
        N = (N[:, None, :] + I[None, :, :]).reshape(-1, 30)
    assert N.shape == (810000, 30)
    return N


def main():
    D5 = int(sys.argv[1]); D7 = int(sys.argv[2]); out = sys.argv[3] if len(sys.argv) > 3 else None
    t0 = time.time()
    fams = stage5_tables(); G = group_elements(); assert len(G) == 240
    reps = {}
    for o2, o3, o6, cov in fams:
        key = canon(cov, G)
        reps.setdefault(key, (o2, o3, o6))
    print(f"{len(fams)} stage-5 families, {len(reps)} orbits of covered-set tables ({time.time()-t0:.1f}s)", flush=True)
    N = stage7_counts()
    # g(n/7) = (n/7 - d7)_+/(1-d7) = (n L - 7 D7)_+ / (7 (L - D7))
    gnum = np.array([max(0, n * L - 7 * D7) for n in range(8)], dtype=np.int64)
    Gm = gnum[N.astype(np.int64)]                                   # (810000, 30)
    den7 = 7 * (L - D7)
    results = []
    best = None
    for tab, fam in reps.items():
        w = p5_weights(tab, D5)
        # common denominator for the weights
        den = 1
        for v in w: den = den * v.denominator // np.gcd(den, v.denominator)
        ws = np.array([int(v * den) for v in w], dtype=object)
        assert all(int(v * den) == v * den for v in w)
        wmax = max(int(v) for v in ws)
        assert wmax * int(gnum.max()) * 30 < 2 ** 62, "int64 overflow guard"
        vals = Gm @ ws.astype(np.int64)                             # exact integers
        imax = int(np.argmax(vals)); vmax = int(vals[imax])
        l7 = Fraction(vmax, den * den7)
        # independent exact re-evaluation of the argmax configuration with Fractions
        cfg = np.unravel_index(imax, (2, 3, 5, 6, 10, 15, 30))
        chk = Fraction(0)
        for x in range(30):
            n = sum(1 for m, a in zip((2, 3, 5, 6, 10, 15, 30), cfg) if x % m == a)
            assert n == int(N[imax, x])
            chk += w[x] * max(Fraction(0), Fraction(n, 7) - Fraction(D7, L)) / (1 - Fraction(D7, L))
        assert chk == l7
        l5 = loss5(tab, D5)
        tot = l5 + l7
        results.append(dict(table=list(tab), family=fam, loss5=str(l5), loss7_max=str(l7), total=str(tot),
                            total_float=float(tot), loss5_float=float(l5), loss7_float=float(l7),
                            argmax_stage7_residues=[int(c) for c in cfg]))
        if best is None or tot > best[0]: best = (tot, results[-1])
    results.sort(key=lambda r: -r["total_float"])
    summary = dict(C0=8, D5=D5, D7=D7, L=L, n_families=len(fams), n_orbits=len(reps), n_stage7_configs=int(N.shape[0]),
                   two_stage_max=str(best[0]), two_stage_max_float=float(best[0]), best=best[1],
                   top10=results[:10], seconds=round(time.time() - t0, 1))
    print(json.dumps({k: v for k, v in summary.items() if k != "top10"}, indent=1))
    for r in results[:10]:
        print(f"  total={r['total_float']:.9f} loss5={r['loss5_float']:.6f} loss7={r['loss7_float']:.6f} fam={r['family']} tab={r['table']}")
    if out: json.dump(summary, open(out, "w"), indent=1)
    return summary


if __name__ == "__main__":
    main()
