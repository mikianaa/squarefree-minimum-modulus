"""Union accounting at C0 = 8: how much of the stage charges sits on fibres that are already covered.

Identity.  B_k (the set covered at stage k) is determined by the first k coordinates, and every later
stage preserves the law of the first k coordinates, so P_n(B_k) = P_k(B_k) and

    P_n( U_k B_k ) = sum_k P_k( B_k minus U_{l<k} B_l )
                   = sum_k sum_{x : fibre x not yet covered} (alpha_k(x) - delta_k)_+ P_{k-1}(x) / (1 - delta_k).

The certificate charges sum_k P_k(B_k) instead (the union bound); the difference is the charge sitting on
fibres that are already covered ("double counting").  This script measures it exactly for the small stages.

Part A (exhaustive, hence a valid bound over ALL stage-5 and stage-7 families): max of the union objective
    U = loss_5 + sum_{x in Z_30 \\ B_5} loss_7(x)
compared with the max of loss_5 + loss_7 (exact_two_stage.py: 0.211779).  Same enumeration and symmetry
reduction as exact_two_stage.py; exact integer / Fraction arithmetic.

Part B (specific families, NOT a bound over all families): exact stage-11 charge and the part of it on fibres
of Z_210 already covered by stages 5 and 7, for
  (i)  the nested/aligned stage-5 and stage-7 families of stage11.py with the recorded CP-SAT stage-11
       adversary (data/stage11-relaxation-tower-C0-8.txt, 0.055092), and
  (ii) the two-stage maximiser of exact_two_stage.py (canonical table [0,0,1,2,2,19], stage-7 residues
       [1,2,2,1,3,2,2]);
with the best ALIGNED stage-11 adversary for both objectives (all fibres / uncovered fibres only), and, with
--cpsat SECONDS (needs OR-Tools), the CP-SAT relaxation of stage11.py run with both weight vectors.

Usage: python code/union_accounting.py [D5 D7 D11] [--cpsat SECONDS] [out.json]   (deltas D/10^6; default 195000 142500 363750)
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_two_stage import (L, X30, stage5_tables, group_elements, canon, p5_weights, loss5, stage7_counts)
from stage11 import fibre_weights, frac, g_fn, aligned_stage11, cpsat_relaxed_stage11, M7, DIV210

STAGE11_RESIDUES = {1: 0, 2: 0, 3: 0, 5: 4, 6: 3, 7: 4, 10: 4, 14: 4, 15: 9, 21: 18, 30: 9, 35: 4, 42: 39,
                    70: 4, 105: 54, 210: 165}     # recorded CP-SAT adversary, objective 0.055091566626


def part_a(D5, D7):
    t0 = time.time()
    fams = stage5_tables(); G = group_elements()
    reps = {}
    for o2, o3, o6, cov in fams:
        reps.setdefault(canon(cov, G), (o2, o3, o6))
    N = stage7_counts()
    gnum = np.array([max(0, n * L - 7 * D7) for n in range(8)], dtype=np.int64)
    Gm = gnum[N.astype(np.int64)]
    den7 = 7 * (L - D7)
    rows = []
    for tab, fam in reps.items():
        w = p5_weights(tab, D5)
        covered = [bool(tab[X30[x][0] * 3 + X30[x][1]] >> X30[x][2] & 1) for x in range(30)]
        wu = [Fraction(0) if covered[x] else w[x] for x in range(30)]
        den = 1
        for v in w: den = den * v.denominator // np.gcd(den, v.denominator)
        ws = np.array([int(v * den) for v in w], dtype=np.int64)
        wus = np.array([int(v * den) for v in wu], dtype=np.int64)
        assert max(int(v) for v in ws) * int(gnum.max()) * 30 < 2 ** 62
        vals = Gm @ ws; uvals = Gm @ wus
        i_all = int(np.argmax(vals)); i_unc = int(np.argmax(uvals))
        l5 = loss5(tab, D5)
        l7_all = Fraction(int(vals[i_all]), den * den7)
        l7_unc_max = Fraction(int(uvals[i_unc]), den * den7)
        l7_unc_at_allmax = Fraction(int(uvals[i_all]), den * den7)
        cfg = np.unravel_index(i_unc, (2, 3, 5, 6, 10, 15, 30))
        chk = Fraction(0)
        for x in range(30):
            if covered[x]: continue
            n = sum(1 for m, a in zip(M7, cfg) if x % m == a)
            chk += w[x] * max(Fraction(0), Fraction(n, 7) - frac(D7)) / (1 - frac(D7))
        assert chk == l7_unc_max
        rows.append(dict(table=list(tab), family=fam, covered_mass_after_5=float(sum(w[x] for x in range(30) if covered[x])),
                         loss5=float(l5), loss7_max=float(l7_all), per_stage_sum=float(l5 + l7_all),
                         loss7_on_uncovered_at_per_stage_maximiser=float(l7_unc_at_allmax),
                         double_counting_at_per_stage_maximiser=float(l7_all - l7_unc_at_allmax),
                         loss7_uncovered_max=float(l7_unc_max), union_objective=float(l5 + l7_unc_max),
                         union_objective_exact=str(l5 + l7_unc_max),
                         union_argmax_stage7_residues=[int(c) for c in cfg]))
    rows.sort(key=lambda r: -r["union_objective"])
    best_sum = max(r["per_stage_sum"] for r in rows)
    return dict(D5=D5, D7=D7, n_orbits=len(reps), per_stage_max=best_sum, union_max=rows[0]["union_objective"],
                union_max_exact=rows[0]["union_objective_exact"], gain=best_sum - rows[0]["union_objective"],
                best=rows[0], all_orbits=rows, seconds=round(time.time() - t0, 1))


def build_from_table(tab, F7, D5, D7):
    """P_7 on Z_210, covered sets after stage 5 and after stage 7, loss5, loss7, from a stage-5 covered-set
    table (index x2*3+x3 -> bitmask over Z_5) and a stage-7 family F7: m -> (residue mod m, class mod 7)."""
    d7 = frac(D7)
    P5 = p5_weights(tab, D5)
    cov5 = {x: bool(tab[X30[x][0] * 3 + X30[x][1]] >> X30[x][2] & 1) for x in range(30)}
    l5 = loss5(tab, D5)
    P7 = {}; cov7 = {}; l7 = Fraction(0)
    for x in range(30):
        mask = 0
        for m, (a, c) in F7.items():
            if x % m == a: mask |= 1 << c
        a = Fraction(bin(mask).count("1"), 7)
        if a > d7: l7 += P5[x] * (a - d7) / (1 - d7)
        w = fibre_weights(7, mask, D7)
        for y in range(210):
            if y % 30 == x:
                P7[y] = P5[x] * w[y % 7]
                cov7[y] = cov5[x] or bool(mask >> (y % 7) & 1)
    assert sum(P7.values()) == 1
    return P7, cov7, l5, l7


def table_from_F5(F5):
    tab = [0] * 6
    for x2 in range(2):
        for x3 in range(3):
            mask = 0
            for m, (a, c) in F5.items():
                if (x2 if m == 2 else x3 if m == 3 else (x2 if False else [x for x in range(6) if x % 2 == x2 and x % 3 == x3][0])) % m == a:
                    mask |= 1 << c
            tab[x2 * 3 + x3] = mask
    return tuple(tab)


def part_b(label, tab, F7, D5, D7, D11, residues=None, cpsat_seconds=None):
    P7, cov7, l5, l7 = build_from_table(tab, F7, D5, D7)
    P7u = {y: (Fraction(0) if cov7[y] else w) for y, w in P7.items()}
    covered_mass = sum(w for y, w in P7.items() if cov7[y])
    g = g_fn(D11)
    out = dict(label=label, table=list(tab), F7={str(m): v for m, v in F7.items()}, D5=D5, D7=D7, D11=D11,
               loss5=float(l5), loss7=float(l7), per_stage_sum_5_7=float(l5 + l7),
               covered_mass_after_7=float(covered_mass), double_counting_5_7=float(l5 + l7 - covered_mass))
    if residues:
        n = {y: sum(1 for m, a in residues.items() if y % m == a) for y in range(210)}
        ch = {y: P7[y] * g(min(Fraction(1), Fraction(n[y], 11))) for y in range(210)}
        tot = sum(ch.values()); on_cov = sum(c for y, c in ch.items() if cov7[y])
        out["recorded_stage11_adversary"] = dict(residues=residues, charge=float(tot), charge_on_covered=float(on_cov),
                                                 charge_on_uncovered=float(tot - on_cov),
                                                 fraction_double_counted=float(on_cov / tot) if tot else None)
    b_all, t_all, _ = aligned_stage11(P7, D11)
    b_unc, t_unc, _ = aligned_stage11(P7u, D11)
    out["aligned_stage11"] = dict(best_all_fibres=float(b_all), alignment_all=t_all,
                                  best_uncovered_only=float(b_unc), alignment_uncovered=t_unc)
    if cpsat_seconds:
        out["cpsat_stage11_all"] = cpsat_relaxed_stage11(P7, D11, time_limit=cpsat_seconds)
        out["cpsat_stage11_uncovered_only"] = cpsat_relaxed_stage11(P7u, D11, time_limit=cpsat_seconds)
        for k in ("cpsat_stage11_all", "cpsat_stage11_uncovered_only"):
            r = out[k]
            if r.get("residues"):
                n = {y: sum(1 for m, a in r["residues"].items() if y % m == a) for y in range(210)}
                ch = {y: P7[y] * g(min(Fraction(1), Fraction(n[y], 11))) for y in range(210)}
                tot = sum(ch.values()); on_cov = sum(c for y, c in ch.items() if cov7[y])
                r["exact_charge_all_fibres"] = float(tot); r["exact_charge_on_covered"] = float(on_cov)
                r["exact_charge_on_uncovered"] = float(tot - on_cov)
    return out


def main():
    argv = sys.argv[1:]
    cps = None
    if "--cpsat" in argv:
        i = argv.index("--cpsat"); cps = float(argv[i + 1]); del argv[i:i + 2]
    out = [a for a in argv if a.endswith(".json")]
    nums = [int(a) for a in argv if not a.endswith(".json")]
    D5, D7, D11 = nums if len(nums) == 3 else (195000, 142500, 363750)
    A = part_a(D5, D7)
    print(json.dumps({k: v for k, v in A.items() if k not in ("all_orbits", "best")}, indent=1))
    print("best union family:", A["best"]["family"], "table", A["best"]["table"], "stage-7 residues", A["best"]["union_argmax_stage7_residues"])
    for r in A["all_orbits"][:5]:
        print(f"  union={r['union_objective']:.9f} per_stage={r['per_stage_sum']:.9f} loss5={r['loss5']:.6f} "
              f"loss7_all={r['loss7_max']:.6f} loss7_unc={r['loss7_uncovered_max']:.6f} "
              f"dc_at_per_stage_max={r['double_counting_at_per_stage_maximiser']:.6f}")
    # (i) stage11.py's families: nested stage 5 at residue 0, aligned stage 7 at residue 0
    tab_i = table_from_F5({2: (0, 0), 3: (0, 1), 6: (0, 2)})
    F7_i = {m: (0, i) for i, m in enumerate(M7)}
    B1 = part_b("aligned families of stage11.py", tab_i, F7_i, D5, D7, D11, residues=STAGE11_RESIDUES, cpsat_seconds=cps)
    print(json.dumps(B1, indent=1, default=str))
    # (ii) the two-stage maximiser (canonical table and its stage-7 argmax from exact_two_stage.py)
    rec = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data",
                                      "exact-two-stage-C0-8-D5-195000-D7-142500.json")))
    tab_ii = tuple(rec["best"]["table"]); cfg = rec["best"]["argmax_stage7_residues"]
    F7_ii = {m: (int(a), i) for i, (m, a) in enumerate(zip(M7, cfg))}
    B2 = part_b("two-stage maximiser of exact_two_stage.py", tab_ii, F7_ii, D5, D7, D11, cpsat_seconds=cps)
    assert abs(B2["per_stage_sum_5_7"] - rec["two_stage_max_float"]) < 1e-12
    print(json.dumps(B2, indent=1, default=str))
    if out:
        json.dump(dict(part_a=A, part_b=[B1, B2]), open(out[0], "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
