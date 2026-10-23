"""Stage 11 at C0 = 8: exact P_7 for explicit stage-5 / stage-7 families and
 (a) the ALIGNED stage-11 adversary (all 16 moduli 11m, m | 210, cylinders aligned at a point of Z_210,
     classes mod 11 chosen so that alpha_11(y) = min(1, n(y)/11): pairs (1,210),(2,105),(3,70),(5,42),(7,30)
     share a class, the other six moduli get their own) -> exact loss_11, a LOWER bound of the stage-11 max;
 (b) the independent-model (comparison) value of stage 11 with the same deltas;
 (c) optionally a CP-SAT relaxation max_a sum_y P_7(y) g(min(1, n(y)/11)) over all residue configurations
     (an UPPER bound of the exact stage-11 max for this P_7).
Everything in (a), (b) is exact (Fractions).
"""
import itertools, json, sys, time
from fractions import Fraction
import numpy as np

L = 10 ** 6
M7 = (2, 3, 5, 6, 10, 15, 30)
DIV210 = [m for m in range(1, 211) if 210 % m == 0]
assert len(DIV210) == 16


def frac(D): return Fraction(D, L)


def fibre_weights(p, covered_mask, D):
    """BBMST fibre rule on Z_p: list of p Fractions (relative to parent weight), covered set as bitmask."""
    d = frac(D); j = bin(covered_mask).count("1"); a = Fraction(j, p)
    w = []
    for y in range(p):
        cov = covered_mask >> y & 1
        if a <= d:
            w.append(Fraction(0) if cov else 1 / (p * (1 - a)))
        else:
            w.append((a - d) / (a * (1 - d) * p) if cov else 1 / (p * (1 - d)))
    assert sum(w) == 1
    return w


def build_P7(F5, F7, D5, D7):
    """F5: dict m -> (residue mod m, class mod 5) for m in {2,3,6} (modulus 5m); F7: dict m -> (residue mod m,
    class mod 7) for m in M7.  Returns P_7 as dict y in Z_210 -> Fraction, plus loss5, loss7."""
    d5, d7 = frac(D5), frac(D7)
    P5 = {}; loss5 = Fraction(0)
    for x in range(30):
        mask = 0
        for m, (a, c) in F5.items():
            if x % m == a: mask |= 1 << c
        w = fibre_weights(5, mask, D5)
        P5[x] = w[x % 5] / 6
    for x in range(6):
        mask = 0
        for m, (a, c) in F5.items():
            if x % m == a: mask |= 1 << c
        a = Fraction(bin(mask).count("1"), 5)
        if a > d5: loss5 += (a - d5) / (1 - d5) / 6
    assert sum(P5.values()) == 1
    P7 = {}; loss7 = Fraction(0)
    for x in range(30):
        mask = 0
        for m, (a, c) in F7.items():
            if x % m == a: mask |= 1 << c
        a = Fraction(bin(mask).count("1"), 7)
        if a > d7: loss7 += P5[x] * (a - d7) / (1 - d7)
        w = fibre_weights(7, mask, D7)
        for y in range(210):
            if y % 30 == x: P7[y] = P5[x] * w[y % 7]
    assert sum(P7.values()) == 1
    return P7, loss5, loss7


def g_fn(D):
    d = frac(D)
    return lambda a: (a - d) / (1 - d) if a > d else Fraction(0)


def aligned_stage11(P7, D11, alignments=None):
    """Exact loss_11 for aligned stage-11 families; returns (best value, best alignment, all values)."""
    g = g_fn(D11); best = (Fraction(-1), None); vals = {}
    for t in (range(210) if alignments is None else alignments):
        tot = Fraction(0)
        for y, w in P7.items():
            if w == 0: continue
            k = sum(1 for p in (2, 3, 5, 7) if y % p == t % p)        # aligned coordinates
            n = 2 ** k                                                 # active moduli = divisors of the product
            tot += w * g(min(Fraction(1), Fraction(n, 11)))
        vals[t] = tot
        if tot > best[0]: best = (tot, t)
    return best[0], best[1], vals


def model_stage11(D5, D7, D11):
    """Comparison value: independent V_2 ~ B(1/2), V_3 ~ B(1/3), V_5 ~ B(nu5/5), V_7 ~ B(nu7/7); S = 2^N/11."""
    g = g_fn(D11); nu5 = 1 / (1 - frac(D5)); nu7 = 1 / (1 - frac(D7))
    prs = [Fraction(1, 2), Fraction(1, 3), nu5 / 5, nu7 / 7]
    tot = Fraction(0)
    for bits in itertools.product((0, 1), repeat=4):
        pr = Fraction(1)
        for b, q in zip(bits, prs): pr *= q if b else 1 - q
        tot += pr * g(Fraction(2 ** sum(bits), 11))
    return tot


def model_stage11_truncated(D5, D7, D11):
    """Same but with S capped at 1 (NOT a valid bound in general; diagnostic only)."""
    g = g_fn(D11); nu5 = 1 / (1 - frac(D5)); nu7 = 1 / (1 - frac(D7))
    prs = [Fraction(1, 2), Fraction(1, 3), nu5 / 5, nu7 / 7]
    tot = Fraction(0)
    for bits in itertools.product((0, 1), repeat=4):
        pr = Fraction(1)
        for b, q in zip(bits, prs): pr *= q if b else 1 - q
        tot += pr * g(min(Fraction(1), Fraction(2 ** sum(bits), 11)))
    return tot


def cpsat_relaxed_stage11(P7, D11, time_limit=600, workers=8):
    """max over residue configurations a_m (m | 210) of sum_y P7(y) g(min(1, n(y)/11)), n(y) = #{m: y = a_m mod m}.
    Upper bound of the exact stage-11 maximum for this P7 (classes ignored, alpha <= min(1,n/11))."""
    from ortools.sat.python import cp_model
    d = frac(D11)
    # objective scaled to integers: weights P7(y) * (g(j/11) - g((j-1)/11)) for level j
    levels = list(range(1, 12))
    inc = {j: (max(Fraction(0), Fraction(j, 11) - d) - max(Fraction(0), Fraction(j - 1, 11) - d)) / (1 - d) for j in levels}
    den = 1
    for y, w in P7.items():
        for j in levels:
            v = w * inc[j]; den = den * v.denominator // np.gcd(den, v.denominator)
    SCALE = 10 ** 12            # objective coefficients rounded UP to integers at this scale -> upper bound
    model = cp_model.CpModel()
    z = {(m, r): model.NewBoolVar(f"z{m}_{r}") for m in DIV210 for r in range(m)}
    for m in DIV210:
        model.AddExactlyOne(z[(m, r)] for r in range(m))
    obj = []
    for y, w in P7.items():
        if w == 0: continue
        n_y = sum(z[(m, y % m)] for m in DIV210)
        ws = [model.NewBoolVar(f"w{y}_{j}") for j in levels]
        model.Add(sum(ws) <= n_y)
        for j in range(1, len(ws)): model.AddImplication(ws[j], ws[j - 1])
        for j, wv in zip(levels, ws):
            c = w * inc[j]
            ci = -((-c.numerator * SCALE) // c.denominator)       # ceiling
            if ci: obj.append(ci * wv)
    model.Maximize(sum(obj))
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = workers
    st = solver.Solve(model)
    status = solver.StatusName(st)
    best = Fraction(int(solver.ObjectiveValue()), SCALE) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    bound = Fraction(int(solver.BestObjectiveBound()), SCALE)
    a = {m: next(r for r in range(m) if solver.Value(z[(m, r)])) for m in DIV210} if best is not None else None
    # exact re-evaluation of the returned configuration (relaxed alpha = min(1, n/11))
    exact = None
    if a is not None:
        g = g_fn(D11); exact = Fraction(0)
        for y, w in P7.items():
            n = sum(1 for m in DIV210 if y % m == a[m])
            exact += w * g(min(Fraction(1), Fraction(n, 11)))
    return dict(status=status, objective=float(best) if best is not None else None, bound_upper=float(bound),
                exact_value_of_solution=(str(exact), float(exact)) if exact is not None else None, residues=a)


if __name__ == "__main__":
    D5, D7, D11 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    # nested stage-5 family and fully aligned stage-7 family (classes distinct), all at residue 0
    F5 = {2: (0, 0), 3: (0, 1), 6: (0, 2)}
    F7 = {m: (0, i) for i, m in enumerate(M7)}
    P7, l5, l7 = build_P7(F5, F7, D5, D7)
    t0 = time.time()
    b, t, vals = aligned_stage11(P7, D11)
    mod = model_stage11(D5, D7, D11); modt = model_stage11_truncated(D5, D7, D11)
    print(f"loss5={float(l5):.6f} loss7={float(l7):.6f}  aligned stage-11 best={float(b):.6f} at t={t} "
          f"(model {float(mod):.6f}, model truncated {float(modt):.6f})  [{time.time()-t0:.1f}s]")
    if len(sys.argv) > 4:
        res = cpsat_relaxed_stage11(P7, D11, time_limit=float(sys.argv[4]))
        print(json.dumps(res, indent=1, default=str))
