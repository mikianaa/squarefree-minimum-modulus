"""Independent verifier for the squarefree minimum-modulus certificate (agent sqaudit).

Written from the sqfree REPORT.md description only (Sections 3.2/3.3), without reading
verify_squarefree.py.  Standard library only.  All stage bounds are computed as integer
upper bounds at scale SC = 10^60 with ceiling divisions (outward rounding); the finitely
many stages p <= C0 and the tail constants use exact Fractions.

Model (independent comparison model): for stage k (prime p = p_k) the earlier primes q = p_j,
j < k, are present independently with probability pr_q = nu_q / q, nu_q = 1/(1 - delta_q),
D = product of present primes, S = (1/p) #{m | D : m p > C0}.  Stage bound =
min( E S , E (S - delta)_+ / (1 - delta) , E S^2 / (4 delta (1 - delta)) ).

Usage: python -I audit_verify.py opt10.json out.json [--p-exact 1000000] [--n-stages 1000000]
"""
import bisect, json, sys, time
from fractions import Fraction
from itertools import combinations

SC = 10 ** 60          # scale of integer upper bounds
L = 10 ** 6            # delta = D / L
NMAX = 120             # Poisson-binomial truncation


def cdiv(a, b):
    """ceil(a / b) for integers a >= 0, b > 0."""
    return -((-a) // b)


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


def nth_prime_list(n):
    # upper bound for p_n: n (log n + log log n) for n >= 6
    import math
    bound = int(n * (math.log(n) + math.log(math.log(n)))) + 10
    ps = sieve(bound)
    assert len(ps) >= n
    return ps[:n]


def load_schedule(path):
    sch = json.load(open(path))
    C0 = int(sch["C0"])
    edges = [int(e) for e in sch["edges"]]
    params = sch["params"]
    assert len(edges) == len(params)
    assert edges and edges[0] >= 1
    assert all(a < b for a, b in zip(edges, edges[1:]))
    Ds = [int(round(float(d) * L)) for d in params]
    for D in Ds:
        assert 0 <= D <= L // 2, D
    return C0, edges, Ds


def D_for_stage(k, edges, Ds):
    i = bisect.bisect_left(edges, k)
    if i == len(edges):
        raise ValueError("stage %d beyond schedule" % k)
    return Ds[i]


# ----------------------------------------------------------------------------------------
# log and exp enclosures (exact Fractions)

def atanh_enclosure(z, terms=40):
    """[lo, hi] for atanh(z), 0 <= z < 1 (Fraction)."""
    assert 0 <= z < 1
    s = Fraction(0)
    zp = z
    z2 = z * z
    for i in range(terms):
        s += zp / (2 * i + 1)
        zp *= z2
    # remainder sum_{i>=terms} z^(2i+1)/(2i+1) <= z^(2 terms+1) / ((2 terms + 1)(1 - z^2))
    rem = zp / ((2 * terms + 1) * (1 - z2))
    return s, s + rem


def log_enclosure(y):
    """[lo, hi] for log(y), y > 0 Fraction.  Uses log y = a log 2 + log(y / 2^a)."""
    y = Fraction(y)
    a = 0
    while y >= 2:
        y /= 2
        a += 1
    while y < 1:
        y *= 2
        a -= 1
    # now 1 <= y < 2
    l2lo, l2hi = atanh_enclosure(Fraction(1, 3))
    l2lo, l2hi = 2 * l2lo, 2 * l2hi
    z = (y - 1) / (y + 1)
    lo, hi = atanh_enclosure(z)
    lo, hi = 2 * lo, 2 * hi
    if a >= 0:
        return a * l2lo + lo, a * l2hi + hi
    return a * l2hi + lo, a * l2lo + hi


def exp_upper(x):
    """upper bound for e^x, 0 <= x < 1 Fraction."""
    assert 0 <= x < 1
    return 1 + x + x * x / 2 + x ** 3 / (6 * (1 - x))


# ----------------------------------------------------------------------------------------

def main():
    args = sys.argv[1:]
    sched_path, out_path = args[0], args[1]
    p_exact = 10 ** 6
    n_stages = 10 ** 6
    if "--p-exact" in args:
        p_exact = int(args[args.index("--p-exact") + 1])
    if "--n-stages" in args:
        n_stages = int(args[args.index("--n-stages") + 1])
    t0 = time.time()
    C0, edges, Ds = load_schedule(sched_path)
    assert edges[-1] >= n_stages
    primes = nth_prime_list(n_stages)
    pN = primes[n_stages - 1]
    print("C0 =", C0, " n_stages =", n_stages, " p_N =", pN, " primes sieved in %.1fs" % (time.time() - t0))

    # running products (integer upper bounds at scale SC): C = prod(1+pr), T = prod(1+3pr), E3 = prod(1+2pr)
    C = SC
    T = SC
    E3 = SC
    # Poisson-binomial upper bounds for N = number of present primes among the "large" earlier
    # primes.  For stages p > C0 all earlier primes are large.  For stages p <= C0 the small
    # primes (q <= C0 // p) are enumerated exactly and only the large ones go through the
    # recursion; since q <= C0//p < p and the recursion is built incrementally over ALL earlier
    # primes, we handle the stages p <= C0 by a separate exact computation (all earlier primes
    # enumerated exactly when few, else recursion restricted to large primes).
    PN = [SC] + [0] * NMAX     # PN[n] >= SC * P(N = n) over earlier primes processed so far

    rows = []
    total = 0                  # integer at scale SC
    part_small = 0
    part_exact = 0
    part_far = 0
    wins = {"first": 0, "stop": 0, "quad": 0}
    by_range = {}
    range_edges = [10, 30, 100, 300, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]

    def range_label(p):
        lo = 1
        for e in range_edges:
            if p < e:
                return "%d<=p<%d" % (lo, e)
            lo = e
        return "p>=%d" % lo

    prs = []  # exact Fractions pr_q for earlier primes (only needed for p <= C0 stages)
    small_stage_frac = Fraction(0)

    for k in range(1, n_stages + 1):
        p = primes[k - 1]
        D = D_for_stage(k, edges, Ds)
        delta = Fraction(D, L)
        nu_num, nu_den = L, L - D           # nu = L / (L - D)
        if p <= C0:
            # exact stage: enumerate the small primes (q <= C0 // p) exactly; large earlier
            # primes (C0//p < q < p) via an exact Poisson-binomial in Fractions (few primes).
            earlier = primes[: k - 1]
            small = [q for q in earlier if q <= C0 // p]
            large = [q for q in earlier if q > C0 // p]
            pr = {q: prs[i] for i, q in enumerate(earlier)}
            # distribution of N over large primes (exact)
            distN = {0: Fraction(1)}
            for q in large:
                nd = {}
                for n, w in distN.items():
                    nd[n] = nd.get(n, 0) + w * (1 - pr[q])
                    nd[n + 1] = nd.get(n + 1, 0) + w * pr[q]
                distN = nd
            first = Fraction(0)
            stop = Fraction(0)
            ES2 = Fraction(0)
            for r in range(len(small) + 1):
                for Sm in combinations(small, r):
                    w = Fraction(1)
                    for q in small:
                        w *= pr[q] if q in Sm else (1 - pr[q])
                    # forbidden divisors of prod(Sm): m with m p <= C0
                    cnt = 0
                    for r2 in range(len(Sm) + 1):
                        for sub in combinations(Sm, r2):
                            m = 1
                            for q in sub:
                                m *= q
                            if m * p <= C0:
                                cnt += 1
                    for n, wn in distN.items():
                        S = Fraction(2 ** (len(Sm) + n) - cnt, p)
                        first += w * wn * S
                        ES2 += w * wn * S * S
                        if S > delta:
                            stop += w * wn * (S - delta) / (1 - delta)
            cands = {"first": first, "stop": stop}
            if D > 0:
                cands["quad"] = ES2 / (4 * delta * (1 - delta))
            kind = min(cands, key=lambda s: cands[s])
            loss = cands[kind]
            loss_int = cdiv(loss.numerator * SC, loss.denominator)
            part_small += loss_int
            rows.append({"k": k, "p": p, "D": D, "first": float(first), "stop": float(stop),
                         "quad": (float(cands["quad"]) if "quad" in cands else None),
                         "loss": float(loss), "kind": kind, "exact": True})
        elif p <= p_exact:
            # stop-loss via truncated Poisson-binomial upper bounds + remainder
            num = 0
            for n in range(NMAX + 1):
                if PN[n] == 0:
                    continue
                t = (2 ** n) * L - D * p
                if t > 0:
                    num += PN[n] * t
            stop_int = cdiv(num, p * (L - D))
            # remainder E[2^N 1{N > NMAX}] <= (2/3)^NMAX E 3^N ; contributes /(p (1-delta))
            rem_int = cdiv(E3 * (2 ** NMAX) * L, (3 ** NMAX) * p * (L - D))
            stop_int += rem_int
            first_int = cdiv(C, p)
            cands = {"first": first_int, "stop": stop_int}
            if D > 0:
                cands["quad"] = cdiv(T * L * L, 4 * D * (L - D) * p * p)
            kind = min(cands, key=lambda s: cands[s])
            loss_int = cands[kind]
            part_exact += loss_int
            if k <= 60 or p in (23, 101) or k % 5000 == 0:
                rows.append({"k": k, "p": p, "D": D, "first": first_int / SC, "stop": stop_int / SC,
                             "quad": (cands["quad"] / SC if "quad" in cands else None),
                             "loss": loss_int / SC, "kind": kind})
        else:
            first_int = cdiv(C, p)
            cands = {"first": first_int}
            if D > 0:
                cands["quad"] = cdiv(T * L * L, 4 * D * (L - D) * p * p)
            kind = min(cands, key=lambda s: cands[s])
            loss_int = cands[kind]
            part_far += loss_int
        wins[kind] += 1
        total += loss_int
        lab = range_label(p)
        by_range[lab] = by_range.get(lab, 0) + loss_int

        # update the running products and the Poisson-binomial with this stage's prime
        # pr = nu / p = L / ((L - D) p)
        pr_num, pr_den = L, (L - D) * p
        prs.append(Fraction(pr_num, pr_den))
        C = cdiv(C * (pr_den + pr_num), pr_den)
        T = cdiv(T * (pr_den + 3 * pr_num), pr_den)
        if p <= p_exact:   # the Poisson-binomial and E 3^N are only needed for exact stages
            E3 = cdiv(E3 * (pr_den + 2 * pr_num), pr_den)
            newPN = [0] * (NMAX + 1)
            for n in range(NMAX + 1):
                v = PN[n] * (pr_den - pr_num)
                if n >= 1:
                    v += PN[n - 1] * pr_num
                newPN[n] = cdiv(v, pr_den)
            PN = newPN
        if k % 100000 == 0:
            print("  k=%d p=%d partial=%.12f  %.1fs" % (k, p, total / SC, time.time() - t0))

    # ---- tail: k > n_stages, delta = 1/2, nu = 2 (I8)-(I11)
    M0 = Fraction(T, SC)                      # upper bound for T_N = prod_{j<=N}(1+3 pr_j)
    llo, lhi = log_enclosure(Fraction(pN))
    assert llo > 3
    c2 = 1 + Fraction(3, 2) / llo             # upper bound (decreasing in log)
    x = 6 / (llo * llo)                       # upper bound of 6/log^2
    e6c1 = exp_upper(x) / llo ** 6            # e^{6 c1} = e^{6/log^2} / log^6
    l = lhi
    Q = l ** 5 + 5 * l ** 4 + 20 * l ** 3 + 60 * l ** 2 + 120 * l + 120
    tail = 2 * c2 * M0 * e6c1 * Q / pN
    tail_int = cdiv(tail.numerator * SC, tail.denominator)
    total += tail_int
    by_range["tail(k>N)"] = tail_int

    out = {
        "status": "CERTIFIED" if total < SC else "NOT CERTIFIED",
        "C0": C0,
        "total_loss_upper": total / SC,
        "total_loss_upper_integer": str(total),
        "scale": str(SC),
        "margin": (SC - total) / SC,
        "parts": {"p_le_C0": part_small / SC, "exact_stop": part_exact / SC,
                  "far_closed_form": part_far / SC, "tail": tail_int / SC},
        "parameters": {"L": L, "NMAX": NMAX, "p_exact": p_exact, "n_stages": n_stages, "p_N": pN,
                       "bands": [{"k_hi": e, "D": D} for e, D in zip(edges, Ds)]},
        "tail_details": {"log_pN_lo": float(llo), "log_pN_hi": float(lhi), "M0_upper": float(M0),
                         "c2_upper": float(c2), "exp6c1_upper": float(e6c1), "Q_upper": float(Q)},
        "wins": wins,
        "loss_by_range": {kk: v / SC for kk, v in by_range.items()},
        "rows": rows,
        "seconds": round(time.time() - t0, 1),
        "schedule_file": sched_path,
    }
    json.dump(out, open(out_path, "w"), indent=1)
    print("status:", out["status"], " total = %.12f  margin = %.12f" % (total / SC, (SC - total) / SC))
    print("parts:", out["parts"])
    print("wins:", wins)
    print("tail:", out["tail_details"])


if __name__ == "__main__":
    main()
