"""Integer certificate for the squarefree minimum-modulus bound via the convex stop-loss comparison.

Claim certified: every finite family of arithmetic progressions with DISTINCT SQUAREFREE moduli,
all > C0, fails to cover the integers; hence every covering system with distinct squarefree moduli
has a modulus <= C0.

Method (all inequalities listed in REPORT.md):
  * BBMST/CFT distorted measure, stage k = k-th prime p_k, distortion delta_k = D_k/L (rational).
  * loss_k = E(alpha_k - delta_k)_+/(1-delta_k) (exact, BBMST eq. (loss)), alpha_k = covered fibre
    proportion.  alpha_k(x) <= S_k(x) := (1/p_k) #{m : m p_k in the family} 1[x in cyl(m)] (union bound,
    distinct moduli => each m once); the comparison lemma (repository Lemma 'Product comparison' with
    gamma_j = 1, i.e. Bernoulli(nu_j/p_j) coordinates) gives E phi(S_k) <= E phi(S_k(D)) for convex
    nondecreasing phi, where D = product of independently present earlier primes, prime q present with
    probability nu_q/q, nu_q = 1/(1-delta_q), and S_k(D) = (1/p_k) #{m | D : m p_k > C0}.
  * Three valid stage bounds, the minimum is used:
      first : loss_k <= E alpha_k <= E S_k(D)                 (alpha <= 1, phi = id)
      stop  : loss_k <= E(S_k(D) - delta_k)_+/(1-delta_k)     (phi = (x-delta)_+)
      quad  : loss_k <= E S_k(D)^2/(4 delta_k (1-delta_k))    ((x-d)_+/(1-d) <= x^2/(4d(1-d)))
  * p_k <= C0: exact enumeration of all 2^(k-1) configurations of the earlier primes (Fractions).
  * C0 < p_k <= P_EXACT: S = 2^N/p_k, N = number of present earlier primes; Poisson-binomial law of N
    with OUTWARD rounded integer upper bounds (scale SC), truncated at NMAX with the explicit remainder
    E[2^N 1{N > NMAX}] <= (2/3)^NMAX E 3^N, E 3^N = prod (1 + 2 nu_q/q).
  * P_EXACT < p_k, k <= N_STAGES: closed forms E S = C/p, E S^2 = T/p^2 with C = prod(1+nu/q),
    T = prod(1+3nu/q) (integer upper bounds).
  * k > N_STAGES with delta = 1/2: CFT Lemma 4.7 (Rosser-Schoenfeld Thm 1 and Thm 5), re-derived in
    the report, evaluated with an interval logarithm:  sum_{k>N} loss_k <= 2 c2 M0 e^{6 c1} Q(log p_N)/p_N,
    c1 = -log log p_N + 1/log^2 p_N, c2 = 1 + 3/(2 log p_N), M0 = T_N,
    Q(l) = l^5 + 5 l^4 + 20 l^3 + 60 l^2 + 120 l + 120.
Everything is integer arithmetic with ceiling divisions in the direction that enlarges the bound.

Usage: python verify_squarefree.py schedule.json out.json [--p-exact P] [--n-stages N]
schedule.json: optimiser output with fields C0, edges (band upper stage indices), params (deltas as
floats); they are converted to D = round(delta * L), L = 10**6, and the D's are recorded.
Standard library only.
"""
import argparse, bisect, json, math, sys, time
from fractions import Fraction

L = 10 ** 6          # delta denominators
SC = 10 ** 60        # probability / loss scale
NMAX = 120           # truncation of the present-prime count


def ceildiv(a, b):
    assert b > 0
    return -((-a) // b)


def primes_to(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


def nth_primes(n):
    limit = max(100, int(n * (math.log(n) + math.log(math.log(n)))) + 100)
    ps = primes_to(limit)
    assert len(ps) >= n
    return ps


# ---------- interval logarithm (Fractions, explicit remainder) ----------
def atanh_series_bounds(z, terms):
    """lower/upper bounds for 2*atanh(z) = log((1+z)/(1-z)), 0 <= z < 1/2, Fractions."""
    assert 0 <= z < Fraction(1, 2)
    s = Fraction(0); zp = z
    for i in range(terms):
        s += zp / (2 * i + 1); zp *= z * z
    # remainder 2*sum_{i>=terms} z^(2i+1)/(2i+1) <= 2 z^(2 terms+1)/((2 terms+1)(1-z^2))
    rem = 2 * zp / ((2 * terms + 1) * (1 - z * z))
    return 2 * s, 2 * s + rem


def log_bounds(x, terms=40):
    """Lower/upper Fractions for log(x), x a positive Fraction or int."""
    x = Fraction(x); assert x > 0
    e = 0
    while x >= 2: x /= 2; e += 1
    while x < 1: x *= 2; e -= 1
    l2lo, l2hi = atanh_series_bounds(Fraction(1, 3), terms)      # log 2 = 2 atanh(1/3)
    ylo, yhi = atanh_series_bounds((x - 1) / (x + 1), terms)     # log y, y in [1,2)
    if e >= 0: return e * l2lo + ylo, e * l2hi + yhi
    return e * l2hi + ylo, e * l2lo + yhi


def exp_upper(x):
    """Upper Fraction for exp(x), 0 <= x < 1: 1 + x + x^2/2 + x^3/(6(1-x))."""
    x = Fraction(x); assert 0 <= x < 1
    return 1 + x + x * x / 2 + x ** 3 / (6 * (1 - x))


# ---------- stage computations ----------
def exact_small_stage(p, earlier, Dq, C0, D):
    """p <= C0.  Forbidden m are the divisors of D with m <= b = floor(C0/p); they are composed of the
    'small' primes q <= b.  Condition exactly (Fractions) on the configuration of the small primes; the
    count N of present 'large' primes (b < q < p) has a Poisson-binomial law, handled with outward
    rounded integer upper bounds and the truncation remainder.  Returns upper integer bounds (scale SC)
    for first = E S and stop = E(S-delta)_+/(1-delta); quad is not used here."""
    b = C0 // p
    small = [q for q in earlier if q <= b]; large = [q for q in earlier if q > b]
    pb = PoissonBinomialUpper(); E3_hi = SC
    for q in large:
        pb.add_prime(q, Dq[q]); den = (L - Dq[q]) * q; E3_hi = ceildiv(E3_hi * (den + 2 * L), den)
    prs = [Fraction(L, (L - Dq[q]) * q) for q in small]
    first = Fraction(0); stop = Fraction(0)
    for mask in range(1 << len(small)):
        prob = Fraction(1); divs = [1]; kk = 0
        for i, q in enumerate(small):
            if mask >> i & 1:
                prob *= prs[i]; divs = divs + [d * q for d in divs]; kk += 1
            else:
                prob *= 1 - prs[i]
        if prob == 0: continue
        cnt = sum(1 for d in divs if d <= b)       # forbidden divisors (m p <= C0)
        f_int = 0; s_int = 0
        for n in range(NMAX + 1):
            if not pb.P[n]: continue
            c = (1 << (kk + n)) - cnt; assert c >= 0
            f_int += pb.P[n] * c
            v = c * L - D * p
            if v > 0: s_int += pb.P[n] * v
        # remainder E[2^(kk+N) 1{N>NMAX}] <= 2^kk (2/3)^NMAX E 3^N   (scale SC)
        rem = ceildiv((1 << kk) * (2 ** NMAX) * E3_hi, 3 ** NMAX)
        first += prob * Fraction(ceildiv(f_int + rem, p))
        stop += prob * Fraction(ceildiv(s_int + rem * L, p * (L - D)))
    up = lambda f: ceildiv(f.numerator, f.denominator)
    return up(first), up(stop), None, float(first) / SC, None


class PoissonBinomialUpper:
    """Upper bounds (scale SC) for P(N = n), n <= NMAX, N = number of present primes."""
    def __init__(self):
        self.P = [SC] + [0] * NMAX

    def add_prime(self, q, D):
        den = (L - D) * q; stay = den - L           # 1 - pr = stay/den, pr = L/den
        assert stay >= 0 and den > 0
        P = self.P; new = [0] * (NMAX + 1)
        new[0] = ceildiv(P[0] * stay, den)
        for n in range(1, NMAX + 1):
            v = P[n] * stay + P[n - 1] * L
            if v: new[n] = ceildiv(v, den)
        self.P = new

    def stop_upper(self, p, D):
        """ceil of sum_n P_hi[n] (2^n/p - D/L)_+/(1-D/L), scale SC."""
        tot = 0
        for n in range(NMAX + 1):
            v = (1 << n) * L - D * p
            if v > 0 and self.P[n]: tot += self.P[n] * v
        return ceildiv(tot, p * (L - D))


def verify(sched_path, out_path, p_exact, n_stages, checkpoint_every=5000):
    t0 = time.time()
    w = json.load(open(sched_path))
    C0 = int(w["C0"]); edges = [int(e) for e in w["edges"]]; params = w["params"]
    assert len(edges) == len(params) and edges == sorted(edges) and edges[-1] >= n_stages
    assert 1 <= C0 < 10 ** 6 and p_exact >= C0
    Ds = []
    for x in params:
        D = int(round(float(x) * L)); assert 0 <= D <= L // 2; Ds.append(D)

    def D_of_k(k):          # 1-based stage index; band i covers edges[i-1] < k <= edges[i]
        return Ds[bisect.bisect_left(edges, k)]
    ps = nth_primes(n_stages)
    assert ps[0] == 2 and ps[n_stages - 1] >= p_exact
    # running integer upper bounds, scale SC
    C_hi = SC; T_hi = SC; E3_hi = SC          # prod(1+pr), prod(1+3pr), prod(1+2pr)
    pb = PoissonBinomialUpper()
    Dq = {}
    total = 0                                 # integer upper bound of the sum of losses, scale SC
    sums = {"p_le_C0": 0, "exact": 0, "far": 0}
    wins = {"first": 0, "stop": 0, "quad": 0}
    range_edges = [10, 30, 100, 300, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]
    range_sums = [0] * (len(range_edges) + 1)
    rows = []; two_third = (2 ** NMAX, 3 ** NMAX)
    k = 0; checkpoints = []
    for idx in range(n_stages):
        p = ps[idx]; k = idx + 1; D = D_of_k(k)
        if p <= C0:
            earlier = ps[:idx]
            first, stop, quad, ES, ES2 = exact_small_stage(p, earlier, Dq, C0, D)
            cands = [(first, "first"), (stop, "stop")] + ([(quad, "quad")] if quad is not None else [])
            loss, kind = min(cands); sums["p_le_C0"] += loss
            rows.append(dict(k=k, p=p, D=D, first=first / SC, stop=stop / SC, quad=(quad / SC if quad else None), loss=loss / SC, kind=kind))
        elif p <= p_exact:
            first = ceildiv(C_hi, p)
            # truncation remainder: (2/3)^NMAX * E3 / (p (1-delta))
            rem = ceildiv(two_third[0] * E3_hi * L, two_third[1] * p * (L - D))
            stop = pb.stop_upper(p, D) + rem
            cands = [(first, "first"), (stop, "stop")]
            if D: cands.append((ceildiv(T_hi * L * L, 4 * D * (L - D) * p * p), "quad"))
            loss, kind = min(cands); sums["exact"] += loss
            if k <= 60 or k % 1000 == 0:
                rows.append(dict(k=k, p=p, D=D, first=first / SC, stop=stop / SC, loss=loss / SC, kind=kind))
        else:
            first = ceildiv(C_hi, p)
            cands = [(first, "first")]
            if D: cands.append((ceildiv(T_hi * L * L, 4 * D * (L - D) * p * p), "quad"))
            loss, kind = min(cands); sums["far"] += loss
        wins[kind] += 1; total += loss
        range_sums[bisect.bisect_left(range_edges, p)] += loss
        # update products and the distribution with this stage's prime
        den = (L - D) * p
        C_hi = ceildiv(C_hi * (den + L), den)
        T_hi = ceildiv(T_hi * (den + 3 * L), den)
        E3_hi = ceildiv(E3_hi * (den + 2 * L), den)
        if p <= p_exact: pb.add_prime(p, D)
        Dq[p] = D
        if k % checkpoint_every == 0 or k == n_stages:
            checkpoints.append(dict(k=k, p=p, partial_total=total / SC, seconds=round(time.time() - t0, 1)))
            print(f"stage {k} p={p} partial total <= {total / SC:.9f}  ({time.time() - t0:.0f}s)", flush=True)
    # ---- tail beyond N_STAGES: CFT Lemma 4.7 with delta = 1/2 ----
    N = n_stages; pN = ps[N - 1]
    assert N >= 61 and pN >= 293          # Rosser-Schoenfeld lower bound needs x >= 286
    llo, lhi = log_bounds(pN)
    assert llo > 3                        # log^6 t / t^2 decreasing beyond e^3
    M0_hi = Fraction(T_hi, SC)
    c2_hi = 1 + Fraction(3, 2) / llo
    e6c1_hi = exp_upper(6 / (llo * llo)) / llo ** 6
    Q_hi = lhi ** 5 + 5 * lhi ** 4 + 20 * lhi ** 3 + 60 * lhi ** 2 + 120 * lhi + 120
    tail_frac = 2 * c2_hi * M0_hi * e6c1_hi * Q_hi / pN
    tail = ceildiv(tail_frac.numerator * SC, tail_frac.denominator)
    grand = total + tail
    status = "CERTIFIED" if grand < SC else "NOT_CERTIFIED"
    result = dict(
        status=status, C0=C0, claim=(f"every covering system with distinct squarefree moduli has a modulus <= {C0}"
                                     if status == "CERTIFIED" else "total loss bound not below 1"),
        total_loss_upper=f"{grand / SC:.12f}", margin=f"{(SC - grand) / SC:.12f}",
        total_loss_upper_integer=str(grand), scale=str(SC),
        parts=dict(p_le_C0=f"{sums['p_le_C0'] / SC:.12f}", exact=f"{sums['exact'] / SC:.12f}",
                   far=f"{sums['far'] / SC:.12f}", tail=f"{tail / SC:.12f}"),
        parameters=dict(L=L, NMAX=NMAX, p_exact=p_exact, n_stages=n_stages, p_N=pN,
                        delta_bands=[dict(k_hi=e, D=D) for e, D in zip(edges, Ds)],
                        delta_beyond_N="1/2", tail_method="CFT Lemma 4.7 (Rosser-Schoenfeld Thm 1, Thm 5)"),
        tail_details=dict(log_pN=[str(llo), str(lhi)], M0_upper=f"{float(M0_hi):.6f}", c2_upper=f"{float(c2_hi):.9f}",
                          exp6c1_upper=f"{float(e6c1_hi):.6e}"),
        loss_by_prime_range={(f"p<{range_edges[0]}" if i == 0 else f"{range_edges[i-1]}<=p<{range_edges[i]}" if i < len(range_edges) else f"p>={range_edges[-1]}"): f"{v / SC:.9f}" for i, v in enumerate(range_sums)},
        wins=wins, checkpoints=checkpoints, rows=rows, seconds=round(time.time() - t0, 1),
        schedule_file=sched_path)
    json.dump(result, open(out_path, "w"), indent=1)
    print(json.dumps({k_: v for k_, v in result.items() if k_ not in ("rows", "checkpoints")}, indent=1))
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("schedule"); ap.add_argument("out")
    ap.add_argument("--p-exact", type=int, default=100000)
    ap.add_argument("--n-stages", type=int, default=10 ** 6)
    a = ap.parse_args()
    if not __debug__: raise RuntimeError("assertions must be enabled")
    verify(a.schedule, a.out, a.p_exact, a.n_stages)
