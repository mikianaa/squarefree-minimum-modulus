# Agent "sqfree": certified squarefree minimum-modulus bound via the convex stop-loss comparison

Date: 23 Oct 2026. Time spent: about 2.5 hours. Nothing in the repository was edited.

## 1. Idea tested

The transfer agent found numerically that the repository's convex stop-loss comparison (Lemmas
"Nested-event" and "Product comparison" of src/conditional-stoploss-improvement.tex), applied to the
squarefree problem of Cummings-Filaseta-Trifonov (CFT, arXiv:2211.08548, minimum modulus <= 118),
gives minimum modulus <= 14. Task: turn this into a certified result with an integer-arithmetic
verifier (outward rounding everywhere, standard library only) including an explicit infinite tail,
and push the cutoff down as far as the margin allows.

## 2. What was computed and read

Read: CFT 2211.08548 (Sections 2-4, Lemma 4.7, proof of Theorem 1.1, delta table); BBMST squarefree
paper 1901.11465 (introduction); the repository's src/conditional-stoploss-improvement.tex Sections 1-3
and the fractional-moment subsection; code/verify_sparse_conditioning.py (outward-rounded integer style,
L-scaled deltas, ceiling divisions); the transfer agent's REPORT.md and squarefree_convex.py.
Web: Dalton-Trifonov, Extreme covering systems, J. Integer Seq. 25 (2022) Art. 22.9.1 (Problem 4),
https://cs.uwaterloo.ca/journals/JIS/VOL25/Trifonov/trif3.pdf ; Klein, On a conjecture of Krukenberg
and a problem of Dalton and Trifonov, https://arxiv.org/abs/2508.18062 (general moduli, not squarefree);
CFT https://arxiv.org/abs/2211.08548 ; BBMST squarefree https://arxiv.org/abs/1901.11465 ;
repository notes/codex-readme-log-2026-10-14-15.md (squarefree targets 30 and 3) and literature/SOURCES.md.

Scripts (all in C:\Users\mekak\AppData\Local\Temp\claude\c--Users-mekak-code-minimum-modulus\23ed4d4e-d4c3-4fc1-93e5-481b63d3d735\scratchpad\agents\sqfree\ ; venv python for the optimiser, plain python 3.12 for the verifier):

* sqfree_float_opt.py - floating-point model and coordinate-descent optimiser of the delta schedule
  (62 bands: one per stage for k <= 30, then growing bands up to stage 10^6; a finer variant with
  per-stage bands to k = 80). Not a proof input; it only produces the schedule.
* verify_squarefree.py - the certificate. Integers only (fractions for the finitely many stages
  p <= C0 and for the tail constants). Reads a schedule JSON, converts deltas to D/L with L = 10^6,
  records the D's, and writes a JSON with every parameter, per-stage rows, per-range sums and the
  certified total. Run: python verify_squarefree.py optC0.json certC0.json --p-exact 1000000
* summarise.py - markdown tables from certificate JSONs.
* cft_sched.json / cft118_check.json - CFT's own delta table run through the verifier (validation).
* opt<C0>.json / .log - optimiser output; cert<C0>.json / .log - certificates (the proof objects).

## 3. Results

### 3.1 Literature sanity check (lower-bound side)

* Covering systems with distinct squarefree moduli exist with minimum modulus 2 (Krukenberg 1971:
  the only such systems whose lcm has <= 4 prime factors have lcm 210; the repository's Codex log
  verified a 14-class squarefree system modulo 210). So the true maximum of the minimum modulus is >= 2.
* Whether a covering system with distinct squarefree moduli can avoid the modulus 2 (minimum modulus
  >= 3) is OPEN: Dalton-Trifonov (2022) Problem 4. BBMST (1901.11465) prove that distinct ODD
  squarefree moduli never cover; Codex's notes give finite obstructions for minimum 3 with primes <= 11
  and for minimum 30 with primes <= 37. Klein's 2025 paper (2508.18062) settles different
  Dalton-Trifonov problems (general moduli: largest modulus / lcm when the minimum is 5 or 6).
* Best published upper bound: 118 (CFT). Any certified bound in [2, 118] is consistent with the
  literature, and every step below 118 is new.
* Squarefree bookkeeping: "all moduli > C0" and "all moduli > C0'" are the same hypothesis when no
  squarefree number lies in (C0', C0]. So C0 = 12 certifies "some modulus <= 11", C0 = 8 or 9 would
  certify "<= 7", C0 = 3 or 4 would certify "<= 3", C0 = 2 would answer Problem 4.

### 3.2 The comparison lemma specialised to {0,1} coordinates (the argument)

Setting (CFT Section 3 = BBMST fibre construction written as weights): coordinates X_j in Z_{p_j},
common period prod p_j, so every exponent gamma_j = 1. A progression with squarefree modulus m p_k,
m = prod_{j in J} p_j, J subset of {1..k-1}, is the cylinder fixing one residue in each coordinate
j in J and one residue y in coordinate k.

1. Loss formula: P_k(B_k) = E_{k-1}(alpha_k - delta_k)_+/(1-delta_k) exactly (CFT weights: on a fibre
   with alpha > delta the covered points carry total weight (alpha-delta) w_{k-1}(x)/(1-delta), else
   0). This is eq. (loss) of the tex with gamma = 1.
2. Cap: CFT (3) is w_k(x,y) <= w_{k-1}(x)/((1-delta_k) p_k), i.e. P(X_k = y | X_1..X_{k-1}) <= nu_k/p_k;
   sum_y w_k(x,y) = w_{k-1}(x) preserves every earlier prefix law, so the cap holds at all later
   stages: eq. (cap) with gamma_j = 1.
3. Union bound: alpha_k(x) = |F_x cap B_k|/p_k <= (1/p_k) sum_{A in A_k} 1[x in cyl_{J(A)}], because a
   hyperplane with fixed set J cup {k} meets the fibre F_x in exactly one point when x lies in its
   earlier-coordinate cylinder and in none otherwise. Distinct moduli = pairwise non-parallel: each J
   occurs at most once, so every coefficient is exactly 1/p_k. All J with m_J p_k <= C0 are absent.
4. Lemma "Product comparison" with gamma_j = 1: V_j in {0,1} with P(V_j >= 1) = nu_j/p_j <= 1 (as
   nu_j <= 2 <= p_j), i.e. independent Bernoulli(nu_j/p_j); prod_j 1[V_j >= a_{lj}] = 1[J_l subset
   {j : V_j = 1}] = 1[m_J | D], D = prod_{V_j = 1} p_j. The lemma's proof (condition on the preceding
   coordinates, apply the nested-event lemma to the last coordinate, induct) uses only (cap) and never
   the value of gamma, so for every convex nondecreasing phi on [0, infinity):
      E phi( sum_{A in A_k} (1/p_k) 1[cyl_{J(A)}] ) <= E phi( sum_{A in A_k} (1/p_k) 1[m_{J(A)} | D] ).
5. Adding the absent allowed J's adds nonnegative terms inside a nondecreasing phi, so the right side
   is at most E phi(S_k(D)), S_k(D) = (1/p_k) #{m | D : m p_k > C0}. With phi(x) = (x-delta)_+ and
   steps 1, 3: P_k(B_k) <= E(S_k(D) - delta_k)_+/(1-delta_k). With phi = x^2 and
   (x-d)_+/(1-d) <= x^2/(4d(1-d)) (pointwise: x^2 - 4d(x-d) = (x-2d)^2 >= 0) the quadratic bound; with
   alpha <= 1, (alpha-d)_+/(1-d) <= alpha and phi = x the first-moment bound.

### 3.3 Every inequality used by the certificate

Stage k, prime p = p_k, delta = D/L (L = 10^6, 0 <= D <= L/2), pr_q = nu_q/q = L/((L-D_q) q).

(I1)  P_k(B_k) = E(alpha_k - delta)_+/(1-delta)                        [exact, BBMST/CFT]
(I2)  P_k(B_k) <= min( E S, E(S-delta)_+/(1-delta), E S^2/(4 delta(1-delta)) ), S = S_k(D)   [3.2]
(I3)  For p > C0: S = 2^N/p, N = #present primes q < p; E S = C_{k-1}/p, E S^2 = T_{k-1}/p^2,
      C = prod(1+pr_q), T = prod(1+3pr_q) (because (1+B)^2 = 1+3B for B in {0,1}).
(I4)  Poisson-binomial recursion P(N'=n) = P(N=n)(1-pr) + P(N=n-1) pr; integer upper bounds at
      scale 10^60 with ceiling division; mass moving above NMAX = 120 is dropped (the upper bounds of
      P(N=n), n <= NMAX, remain valid; all stage functionals are nonnegative).
(I5)  Truncation remainder: E[2^N 1{N > NMAX}] <= E[2^N (3/2)^{N-NMAX}] = (2/3)^NMAX E 3^N,
      E 3^N = prod(1+2pr_q); added (times 1/(p(1-delta))) to the stop-loss bound.
(I6)  For p <= C0: forbidden m <= floor(C0/p) consist of the primes q <= floor(C0/p); the
      configuration of these primes is enumerated exactly (Fractions), the count N of the remaining
      earlier primes by (I4)-(I5); S = (2^{kk+N} - cnt)/p, cnt = #forbidden divisors of the present
      small primes. Only first and stop are used at these stages.
(I7)  Stages k <= N = 10^6 with p > P_EXACT = 10^6: min(C/p, T/(4 delta(1-delta) p^2)) by (I3).
(I8)  Tail k > N, delta_k = 1/2 (nu = 2): P_k(B_k) <= T_{k-1}/p_k^2 and
      T_{k-1} = M0 prod_{N<j<k}(1+6/p_j) <= M0 exp(6 sum_{p_N < q <= p_{k-1}} 1/q), M0 = T_N.
(I9)  Rosser-Schoenfeld (Illinois J. Math. 6 (1962) 64-94) Theorem 5, (3.17)-(3.18):
      sum_{q <= x} 1/q < log log x + B + 1/(2 log^2 x) for x > 1 and > log log x + B - 1/(2 log^2 x)
      for x >= 286. Hence (p_N >= 286) sum_{p_N<q<=p_{k-1}} 1/q <= log log p_{k-1} - log log p_N +
      1/log^2 p_N, so P_k(B_k) <= M0 e^{6c1} log^6 p_k / p_k^2, c1 = -log log p_N + 1/log^2 p_N.
(I10) f(t) = log^6 t / t^2 decreases for t > e^3; sum_{k>N} f(p_k) = int_{p_N}^inf (-f'(t))
      (pi(t) - N) dt <= int_{p_N}^inf pi(t) 2 log^6 t / t^3 dt. Rosser-Schoenfeld Theorem 1 (3.6):
      pi(x) < (x/log x)(1 + 3/(2 log x)) for x > 1, so pi(t) <= c2 t/log t for t >= p_N with
      c2 = 1 + 3/(2 log p_N). Therefore sum_{k>N} f(p_k) <= (2 c2/p_N) Q(log p_N),
      Q(l) = l^5 + 5l^4 + 20l^3 + 60l^2 + 120l + 120 (antiderivative of log^5 t / t^2).
      Total tail <= 2 c2 M0 e^{6 c1} Q(log p_N)/p_N. This is CFT's Lemma 4.7, re-derived.
(I11) log p_N enclosed via log 2 = 2 atanh(1/3) and log y = 2 atanh((y-1)/(y+1)), 40 terms, remainder
      2 z^{2n+1}/((2n+1)(1-z^2)); e^x <= 1 + x + x^2/2 + x^3/(6(1-x)) for 0 <= x < 1; Q with the
      upper log, 1/log^6, c2 and e^{6/log^2} with the lower log.
(I12) Conclusion: sum_{k>=1} (stage bounds) < 1 => P_n(union B_k) <= sum_k P_k(B_k) < 1 (prefix laws
      preserved) => an uncovered point of Q_n, i.e. an uncovered residue class modulo prod_{j<=n} p_j
      (CFT Lemma 3.1 / BBMST). The bounds are uniform in the family, so stages beyond the family's
      largest prime (which contribute 0) are harmless.

### 3.4 Validation of the verifier

* CFT's delta table at C0 = 118 (cft_sched.json): stages 1-7 (delta = 0, first moment) give exactly
  194/1365 = 0.142124542124...; the tail at N = 10^6 is 0.000402960697 against CFT's (11)
  0.000402960685 (mine is an outward-rounded upper bound; agreement to 10 digits confirms M0, c1, c2,
  Q and the interval logarithm). Stages 8..10^6 come out at 0.64 instead of CFT's 0.857 because the
  verifier may use the stop-loss there; total 0.7831 < 1 at C0 = 118 with CFT's own deltas.
* Stage-by-stage agreement with the floating-point model to 9 digits at every printed stage (C0 = 14).
* Closed-form / Monte Carlo checks of the exact-distribution routine were done by the transfer agent.

### 3.5 Certified results (cert14.json, cert12.json, cert10.json; ~20 s each)

| C0 (all moduli > C0) | bound | status | certified total | margin 1 - total | p<=C0 | C0<p<=1e6 exact | far (quad) | tail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 14 | <= 14 | CERTIFIED | 0.823434461 | 0.176565539 | 0.269689 | 0.551330 | 0.001828 | 0.000587 |
| 12 | <= 11 | CERTIFIED | 0.885626440 | 0.114373560 | 0.253381 | 0.629780 | 0.001865 | 0.000600 |
| 10 | <= 10 | CERTIFIED | 0.924567867 | 0.075432133 | 0.211370 | 0.710628 | 0.001945 | 0.000624 |
| 10 | <= 10 | CERTIFIED | 0.922645169 | 0.077354831 | 0.211370 | 0.710651 | 0.000000 | 0.000624 |

(cert10full.json = the C0 = 10 schedule with the exact range extended to all 10^6 stages, 120 s; the far part vanishes and the total drops by 0.0019.)

Loss by prime range (certified upper bounds):

| prime range | C0=14 | C0=12 | C0=10 | C0=10 |
| --- | --- | --- | --- | --- |
| p<10 | 0.172387 | 0.211370 | 0.211370 | 0.211370 |
| 10<=p<30 | 0.273546 | 0.285059 | 0.303481 | 0.303481 |
| 30<=p<100 | 0.181694 | 0.186739 | 0.196443 | 0.196443 |
| 100<=p<300 | 0.107294 | 0.110896 | 0.116661 | 0.116661 |
| 300<=p<1000 | 0.050499 | 0.052358 | 0.055283 | 0.055283 |
| 1000<=p<10000 | 0.030286 | 0.031264 | 0.032961 | 0.032961 |
| 10000<=p<100000 | 0.004977 | 0.005128 | 0.005430 | 0.005430 |
| 100000<=p<1000000 | 0.000335 | 0.000347 | 0.000369 | 0.000369 |
| 1000000<=p<10000000 | 0.001721 | 0.001756 | 0.001832 | 0.000022 |
| 10000000<=p<100000000 | 0.000107 | 0.000109 | 0.000113 | 0.000000 |

Where the loss sits (C0 = 14): 0.17 for p < 10 (p = 5 first moment 0.1000, p = 7 stop-loss 0.0724),
0.27 for 10 <= p < 30, 0.18 for 30 <= p < 100, 0.11 for 100-300, 0.05 for 300-1000, 0.035 beyond
1000 including the far part and the tail. Stage bound used: first moment at p = 2, 3, 5 (delta = 0),
stop-loss at every other stage up to p = 10^6, quadratic bound on the 921,502 stages between
p = 10^6 and stage 10^6 (sum 0.0018), tail 0.0006.

Delta schedule of the C0 = 14 certificate (k_hi: delta; a band is the stages up to k_hi):
1:0, 2:0, 3:0, 4:0.1425, 5:0.27, 6:0.2325, 7:0.2375, 8:0.215, 9:0.35, 10:0.275, 11:0.2625, 12:0.2175, 13:0.39, 14:0.3675, 15:0.3375, 16:0.3, 17:0.27, 18:0.2625, 19:0.24, 20:0.2325, 21:0.435, 22:0.405, 23:0.3825, 24:0.36, 25:0.33, 26:0.3175, 27:0.31, 28:0.3025, 29:0.295, 30:0.2875, 35:0.2575, 40:0.385, 50:0.3325, 60:0.2725, 75:0.3775, 100:0.325, 130:0.4, 170:0.335, 230:0.4025, 300:0.3425, 400:0.4175, 500:0.35, 700:0.41, 1000:0.35, 1400:0.3925, 2000:0.325, 3000:0.34, 5000:0.28, 7000:0.28, 10000:0.25, 15000:0.385, 20000:0.385, 30000:0.385, 50000:0.3925, 80000:0.3925, 100000:0.3925, 150000:0.3925, 200000:0.385, 300000:0.385, 500000:0.3775, 700000:0.37, 1000000:0.3625

C0 = 8 (would give "<= 7"): NOT reached. Floating-point optimum of the same model: 1.020217 with the standard 62-band schedule (10 sweeps, opt8.json) and 1.017877 with per-stage bands to k = 80 (opt8fine.json, last value in opt8fine.log at report time). Full exactness through 10^6 stages would subtract about 0.002 and a larger N before the tail about 0.0003, so the convex comparison misses '<= 7' by roughly 1.5-2 % of the budget. The extra loss comes from the single new modulus 10 (p = 5: E S rises from 0.100 to 0.200; with delta_5 = 0.2 the stop-loss is 0.084, and nu_5 = 1.25 then inflates every later stage). Stage losses at C0 = 8: p = 5: 0.084, 7: 0.125, 11: 0.055, 13: 0.073, 17: 0.071, 19: 0.075 (first moments 0.20, 0.21, 0.27, 0.26, 0.22, 0.21).

### 3.6 What is rigorous and what is numerical

Rigorous: everything in cert*.json: integer arithmetic at scale 10^60 with ceiling divisions in the
enlarging direction, exact rational deltas D/10^6, exact Fractions for p <= C0 and the tail
constants, explicit truncation remainder, explicit tail with the two cited Rosser-Schoenfeld
inequalities; the mathematical inputs are BBMST/CFT's measure construction and the repository's two
comparison lemmas, specialised as in 3.2 (that specialisation is a 5-line argument, written above,
not a new lemma). The verifier is standard library only (no numpy).
Numerical only: the choice of the delta schedule (floating-point coordinate descent); it is an input,
not a proof step. Not done: no proof-assistant formalisation; the Poisson-binomial recursion has an
independent non-rigorous (float) implementation that agrees, not a second rigorous one.

## 4. Verdict

LEAD (side result, certified; nothing toward 10,000 for the general problem).

Every covering system with distinct squarefree moduli has a modulus <= 10, i.e. one of 2, 3, 5, 6, 7, 10
(certified total 0.924567867 < 1 at C0 = 10, margin 0.0754; 0.922645 with full exactness). Intermediate
certificates: C0 = 14 (total 0.823434, margin 0.1766) and C0 = 12, i.e. "<= 11" (total 0.885626, margin
0.1144). This improves Cummings-Filaseta-Trifonov's 118 and uses only ingredients already proved in the
repository (comparison lemmas) plus BBMST/CFT's construction and two classical Rosser-Schoenfeld
inequalities; the integer verifier is self-contained (standard library, 20 s per certificate, 120 s with
full exactness). The result is rigorous modulo the correctness of verify_squarefree.py (validated against
CFT's published numbers (9) and (11) and against an independent float implementation) and of the
repository's two lemmas as specialised in 3.2. The next squarefree cutoff ("<= 7", C0 = 8) fails
numerically by 1.5-2 %, so the method as it stands saturates at 10.

Why this does not help the 10,000 target: the gain over CFT comes from the stop-loss replacing the
quadratic bound stage by stage, exactly what the repository's 14,484 certificate already exploits; the
squarefree problem is simply much easier (one modulus per subset, no prime powers, Euler factors 1+nu/p
instead of 1+nu/(p-1)).

## 5. Concrete next step for a follow-up agent

1. Write-up: a short section/note "squarefree minimum modulus <= 10" with the argument of 3.2, the
   inequality list 3.3, the tail (I8)-(I11), and the three certificate JSONs; copy verify_squarefree.py
   to code/ and cert*.json to data/ (main agent; this agent did not touch the repository). Worth a
   standalone arXiv note: 118 -> 10 for a problem CFT published in 2022/25 and that Dalton-Trifonov's
   Problem 4 asks about (whether the answer is 2).
2. Toward "<= 7" (C0 = 8, gap 1.5-2 %): the whole gap sits in the stages p = 5 and p = 7 (0.21 of the
   budget, moduli 10, 15, 30 and 14, 21, 35, 42, 70, 105, 210), where delta_2 = delta_3 = 0 makes the
   measure on (X_2, X_3) exactly uniform. An exact adversary over the residue choices of these <= 10
   moduli (class collisions inside the fibre make alpha < S; the repository's fable_adversary_toy*.py
   style, or an LP over the 30 x 5 / 210 x 7 cells), continued with the model from p = 11 on, would
   quantify whether 10 % of those two stages (0.02) is realisable. Cost: half a day. If it works, the
   same device at p = 5, 7, 11 might reach C0 = 6 ("<= 6"); C0 <= 4 ("<= 3") looks out of reach of
   moment methods (first moments at p = 3, 5 alone exceed 0.5).
3. Cheap polish if a margin is ever needed: full exactness (-0.002, 120 s), N = 4 x 10^6 before the
   tail (-0.0003), true per-stage delta optimisation (-0.003 estimated from the 62- vs 107-band runs).

## 6. Files produced

All in C:\Users\mekak\AppData\Local\Temp\claude\c--Users-mekak-code-minimum-modulus\23ed4d4e-d4c3-4fc1-93e5-481b63d3d735\scratchpad\agents\sqfree\ :
REPORT.md; sqfree_float_opt.py; verify_squarefree.py; summarise.py; cft_sched.json, cft118_check.json
(validation); opt14.json, opt13.json (checkpoint), opt12.json, opt10.json, opt8.json, opt8fine.json,
opt6.json and their .log files (schedules); cert14.json, cert12.json, cert10.json (certificates,
the proof objects, with cert*.log), cert10full.json (exact range through all 10^6 stages);
test14*.json, t.json (smoke tests).
