# Stochastic orders behind the comparison lemma, and what is left to gain

Historical note: the blanket sharpness statements below refer to the old
comparison information and should not be read as excluding geometry-aware
refinements. The later core comparison proves a new no-modulus-10 branch
at cutoff 6; see `core30-overlap-refinement-2026-10-26.md`. Claims below
that union accounting is the only discarded information are too broad:
actual conditional laws and last-coordinate collisions can also matter.

Date: 25 Oct 2026 (Fable). Question from the owner: the nested-event lemma is a special case of
the comonotonic convex-order theorem; does the stochastic-orders / dependence / actuarial
literature contain sharper results that would strengthen the certified bound M_sf <= 10
(and, in the parent repository, M <= 14,479)?

Short answer: no result of that literature sharpens the per-stage inequality, because the
inequality is already attained within the information it uses (Section 2). The only
information the certificate discards is cross-stage (a fibre already covered is charged
again). Section 3 gives the exact identity for that slack and Section 4 measures it at
C0 = 8: exhaustively over stages 5 and 7 it is worth 0.0015 against a deficit of 0.0179, and
at stage 11 an adapted adversary gives back only 0.001 to 0.003 of 0.055. The direction is a
lead for a multi-stage lemma, not a route that closes C0 = 8 with present tools.

## 1. Where the two lemmas sit in the literature

The nested-event lemma (`minimum-modulus/src/conditional-stoploss-improvement.tex`, Lemma
"Nested-event comparison") says: for events A_l with P(A_l) <= u_l and weights c_l >= 0,
E phi(sum c_l 1_{A_l}) <= E phi(sum c_l 1_{U <= u_l}) for every convex nondecreasing phi. It is
the Bernoulli case of two classical facts.

* Comonotonic sums are convex-largest. For X_l with fixed marginals F_l, sum_l X_l <=_cx
  sum_l F_l^{-1}(U). Meilijson and Nadas, "Convex majorization with an application to the
  length of critical paths", J. Appl. Probab. 16 (1979) 671-677 (their Theorem 2 gives it);
  Tchen, "Inequalities for distributions with given marginals", Ann. Probab. 8 (1980) 814-827
  (supermodular form: the upper Frechet bound maximises E f(X) for supermodular f, and
  phi(sum) is supermodular when phi is convex); Kaas, Dhaene, Vyncke, Goovaerts and Denuit,
  "A simple geometric proof that comonotonic risks have the convex-largest sum", ASTIN Bull.
  32 (2002) 71-80 (the elementary proof via stop-loss transforms). Our proof is the
  geometric one: E(S-t)_+ <= max_v (sum c_l min(v,u_l) - t v) = int (S*(u)-t)_+ du.
* Monotonicity in the marginals. Replacing P(A_l) = r_l by the bound u_l enlarges the
  comonotonic sum pointwise, so the comparison is in increasing convex order (stop-loss order),
  characterised by E(S-t)_+ <= E(T-t)_+ for all t: Shaked and Shanthikumar, Stochastic Orders
  (2007), Section 4.A; Dhaene, Denuit, Goovaerts, Kaas and Vyncke, "The concept of
  comonotonicity in actuarial science and finance: theory", Insur. Math. Econ. 31 (2002) 3-33,
  Section 2.

The product-comparison lemma conditions on the earlier coordinates, applies the nested-event
lemma to the last coordinate with the cap nu_j/p_j as the marginal bound, and inducts. In
actuarial terms this is the "improved comonotonic upper bound" of Kaas, Dhaene and Goovaerts,
"Upper and lower bounds for sums of random variables", Insur. Math. Econ. 27 (2000) 151-168,
with the conditioning variable Lambda = the earlier coordinates, i.e. conditional
comonotonicity in the sense of Cheung, "Improved convex upper bound via conditional
comonotonicity", Insur. Math. Econ. 42 (2008) 651-655. The parent repository's conditioning on
configurations of the primes 2 and 3 (and later through 53) is the same device with a coarser
Lambda, evaluated exactly on the retained states.

Attribution for the write-up: the nested-event lemma should be stated as the Bernoulli
specialisation of the comonotonic convex-order theorem (Meilijson-Nadas; Kaas et al. 2002),
and the product lemma as its conditional form. What is specific to the covering problem is
the identification of the conditional caps (BBMST's fibre rule) and of the cylinder structure
that makes the comonotonic configuration realisable (Section 2).

## 2. Why no stochastic-order result can sharpen the per-stage bound

The literature on sharper bounds under additional information falls into four groups. For
each, the information it needs is either already used or not available in the sieve.

| Result | Extra information used | In the covering sieve |
| --- | --- | --- |
| Improved comonotonic upper bound (Kaas-Dhaene-Goovaerts 2000; Cheung 2008) | conditional laws given Lambda | Used: Lambda = earlier coordinates (product lemma), and the exact small-prime configurations. |
| Improved Hoeffding-Frechet bounds (Nelsen; Rachev-Ruschendorf; Ruschendorf 2017 survey; Lux-Papapantoleon, Ann. Appl. Probab. 27 (2017)) | the copula on a subset, lower-dimensional marginals, or a measure of association | Pairwise and higher joint probabilities of cylinder events are known only as upper bounds (products of caps), and these are attained by the nested configuration, so the Frechet class is not reduced. |
| Variance or moment constraint on the sum (Cheung-Vanduffel, Scand. Actuar. J. 2013 (2) 103-118; Bernard-Ruschendorf-Vanduffel VaR bounds) | Var(sum) fixed or bounded | Cheung-Vanduffel prove that with Var(sum) fixed no convex-order upper bound exists; only tail (VaR) bounds improve. Our "truncation alpha <= 1" and modulated-cap routes were the covering analogues and gave 0.4-5 % (PROGRESS.md, closed routes). |
| Negative dependence (Christofides-Vaggelatou, J. Multivariate Anal. 88 (2004) 138-151: negatively associated vectors lie below independent ones in supermodular order; Puccetti-Wang, Statist. Sci. 30 (2015) 485-517; tree bounds for Bernoulli sums, Padmanabhan-Natarajan) | known negative dependence between summands | Within one coordinate, indicators of distinct residues are mutually exclusive, but the adversary is free to use one residue per coordinate, so the class contains the comonotonic configuration. |

Sharpness. At stage p the comparison variable is S = (1/p) #{J present : m_J p > C0} with the
primes q < p "present" independently with probability nu_q/q. The extremal configuration of the
Frechet class, nested across the earlier coordinates (all progressions through coordinate q
use the same residue) and distinct in the last coordinate (one class mod p per modulus), is a
legitimate family of progressions with distinct squarefree moduli, and the caps are attained
with equality on the uncovered points of every active fibre. Hence the bound is the exact
maximum over families of the per-stage charge, up to the fibres where the cap is not
attained. The measurements confirm it: the exhaustive stage-5/7 maximum is 0.211779 against
the comparison 0.211786 (`data/exact-two-stage-C0-8-D5-195000-D7-142500.json`); at stage 11
a CP-SAT adversary reaches 0.055092 against 0.055291 (`data/stage11-relaxation-tower-C0-8.txt`).
No inequality that only sees marginals, conditional caps, or dependence between the
indicators of one stage can beat an attained bound.

## 3. The information the certificate discards: an exact union identity

Let B_k be the set covered at stage k (prime p_k). B_k is determined by the first k
coordinates, and every later stage redistributes weight inside fibres of the first l-1 >= k
coordinates while preserving their totals, so P_n(B_k) = P_k(B_k). Therefore

    P_n( union_k B_k ) = sum_k P_k( B_k minus union_{l<k} B_l )
                       = sum_k sum_{x fibre, x not yet covered} (alpha_k(x) - delta_k)_+ P_{k-1}(x) / (1 - delta_k),

because a fibre at stage k is contained in, or disjoint from, every earlier B_l. The
certificate's total sum_k P_k(B_k) differs from this exactly by the charge sitting on fibres
that are already covered. This is the Hunter-Worsley (1976, 1982) improvement of the union
bound, written with the stage structure making the intersections explicit, and it is the
only place where cross-stage information enters.

Two facts control the size of the discarded charge.

* Covered fibres are light. A fibre covered at stage l with fraction alpha keeps the fraction
  (alpha - delta_l)/alpha of its covered weight, so a fibre covered just above the threshold is
  almost weightless, while a fully covered fibre (alpha = 1) keeps everything.
* Uncovered fibres are heavy only near earlier towers. The weight of an uncovered point is the
  cap product only if at every earlier stage its fibre was active (alpha_l > delta_l); where a
  stage did not act the factor is 1/p_l instead of 1/(p_l(1-delta_l)). The comparison model
  gives every point the full cap product.

The adversary therefore faces a trade-off the comparison does not see: heavy charge at stage k
needs the stage-k tower to sit on points that were near every earlier tower, and such points
are covered with substantial probability. The coupling toys of the parent repository (notes
`fable-notes-2026-10-23-multiagent.md`, agent "coupling") measured this trade-off at 2-20 %
of the charge depending on delta and depth, and found that the adversary who optimises the
union directly keeps 93-99 % of the per-stage sum.

## 4. Measurements at C0 = 8 (`code/union_accounting.py`, `data/union-accounting-C0-8*.json`)

Schedule: the tested float schedule `data/schedule8-float-only.json` (delta_5 = 0.195,
delta_7 = 0.1425, delta_11 = 0.36375), total 1.017877; deficit 0.017877. Stage charges: 5:
0.086957, 7: 0.124829, 11: 0.055291, 13: 0.072576, 17: 0.070122, 19: 0.075299, stages >= 13
together 0.806.

Part A, exhaustive and therefore a valid bound over all stage-5 and stage-7 families (same
enumeration as `exact_two_stage.py`: 32 orbits of stage-5 covered-set tables under
Aff(Z_2) x Aff(Z_3) x Aff(Z_5), all 810,000 stage-7 residue configurations, exact arithmetic):

| Quantity | Value |
| --- | --- |
| max over families of loss_5 + loss_7 (per-stage sum) | 0.211779 |
| max over families of loss_5 + (stage-7 charge on fibres uncovered after stage 5) | 0.210246 = 1219093/5798415 |
| rigorous gain of the union accounting at stages 5 and 7 | 0.001533 |
| double counting of the per-stage maximiser itself | 0.007280 (the union maximiser is a different family) |

So the union-aware certificate at stages 5 and 7 would lower the C0 = 8 total from 1.017877 to
1.016344. The adversary gives up 0.7 % of the two-stage charge to make the stage-7 cover
almost disjoint from the stage-5 cover.

Part B, specific families (lower bounds on what an adapted adversary keeps; not bounds over
all families):

| Family | per-stage 5+7 | covered mass after 7 | stage-11 adversary (all fibres) | same objective on uncovered fibres only |
| --- | --- | --- | --- | --- |
| nested/aligned stage 5 and 7 at residue 0 (stage11.py) | 0.175907 | 0.128616 | recorded CP-SAT 0.055092 (0.007799 of it on covered fibres, 14 %) ; best aligned 0.049984 | best aligned 0.049034 |
| two-stage maximiser | 0.211779 | 0.197601 | best aligned 0.050363 | best aligned 0.047813 |

CP-SAT relaxation of stage 11 (`stage11.py`, alpha <= min(1, n/11), 300 s per run, all FEASIBLE,
so each value is what an adversary found, not a proven maximum), `data/union-accounting-C0-8-cpsat.json`:

| Family (P_7) | objective = all fibres | of which on covered fibres | objective = uncovered fibres only | loss from being paid once |
| --- | --- | --- | --- | --- |
| nested/aligned stage 5 and 7 | 0.055092 | 0.007799 (14 %) | 0.054260 (its charge on covered fibres: 0.000076) | 0.000832 = 1.5 % |
| two-stage maximiser | 0.055280 | 0.011929 (22 %) | 0.052322 (0.001757 on covered fibres) | 0.002958 = 5.4 % |

Reading. The per-stage adversaries double count 14-22 % of their stage-11 charge, but an
adversary that is told it is only paid on uncovered fibres moves its residues and loses only
1.5-5.4 % (CP-SAT) or 2-5 % (aligned families). Summing such losses over all stages >= 7 would
give of the order of 1-5 % of 0.93, i.e. 0.01-0.05, which brackets the deficit 0.018 from both
sides; the lower end, reached by the family that is otherwise weaker, is below it. The direction is therefore not dead, but it needs two things that do not exist yet:

1. A rigorous multi-stage comparison that charges a fibre only once. The natural statement
   bounds sum_k E_{k-1}[1(x uncovered) (alpha_k - delta_k)_+] over all families; the
   coordinate induction of the product lemma breaks because "uncovered" depends on all
   earlier residues, not on a cap. The parent repository recorded the same obstacle (entries
   15 and 26 of its meta log).
2. A schedule re-optimised for the union objective; the present deltas equalise per-stage
   marginal losses, not union ones.

Exhaustive enumeration cannot replace the lemma beyond stage 7: stage 11 alone has 16
moduli over Z_210, and the stage-5/7/11 joint maximum of the union objective is out of reach
of exact enumeration (the stage-11 CP-SAT relaxation already needs minutes per P_7).

## 5. Verdict

* Attribution: the nested-event lemma is the Bernoulli case of the comonotonic convex-order
  theorem (Meilijson-Nadas 1979; Tchen 1980; Kaas et al. 2002); the product lemma is its
  conditional version (Kaas-Dhaene-Goovaerts 2000; Cheung 2008). This should be said in the
  write-up.
* No sharper stochastic-order result applies: the per-stage bound is attained within its
  information class, variance constraints give no convex-order improvement (Cheung-Vanduffel),
  and the available negative dependence is defeated by residue alignment.
* The one discarded piece of information is the union identity of Section 3. It is
  rigorously worth 0.0015 at C0 = 8 (stages 5-7) and measured at 2-5 % per later stage for
  adapted adversaries, against a deficit of 0.018. Open as a lead for a multi-stage lemma; the
  numbers do not justify expecting M_sf <= 7 from it without new ideas.
