# Progress

Status of the determination of M_sf, the largest possible minimum modulus of a finite
covering system with distinct squarefree moduli. Published state of the art: 2 <= M_sf <= 118.
Last updated 24 October 2026 (10:10 JST).

## Current state

| Side | Result | Status | Where |
| --- | --- | --- | --- |
| Upper | M_sf <= 10 | Certified 24 Oct 00:40 by the integer verifier (total 0.92457 < 1, margin 0.075; full exactness 0.92265); independently re-implemented verifier agrees to 12 decimal digits (24 Oct 00:55); re-run from this repository reproduces the integer total exactly | `data/cert10.json`, `data/cert10-full.json`, `data/audit-cert10.json`, `code/verify_squarefree.py`, `code/audit_verify.py` |
| Upper target | M_sf <= 7 (C0 = 8) | Not reached: the comparison gives 1.0078 with full exactness (1.018 in the float optimiser). The exact worst case over all families at p = 5, 7 equals the comparison value (0.211779 vs 0.211786), so the comparison is tight there and exact small-stage handling gains nothing at 5, 7. Truncation alpha <= 1 at p >= 11 has a numerical ceiling of 0.905 (saving 0.10) but is NOT realisable: for the tower measure P_7 a spread stage-11 family (CP-SAT over the 16 residue choices on Z_210) reaches 0.05509 of the comparison's 0.05529 (99.6 %), so the cap saves 0.0002 of the 0.0049 ceiling. The method saturates at 10 | `data/schedule8-float-only.json`, `data/exact-two-stage-C0-8-D5-195000-D7-142500.json`, `code/exact_two_stage.py`, `code/truncation_ceiling.py` |
| Lower | M_sf >= 2 | Krukenberg 1971 (lcm 210); 14- and 15-class systems mod 210 verified | literature, `data/cpsat-lcm210-min2-sat.json` |
| Lower (obstructions) | least modulus 3 needs a prime factor >= 13 | Codex (15 Oct 2026, parent repository): rigorous counting obstruction, every distinct squarefree system with moduli >= 3 and primes <= 11 misses >= 22 residues mod 2310; best heuristic through prime 17 misses 1.575 % (126 congruences); minimum 30 needs a prime >= 37 (density) | `minimum-modulus`: `src/construction-search.tex` (sections on squarefree 30 and 3), `code/squarefree_parity_search.py`, `code/squarefree_parity_repair.py`, `data/squarefree-parity-3-*.json` |
| Lower target | a system with least modulus >= 3 | Open (Dalton-Trifonov Problem 4). CP-SAT with the lossless normal form (all divisors >= 3 used, class of each odd prime p at 0, class of 2p at 1 mod p, parity of the smallest 2p fixed): lcm 210 UNSAT 0.1 s, lcm 2310 UNSAT 27 s (agrees with Codex's counting obstruction); lcm 30030 running since 24 Oct 09:05 (limit 2 h). Without the normal form 2310 and 30030 were UNKNOWN after 15 min | `code/cpsat_cover.py` |

## Timeline

| Date (2026) | Side | Result | Files |
| --- | --- | --- | --- |
| 23 Oct | upper | The convex stop-loss comparison of the general problem, applied to squarefree moduli, gives minimum modulus <= 14 numerically (published 118). Found by the "transfer" agent of the multi-agent run in the `minimum-modulus` repository. | `notes/sqfree-agent-report-2026-10-23.md` |
| 24 Oct | upper | Certificates C0 = 14 (total 0.8234), C0 = 12 ("<= 11", 0.8856), C0 = 10 (0.9246). Verifier validated on CFT's own delta table at C0 = 118 (stages 1-7 exactly 194/1365, tail 0.000402960697 vs CFT 0.000402960685). C0 = 8 fails numerically at 1.018-1.020. | `data/cert14.json`, `data/cert12.json`, `data/cert10.json`, `data/cft118-check.json` |
| 24 Oct | upper | Independent audit: a second verifier written from the report alone reproduces total 0.924567866690 (12 digits) and the part sums; outward roundings differ beyond the 20th digit. | `code/audit_verify.py`, `data/audit-cert10.json` |
| 24 Oct | upper | Exact adversary over all stage-5 families (orbits under Aff(Z_2) x Aff(Z_3) x Aff(Z_5)) and all 810,000 stage-7 residue configurations at C0 = 8: max loss_5 + loss_7 = 409327/1932805 = 0.2118 at (D5, D7) = (0.195, 0.1425). Run interrupted (machine out of memory) before the stage-11 continuation. | `code/exact_two_stage.py`, `code/stage11.py`, `data/exact-two-stage-*.json` |
| 24 Oct | lower | Brute CP-SAT over squarefree lcm with minimum modulus 3: 210 UNSAT (0.6 s; min 2 SAT with all 15 moduli, Krukenberg type), 2310 and 30030 UNKNOWN after 900 s with 8 workers. | `data/cpsat-lcm210-*.json`, `code/cpsat_cover.py` |
| 24 Oct | repo | Project split out of `minimum-modulus` into this repository. | |
| 24 Oct | lower | Slice extension: a cover with lcm Lq (q a new prime) consists of a q-free family over L with residual R plus the shapes qd (d | L) distributed over the q slices; one slice is covered by q itself, each other slice must cover R. Minimal residuals at 210 (22 points) and 2310 (121 points) are spread, so Sum_d max_a |R cap (a + dZ)| is about 4.5 |R|, far below the 12 |R| needed at q = 13. The most concentrated residual with >= 22 points at L = 2310 is a full class mod 105 (Codex's bound 22 = 2310/105 is exactly this size), and it is the ONLY single-class residual that passes the count Sum_d gcd(m,d)/d >= 12 (12.09); the slice assignment for it exists (CP-SAT, 4 s). Hence a cover with least modulus 3 and lcm 30030 exists if Z_2310 can be covered by all divisors >= 3 with the modulus 105 used twice. That test is UNKNOWN after 300 s; 1 h run started 10:05. | `code/slice_extend.py`, `code/duplicate_test.py`, `data/duplicate-2310-105-300s-unknown.json` |
| 24 Oct | upper | Exact worst case at p = 5, 7 (C0 = 8) equals the comparison value: stage 5 exactly 2/23 both, stage 7 0.1248222 exact vs 0.1248289 comparison. Route closed. Truncation ceiling at C0 = 8 computed: 0.905 total, savings per stage 0.005 (11), 0.008 (13), 0.004 (17), 0.005 (19). | `code/truncation_ceiling.py` |
| 24 Oct | upper | Stage-11 exact adversary (relaxation alpha <= min(1, n/11), CP-SAT 540 s) for the aligned stage-5/7 measure at C0 = 8: a spread family realises 0.055092 against the comparison 0.055291 and the truncated model 0.050363. Truncation is dead for squarefree moduli too; the upper side saturates at 10. | `data/stage11-relaxation-tower-C0-8.txt`, `code/stage11.py` |
| 24 Oct | lower | Minimum uncovered count for least modulus 3: lcm 210 exactly 22 of 210 (10.5 %, OPTIMAL 0.3 s); lcm 2310 at most 121 of 2310 (5.2 %, FEASIBLE after 540 s, solver bound 19; Codex's rigorous bound 22). | `code/cpsat_min_uncovered.py`, `data/min-uncovered-lcm2310-feasible121.json` |
| 24 Oct | lower | Lossless normal form added to the CP-SAT search; lcm 2310 with least modulus 3 proved infeasible in 27 s. | `code/cpsat_cover.py`, `data/cpsat-lcm2310-min3-normalized-unsat.json` |

## Closed routes (do not retry)

* Exact adversary at the stages p = 5 and 7 for C0 = 8: the comparison is already tight there (the aligned tower realises it), gain 7e-6.
* Truncation alpha <= 1 at p >= 11: spread families realise 99.6 % of the untruncated comparison at stage 11 (same phenomenon as in the general problem). No moment-method route below 10 remains.

* Parameter tuning of the convex comparison below C0 = 10: the floating-point optimum at C0 = 8 is 1.018 with per-stage bands; full exactness (-0.002) and a larger N before the tail (-0.0003) do not close the gap.
* Moment methods alone for C0 <= 4 ("<= 3" or Problem 4): first moments at p = 3 and 5 already exceed 0.5.

## Open directions

* Upper: only arguments outside the stage-wise measure paradigm (coverage itself, global residue consistency) could go below 10. Within the paradigm the adversary realises the comparison to within 0.4 % at every examined stage.
* Lower: exact search for a distinct squarefree system with least modulus >= 3 over squarefree lcm. Brute CP-SAT stops at lcm 2310; Codex's parity-layer normalisation (class mod p at 0, class mod 2p at 1 mod p) plus counting gives the obstruction at 2310 and should be the basis for 30030 and 510510. Finite obstructions by lcm; the first lcm that is not excluded tells where a construction could live.
* Write-up: a short note (118 -> 10) with the argument, the inequality list and the certificates.
