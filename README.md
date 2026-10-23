# Minimum modulus of distinct squarefree covering systems

How large can the smallest modulus of a finite covering system be when the moduli are
distinct and squarefree? Write M_sf for the supremum. Known before this project:
M_sf >= 2 (Krukenberg 1971: explicit systems with least modulus 2 and lcm 210) and
M_sf <= 118 (Cummings, Filaseta, Trifonov 2025, building on the method of Balister,
Bollobas, Morris, Sahasrabudhe and Tiba, who proved that some modulus must be even).
Dalton and Trifonov (2022, Problem 4) ask whether M_sf = 2, i.e. whether every distinct
squarefree covering system contains the modulus 2. This repository aims to determine M_sf.
Progress and status live in [PROGRESS.md](PROGRESS.md); this file only says what the
project is and how it is laid out.

## What is being done

* **Upper side.** The convex stop-loss comparison developed for the general minimum-modulus
  problem (repository `minimum-modulus`) is specialised to squarefree moduli: coordinates
  Z_p, one modulus per subset of primes, Euler factors 1 + nu/p. An integer-arithmetic
  verifier certifies cutoffs C0 ("all moduli > C0 cannot cover"). Certified: C0 = 10, i.e.
  M_sf <= 10 (the only candidates left are 2, 3, 5, 6, 7, 10). The same method fails at
  C0 = 8 by about 2 %, so smaller cutoffs need exact arguments at the primes 5 and 7.
* **Lower side.** Exact SAT / CP-SAT searches for distinct squarefree covering systems with
  least modulus >= 3 (moduli dividing a squarefree lcm); nonexistence proofs for small lcm.
* **Determination.** M_sf = 2 if no system avoids the modulus 2 (Dalton-Trifonov Problem 4);
  otherwise an explicit system with least modulus in {3, 5, 6, 7, 10} must be found.

## Layout

| Path | Contents |
| --- | --- |
| `PROGRESS.md` | Current state, timeline, closed routes, open directions. |
| `code/verify_squarefree.py` | The certificate verifier (standard library only, about 20 s; 120 s with `--p-exact 1000000`). Usage: `python code/verify_squarefree.py data/schedule10.json out.json`. |
| `code/audit_verify.py` | Independent re-implementation of the verifier, written from the report alone; reproduces the C0 = 10 total exactly. |
| `code/schedule_opt.py` | Floating-point optimiser of the distortion schedule (input to the verifier, not a proof step). |
| `code/exact_two_stage.py`, `code/stage11.py` | Exact worst case of the stage-5 and stage-7 losses at C0 = 8 over all families, and the stage-11 adversary (toward "<= 7"). |
| `code/cpsat_cover.py` | CP-SAT search for distinct covering systems with moduli dividing N and minimum modulus >= M (squarefree N gives squarefree moduli). |
| `code/mc_stage.py` | Monte Carlo check of single stage bounds. |
| `data/cert*.json`, `data/schedule*.json` | Certificates (proof objects) and their schedules for C0 = 14, 12, 10; `cert10-full.json` uses exact stages through p = 15,485,863. |
| `data/cft118-*.json`, `data/audit-*.json` | Validation against Cummings-Filaseta-Trifonov's published numbers at C0 = 118, and the independent audit outputs. |
| `notes/` | Dated reports. |
| `literature/SOURCES.md` | References. |

## Provenance

The method, code and certificates were produced in October 2026 by Claude (Fable 5.1)
agents working with the repository owner; the verifier was audited by an independent
re-implementation but not by a proof assistant. The mathematical inputs are the
Balister-Bollobas-Morris-Sahasrabudhe-Tiba / Cummings-Filaseta-Trifonov distorted measure,
two comparison lemmas proved in the `minimum-modulus` repository
(`src/conditional-stoploss-improvement.tex`), and two inequalities of Rosser and Schoenfeld (1962).
