# Squarefree minimum modulus: working rules for Codex

Goal: determine M_sf, the largest possible minimum modulus of a finite covering system
with distinct squarefree moduli. State: `PROGRESS.md`; details: `notes/`; agent behaviour
log: `meta/` (see `meta/README.md`). Claude Code works in this tree at the same time, so do
not rewrite or revert files you did not create or change for your task.

## Context discipline

Most of this project's tokens are cache reads: one very long session re-reads its whole
context every turn (the 25 Oct Codex session was 108 M tokens, 99 % cache reads).

* Start a new session at each milestone or closed attempt instead of continuing one long
  session.
* Carry state in files: update `PROGRESS.md`, write the dated note in `notes/`, and append
  a record to `meta/agent-progress.jsonl` before you stop. A fresh session must be able to
  continue from those files alone.
* Delegate reading-heavy checks (independent re-implementation, regression sweeps) to a
  short fresh sub-task that is given the note and the input data only, not the original
  verifier and not the expected value.

## Independent audits

An audit is independent only if it is rebuilt from the written argument. Do not tune the
audit until it matches the reference. If it disagrees, record the disagreement in
`meta/audit-comparison.jsonl` (schema in `meta/README.md`); `python meta/compare_audits.py`
computes the agreement and the cost.

## Logging

Computed columns of the meta log (`tokens`, `cost`) are produced by
`python meta/meta_stats.py` and are never written by hand.
