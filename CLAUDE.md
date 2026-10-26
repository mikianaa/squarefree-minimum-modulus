# Squarefree minimum modulus: working rules for Claude Code

Goal: determine M_sf, the largest possible minimum modulus of a finite covering system
with distinct squarefree moduli. Mathematics: `PROGRESS.md` (state) and `notes/` (details).
Agent behaviour log: `meta/` (see `meta/README.md`). Codex works in this tree at the same
time and leaves large uncommitted changes; commit only your own files.

## Model and context roles

| Role | Model | Context |
| --- | --- | --- |
| Discovery: new lemmas, direction, deciding what a result means | main session (Fable) | keep short; start a new session at each milestone |
| Independent audit, integer re-implementation, regression checks | `auditor` subagent (Sonnet) | fresh each time; input is the note and data only |

Why: most tokens in this project are cache reads (over 90 % of the total), i.e. a long
context re-read every turn. Splitting sessions and pushing reading-heavy checks into short
fresh contexts reduces that.

* Hand state over through files, not through a long session: `PROGRESS.md`, `notes/`,
  `meta/agent-progress.jsonl`. Before ending a session, make sure they say what the next
  session needs.
* Delegate audits with the `auditor` subagent. Give it the note path and the input data
  path, not the original verifier and not the expected value.
* A Sonnet audit only counts toward a certified claim after it has been compared with an
  existing audit (see below). Until a task type has been compared, also run the reference
  audit.

## Comparison protocol (cheap vs reference audit)

1. Pick a task that already has a reference result (e.g. the C0 = 10 certificate).
2. Run the `auditor` on it. State every run parameter (p_exact, C0, schedule file) in the prompt, but not the reference value.
3. Record it in `meta/audit-comparison.jsonl` (schema in `meta/README.md`) with the
   reference and the candidate values as exact strings and the time window of the run.
4. Run `python meta/compare_audits.py`. It computes the agreeing digits and the tokens and
   list-price cost of the window; never write those by hand.
5. Disagreement is a finding, not a failure to hide. Log it too.

## Meta log

On each milestone or closed attempt, append a record to `meta/agent-progress.jsonl` and
run `python meta/meta_stats.py`. Computed columns are never written by hand.

## Git

Commit only files you created or changed for the task. For README / PROGRESS hunks that
Codex also touched, stage through the index rather than `git add <file>`. Ask before
changing whether `meta/` is tracked.
