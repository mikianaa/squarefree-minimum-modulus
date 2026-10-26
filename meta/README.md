# Meta-analysis log of the agents' research progress

Records how the two agents (Codex, Claude Code "Fable") advance the determination of
M_sf, the largest possible minimum modulus of a distinct squarefree covering system, for
later analysis of agent research behaviour. Separate from `PROGRESS.md` (which states the
mathematics) and from the dated notes (which give the details). Same layout and schema as
`meta/` in the parent repository `minimum-modulus`, with the fields adapted to this problem.

## Files

| File | Role |
| --- | --- |
| `agent-progress.jsonl` | Source of truth. One JSON record per milestone or attempt, hand-written at the time. Lines starting with `#` are comments. |
| `meta_stats.py` | Reads the local Claude Code and Codex transcripts, aggregates token usage, and renders the log. Standard library only. |
| `token-usage.json` | Output: tokens, API-list-price cost and wall-clock per session and per day, per agent. |
| `agent-progress.md` | Output: the log as a table with the cumulative tokens spent at each milestone. |

Regenerate the two outputs with `python meta/meta_stats.py` from the repository root.

## Which transcripts count

Transcripts are read from `~/.claude/projects/<project>/**/*.jsonl` (usage per API call,
subagent transcripts included) and `~/.codex/sessions/**/rollout-*.jsonl` (per-turn
`token_count` events), from 23 Oct 2026. The project was split out of `minimum-modulus` on
24 Oct 2026, so its first sessions were started from the parent repository and live under
the parent's transcript folder. `CLAUDE_SOURCES` in `meta_stats.py` lists exactly what is
adopted:

* all transcripts of this repository's own project folder;
* from the parent's multi-agent run of 23 Oct (session `23ed4d4e`), only the three
  squarefree subagents `sqfree`, `sqaudit` and `sqfree7`. The `transfer` agent that found
  the numerical lead worked on the general problem and stays in the parent's log;
* the 24 Oct session `a5a2cf16`, which created this repository, in full (its first quarter
  hour reviewed the parent's progress).

Codex sessions are selected by working directory (`squarefree-minimum-modulus`). This
includes the hourly `codex-auto-review` child threads, Codex's approval reviewer, which are
reported under their own model name.

## Metrics and how each is recorded

| Metric (user's list) | Field(s) in `agent-progress.jsonl` | How it is measured |
| --- | --- | --- |
| Final progress: rigorously proved bounds | `bounds_after.upper_certified`, `.lower`, `.excluded_primorial`, `.primes_needed` | State of knowledge after the entry. `upper_certified` is the best cutoff bound on M_sf passed by the integer verifier (118 published, 10 here); `lower` stays at Krukenberg's 2; `excluded_primorial` is the largest primorial N for which "the divisors >= 3 of N cannot cover" is established (by a solver verdict for 210 and 2310, by audited integer certificates from 30030 on); `primes_needed` is the least number of distinct prime factors a counterexample must have, once the BBMST even-modulus theorem and prime compression are invoked (null before that). |
| Compute resources: tokens, cost, time | `compute.tokens` (computed), `compute.cost_usd` (computed), `compute.wall_clock_hours` (hand-entered, null if unknown) | Tokens and cost come from the transcripts via `meta_stats.py`; cost is API list price (Fable 5.1: 10 / 50 USD per MTok in / out, cache read 0.25; Codex rate unset), so it is an estimate, not a bill. |
| Path to discovery: at which token the improvement appeared | computed column "tokens at discovery" in `agent-progress.md` | Cumulative tokens of that agent (and of both agents) on this project up to the entry's `timestamp`. Entry 1 was found by an agent counted in the parent repository, so its cumulative count here is zero. |
| Search efficiency: candidates and success rate | `search.candidates`, `search.successes` | Hand-entered: number of distinct candidates / designs / runs / cases tried in the attempt, and how many led to a verified result. Null when the attempt had no enumerable candidates. |
| Mathematical method | `method` (list of tags) | `parameter-optimization`, `inequality-improvement`, `new-lemma`, `conditioning`, `construction-patch`, `exact-search` (SAT / CP-SAT / exhaustive enumeration), `capacity-bound` (counting, plane-load, compression and LP relaxations of residual capacity), `density-bound`, `audit`, `numerical-experiment`, `literature`, `tooling`. |
| Cause of failure | `failure_cause` | `numerical` (a numerical estimate was wrong, or an exact adversary shows the estimate is already tight), `logical-gap` (no proof route for an idea), `search-dead-end` (solver / search limits, run inconclusive), `capacity-shortfall` (a counting or density obstruction kills the route), `relaxation-gap` (the relaxation admits a witness that is not a covering system, so it cannot exclude the case), `none`. |

Other fields: `id`, `timestamp` (ISO 8601 with offset; `approx: true` when reconstructed
from file dates or notes rather than logged), `agent` (`codex` or `claude-code`), `model`,
`side` (`upper`, `lower`, `side`), `kind` (`certificate`, `construction`, `obstruction`,
`lead`, `refutation`, `audit`, `tool`), `result.summary`, `result.rigor` (`certified` =
integer verifier or independent integer audit passed; `solver` = exact SAT / CP-SAT verdict
without a checked proof log; `conditional`; `candidate`; `none`), `notes`, `files`.

## Adding an entry

Append one line to `agent-progress.jsonl` when a milestone is reached or an attempt is
closed, then run `python meta/meta_stats.py`. Keep entries factual; the computed columns
should never be written by hand.

## Audit comparison (cheap model vs reference)

Added 26 Oct 2026. Tests whether a short fresh-context audit on a cheaper model (Sonnet,
`.claude/agents/auditor.md`) reproduces a reference audit, and what it costs. Rules for the
agents are in `CLAUDE.md` and `AGENTS.md`.

| File | Role |
| --- | --- |
| `audit-comparison.jsonl` | Source of truth, one hand-written record per comparison. |
| `compare_audits.py` | Computes agreeing digits, difference, and the tokens / list-price cost of the candidate run from the transcripts (by model and time window). Standard library only. |
| `audit-comparison.md` | Output table. |

Record fields: `id`, `task`, `threshold_digits` (default 12), `reference` (`agent`, `model`,
`source` file, `value` as an exact decimal string), `candidate` (`agent`, `model`, `value`,
`window` = [start, end] ISO timestamps of the run, `runtime_s`, `files`, `issues`).
Regenerate with `python meta/compare_audits.py`. The candidate must not be shown the
reference value before its run; record any disagreement as well as agreement.
