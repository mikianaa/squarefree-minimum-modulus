# Cheap-model audit vs reference audit (rendered by meta/compare_audits.py)

Agreement and tokens are computed from the exact value strings and the transcripts.
`match` = at least the threshold number of leading digits agree.

| # | task | reference (agent / model) | candidate (agent / model) | agreeing digits | abs. diff | verdict | candidate tokens (total / output) | candidate cost USD | runtime s | issues |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | C0=10 certificate total, run as first specified (p_exact not named in the prompt, auditor took 100000 from the schedule file) | claude-code / claude-fable-5-1 | claude-code / claude-sonnet-5-5 | 2 | 9.097E-3 | MISMATCH | 260,437 / 14,226 | 0.34 | 3 | Task statement omitted p_exact; the schedule file's own field says 100000 while the reference certificate uses 1000000. Auditor flagged this itself. Not a defect of the audit. |
| 2 | C0=10 certificate total, same auditor code with p_exact = 1000000 (override of one constant by the main session, no new Sonnet run) | claude-code / claude-fable-5-1 | claude-code / claude-sonnet-5-5 | 16 | 3.043E-17 | match | - | - | 6 | Not independent of entry 1: same code, one constant changed after the reference was known. |
| 3 | C0=10 certificate total, independent rerun with all parameters stated (p_exact 1000000, n_stages 1000000, L 1000000), fresh context, no access to existing code or results | claude-code / claude-fable-5-1 | claude-code / claude-sonnet-5-5 | 16 | 3.043E-17 | match | 193,109 / 11,742 | 0.23 | 8.5 | Auditor trusts only 9 digits (cross-checked against the report's printed precision); it did not re-prove the comparison lemma or Rosser-Schoenfeld inequalities. |
