"""Compare cheap-model audits with reference audits.

Reads meta/audit-comparison.jsonl (one hand-written record per comparison) and writes
meta/audit-comparison.md.  Agreement and cost are computed, never written by hand:

  * agreeing digits / absolute difference from the two exact value strings;
  * tokens and API-list-price cost of the candidate run from the Claude Code
    transcripts, selected by model and by the run's time window.

Record schema (one JSON object per line; lines starting with # are comments):

  id, task, note                       free text
  reference: {agent, model, source, value}   value = exact decimal string
  candidate: {agent, model, value, window: [start_iso, end_iso], runtime_s, files, issues}
  threshold_digits                     digits that must agree for "match" (default 12)

Run from the repository root:  python meta/compare_audits.py
Standard library only.
"""
from __future__ import annotations

import json
import os
import sys
from decimal import Decimal, getcontext

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meta_stats as ms  # noqa: E402

getcontext().prec = 80
SRC = os.path.join(ms.META, "audit-comparison.jsonl")
OUT = os.path.join(ms.META, "audit-comparison.md")


def agreeing_digits(a: str, b: str) -> tuple[int, Decimal]:
    x, y = Decimal(a), Decimal(b)
    diff = abs(x - y)
    if diff == 0:
        return 99, diff
    scale = max(abs(x), abs(y))
    digits = 0
    if scale:
        digits = max(0, int((-(diff / scale).log10()).to_integral_value(rounding="ROUND_FLOOR")))
    return digits, diff


def window_usage(model: str, start: str, end: str) -> dict:
    t0, t1 = ms.parse_ts(start), ms.parse_ts(end)
    tot = {"input": 0, "cache_write": 0, "cache_read": 0, "output": 0}
    rates = ms.CLAUDE_RATES.get(model)
    cost = 0.0
    for e in ms.claude_events():
        # subagent transcripts only: the main session may run on the same model
        if e["model"] != model or not (t0 <= e["ts"] <= t1) or "subagents" not in e["file"]:
            continue
        cw = e["cache_write_5m"] + e["cache_write_1h"]
        tot["input"] += e["input"]
        tot["cache_write"] += cw
        tot["cache_read"] += e["cache_read"]
        tot["output"] += e["output"]
        if rates:
            cost += (e["input"] * rates["input"] + e["output"] * rates["output"]
                     + e["cache_read"] * rates["cache_read"]
                     + e["cache_write_5m"] * rates["cache_write_5m"]
                     + e["cache_write_1h"] * rates["cache_write_1h"]) / 1e6
    tot["total"] = sum(tot.values())
    tot["cost_usd"] = round(cost, 2) if rates else None
    return tot


def main() -> None:
    recs = []
    if os.path.exists(SRC):
        with open(SRC, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    recs.append(json.loads(line))
    rows = []
    for r in recs:
        ref, cand = r["reference"], r["candidate"]
        digits, diff = agreeing_digits(ref["value"], cand["value"])
        need = r.get("threshold_digits", 12)
        u = window_usage(cand["model"], *cand["window"]) if cand.get("window") else None
        rows.append((r, digits, diff, need, u))
    lines = [
        "# Cheap-model audit vs reference audit (rendered by meta/compare_audits.py)",
        "",
        "Agreement and tokens are computed from the exact value strings and the transcripts.",
        "`match` = at least the threshold number of leading digits agree.",
        "",
        "| # | task | reference (agent / model) | candidate (agent / model) | agreeing digits | abs. diff | verdict | candidate tokens (total / output) | candidate cost USD | runtime s | issues |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r, digits, diff, need, u in rows:
        ref, cand = r["reference"], r["candidate"]
        dstr = ">= 80" if digits >= 99 else str(digits)
        verdict = "match" if digits >= need else "MISMATCH"
        tok = f"{u['total']:,} / {u['output']:,}" if u else "-"
        cost = f"{u['cost_usd']}" if u and u["cost_usd"] is not None else "-"
        lines.append(
            f"| {r['id']} | {r['task']} | {ref['agent']} / {ref['model']} | {cand['agent']} / {cand['model']} "
            f"| {dstr} | {diff:.3E} | {verdict} | {tok} | {cost} | {cand.get('runtime_s', '-')} "
            f"| {cand.get('issues') or '-'} |")
    if not rows:
        lines.append("| - | no comparisons recorded | | | | | | | | | |")
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {OUT} ({len(rows)} comparison(s))")


if __name__ == "__main__":
    main()
