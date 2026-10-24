"""Case split for lcm N = L q (least modulus 3): WLOG slice 1 holds the fewest shapes, hence at most
k = floor((tau(L)-1)/(q-1)); enumerate its exact shape set S (|S| <= k; |S| = 0 is impossible since the
slice would be uncovered) and solve each case with cpsat_cover.solve.  All cases UNSAT => no cover.
Cases whose union of progressions is smaller than the residual lower bound are skipped as trivially UNSAT.
Usage: python case_split_slice1.py N k rlb [--time 120] [--workers 8] [--output cases.jsonl]
"""
import argparse, itertools, json, time
from cpsat_cover import divisors, solve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("k", type=int); ap.add_argument("rlb", type=int)
    ap.add_argument("--time", type=float, default=120); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output", default="cases.jsonl"); ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    N = a.N
    q = max(p for p in range(2, N + 1) if N % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1)))
    L = N // q
    shapes = [d for d in divisors(L) if d > 1]
    done = set()
    if a.skip_done:
        try:
            for line in open(a.output):
                e = json.loads(line)
                if e["status"] in ("UNSAT", "SAT", "SKIP"): done.add(tuple(e["slice1"]))
        except FileNotFoundError: pass
    cases = [c for r in range(1, a.k + 1) for c in itertools.combinations(shapes, r)]
    print(f"N={N} q={q} L={L}: {len(cases)} cases", flush=True)
    with open(a.output, "a") as out:
        for S in cases:
            if S in done: continue
            if sum(L // d for d in S) < a.rlb:
                e = dict(slice1=list(S), status="SKIP", reason="union smaller than residual bound")
            else:
                r = solve(N, 3, a.time, a.workers, normalize=True, slice_cuts=a.rlb, small_slice=a.k, slice1_shapes=S)
                e = dict(slice1=list(S), status=r["status"], seconds=r["seconds"])
                if "system" in r: e["system"] = r["system"]
            out.write(json.dumps(e) + "\n"); out.flush()
            print(json.dumps({k: v for k, v in e.items() if k != "system"}), flush=True)


if __name__ == "__main__":
    main()
