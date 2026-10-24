"""Slice extension: from a family over squarefree L whose distinct part misses a residual R, build a
covering system over L*q (q a new prime) in which the shapes q*d (d | L) are distributed over the
q slices x = c mod q: one slice is covered by the modulus q itself, every other slice must cover R
with the shapes assigned to it.  Also verifies any produced system directly on Z_{Lq}.

Modes:
  python slice_extend.py class L m q            R = a full class mod m (slice-assignment test only)
  python slice_extend.py file dup.json q out.json   family from duplicate_test.py (uses the first SAT entry),
                                                  R = residual of the distinct part, builds and verifies the cover
"""
import json, sys, time
from ortools.sat.python import cp_model
from cpsat_cover import divisors


def assign_slices(L, R, q, time_limit=600, workers=8):
    """R: list of residues mod L.  Returns dict (d -> (slice c, residue a mod d)) covering R in every slice
    c != c0 (c0 gets the modulus q, shape d = 1), or None."""
    shapes = [d for d in divisors(L)]                     # d = 1 is the modulus q itself
    mdl = cp_model.CpModel()
    y = {(d, c, a): mdl.NewBoolVar(f"y{d}_{c}_{a}") for d in shapes for c in range(q) for a in range(d)}
    for d in shapes:
        mdl.AddExactlyOne(y[(d, c, a)] for c in range(q) for a in range(d))
    mdl.Add(y[(1, 0, 0)] == 1)                            # WLOG the modulus q covers slice 0
    for c in range(1, q):
        for x in R:
            mdl.AddBoolOr([y[(d, c, x % d)] for d in shapes if d > 1])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    st = s.Solve(mdl)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, s.StatusName(st)
    return {d: next((c, a) for c in range(q) for a in range(d) if s.Value(y[(d, c, a)])) for d in shapes}, s.StatusName(st)


def crt(L, r, q, c):
    """x mod Lq with x = r mod L and x = c mod q."""
    for k in range(q):
        if (r + k * L) % q == c: return r + k * L
    raise ValueError


def verify(N, system):
    cov = bytearray(N)
    for d, r in system:
        for x in range(r % d, N, d): cov[x] = 1
    assert len({d for d, _ in system}) == len(system), "moduli not distinct"
    return all(cov)


def main():
    mode = sys.argv[1]
    if mode == "class":
        L, m, q = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        R = [x for x in range(L) if x % m == 0]
        t = time.time(); asg, st = assign_slices(L, R, q)
        print(f"L={L} R=class 0 mod {m} (|R|={len(R)}) q={q}: {st} ({time.time()-t:.1f}s)")
        if asg:
            per = {}
            for d, (c, a) in asg.items(): per.setdefault(c, []).append(d)
            for c in sorted(per): print(f"  slice {c}: shapes {sorted(per[c])}")
        return
    dup = json.load(open(sys.argv[2])); q = int(sys.argv[3]); out = sys.argv[4] if len(sys.argv) > 4 else None
    entries = dup if isinstance(dup, list) else [dup]
    e = next((e for e in entries if e["status"] == "SAT"), None)
    if e is None: print("no SAT entry"); return
    L, m = e["L"], e["m_dup"]
    fam = [(s["modulus"], s["copy"], s["residue"]) for s in e["system"]]
    distinct = [(d, r) for d, i, r in fam if i == 0]
    assert len({d for d, _ in distinct}) == len(distinct)
    R = [x for x in range(L) if not any(x % d == r for d, r in distinct)]
    print(f"L={L} m_dup={m}: distinct part misses |R|={len(R)}")
    asg, st = assign_slices(L, R, q)
    print("slice assignment:", st)
    if not asg: return
    system = [(d, r) for d, r in distinct]
    for d, (c, a) in asg.items():
        system.append((q * d, crt(d, a, q, c) if d > 1 else c))
    N = L * q
    ok = verify(N, system)
    print(f"system over {N}: {len(system)} moduli, min modulus {min(d for d, _ in system)}, covers: {ok}")
    if out: json.dump(dict(N=N, min_modulus=min(d for d, _ in system), covers=ok,
                           system=[dict(modulus=d, residue=r) for d, r in sorted(system)]), open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
