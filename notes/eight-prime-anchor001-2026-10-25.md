# The entire 001 branch is excluded

## Result

At N=9699690=2*3*5*7*11*13*17*19, every normalized family in low-anchor
pattern (b_3,b_5,b_7)=001 leaves at least **305** residues uncovered.
This includes every possible 210-core and all sixteen higher-anchor
patterns (b_11,b_13,b_17,b_19). The proof uses only the core count frontier
and the arbitrary-subset compression bound. It does not use the orbit or
inclusion reduction, nor an enumeration of actual final core residuals.

Together with the audited 000 result, this excludes b_3=b_5. Consequently,
any hypothetical covering system with distinct moduli greater than 1
dividing N and omitting modulus 2 must contain both **6 and 10**, and the
classes for those two moduli must have opposite parities.

The presence assertion follows by adding a missing class of modulus 6 or
10 in the same parity as the other one. If both are absent, add them in
the same parity. This preserves coverage and still omits modulus 2,
contradicting the excluded branches after normalization. Parity reversal
handles b_3=b_5=1. The statement concerns divisors of N; for an arbitrary
eight-prime support it applies to the compressed system, without claiming
that its original moduli literally include 6 and 10.

The 010 and 011 branches remain open. Larger prime supports remain open.
The global conjecture and the global lower bound of eight primes for a
counterexample are unchanged.

## Count-based certificate

The two anchor boxes for 001 have dimensions (1,3,6) and (2,4,5). Every
210-core residual has counts componentwise above one of

```
(0,29), (1,25), (2,23), (3,22), (4,21),
(7,19), (8,18), (9,17), (10,16).
```

As in `notes/eight-prime-anchor-branch-2026-10-25.md`, choose a subset of
the actual residual with one of these exact count pairs. The number of
uncovered points on this subset is a lower bound for the original family.
This step does not require monotonicity of the final dual expression.

Let C(m0,m1) be the sharp compression upper bound on the sum of the 14
non-unit core projection capacities for arbitrary subsets of these two
boxes. Tail groups, which also have shapes 1 and 2, use

    C_16(m0,m1) = C(m0,m1) + m0 + m1 + max(m0,m1).

The conditional-tail local-state inequality and its integer dual are the
same as in `notes/conditional-tail-obstruction-2026-10-25.md`. The local
point inequalities depend only on the higher-prime pattern, not the core.
Their concentration right-hand sides are recomputed from C on the
rectangle 0<=m_e<=r_e for each frontier pair (r0,r1).

New duals were optimized directly for these abstract count/capacity
inputs. A finite pool supplies a positive lower bound for all 9*16 cases.
Fifteen selected duals suffice; only seventeen distinct dual/higher-pattern
point inequalities are needed after selecting a certificate for each case.
Their minimum bound is 305, for counts (8,18) and higher pattern 1111.
No optimality claim is made for these bounds.

## Independent audit

`code/audit_compressed_tail_branch.py --pattern 001` passed using only
Python's standard library and integer arithmetic. It checks:

* 12495350 ordinary core pair-union/final-even-line combinations, without
  label symmetry, reproducing the nine-pair frontier;
* the compressed capacity table using explicit integers modulo 210 and
  all 128 allocations of the seven even shapes;
* 1004175 local plane-count tuples for the seventeen selected point
  inequalities, minimizing each of the eleven tail counts independently
  over all seventeen possibilities;
* every one of the 144 final count/higher-pattern inequalities.

The audit passed in about 14 seconds. It does not rely on solver
optimality, a floating-point tolerance, or the earlier 167-input screening
catalogue. All fifty regression tests also passed.

Reproduce the theorem:

```
python code/audit_compressed_tail_branch.py --pattern 001
```

Main files:

* `code/compressed_tail_branch.py`: count-frontier certificate generator,
  including resume support and reduction of the dual pool.
* `code/audit_compressed_tail_branch.py`: independent whole-branch audit.
* `data/compressed-tail-branch-001.json`: selected duals, assignments to
  all 144 cases, and their integer bounds.
* `data/audit-compressed-tail-branch-001.json`: independent audit report.

This is not proof-assistant verification and no novelty claim is made.

## Research state for continuation

The earlier joint-slice screening's 167 unresolved 001 inputs all now
have positive generated bounds, after extending its dual pool to sixteen
members. That work is superseded for the whole-branch theorem by the
simpler, independently audited count-only proof above. Do not spend time
finishing an orbit-completeness audit merely to prove 001 again.

The generator can also be applied to 010 and 011. A first 010 run is saved
in `data/compressed-tail-branch-010.json`: 99 of its 160 count/pattern
cells remain without positive bounds. The newly tried worst-case input,
counts (7,16) with higher pattern 1111, produced only the inconclusive
dual lower bound -10383. This does not prove that the relaxation cannot
work; solver termination information was not retained for that first
probe. Subsequent generated duals now preserve that diagnostic metadata.

If count-only bounds remain insufficient, the joint slice and exact
subset geometry tools are still available for 010/011. Even excluding
both branches would establish only the eight-prime obstruction, not the
global minimum-modulus-2 conjecture.
