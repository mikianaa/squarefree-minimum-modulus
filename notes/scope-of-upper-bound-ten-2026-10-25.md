# Historical: why the tested upper-bound comparison stopped at 10

**Superseded on 26 October 2026:** the same comparison now certifies minimum
modulus at most 7, with total loss <0.999087538091. Per-prime breakpoint
optimization and a longer distribution calculation overcome the numerical
gap discussed below. See `upper-bound-seven-2026-10-26.md`. The remaining
text records the evidence available on 25 October, not the current bound.

The saved evidence proves M_sf<=10. It does **not** prove that 10 is the
best possible outcome of every distortion, stop-loss, or moment argument.
The accurate statement is that 10 is the best certified cutoff reached with
the particular estimates and refinements tested here.

## The success criterion

Assume every modulus exceeds C0. Group congruences by their largest prime,
and let B_p be the set covered at the stage for p. The distortion construction
uses a probability measure whose earlier marginals are preserved. If the
sum of the upper bounds for these stage losses is less than 1, the family
cannot cover. Thus a schedule giving loss <1 proves M_sf<=C0.

The fully evaluated certificate `data/cert10-full.json` gives

    total loss upper bound = 0.922645168979... < 1.

To improve the conclusion below 10, it suffices and is necessary to target
M_sf<=7: the integers 8 and 9 are not squarefree. The hypothesis that all
moduli exceed 7, 8, or 9 is the same in this problem. The saved numerical
schedule `data/schedule8-float-only.json` gives 1.0178771888087672, above 1.
That value is a numerical schedule evaluation, not a globally certified
minimum over all schedules. More exact evaluation and the tested refinements
have not supplied a certificate below 1 at this cutoff.

An older progress-table figure 1.0078 was not supported by a saved cutoff-8
certificate found in the current repository. It has been removed from the
current status rather than treated as authoritative.

## Why allowing just modulus 10 matters

At the prime-5 stage, the allowed moduli dividing 2*3*5 are

* for minimum modulus >10: 15 and 30;
* for minimum modulus >7: 10, 15 and 30.

Earlier stages have no allowed moduli, so the distribution modulo 6 is uniform.
In the aligned comparison, let X and Y be independent Bernoulli variables of
parameters 1/2 and 1/3, indicating the chosen residue modulo 2 and modulo 3.
The comparison loads at prime 5 are respectively

    S_old = (Y+XY)/5,       S_new = (X+Y+XY)/5.

Their means are 1/10 and 1/5. With distortion delta=1/5, the stop-loss values
E[(S-delta)_+]/(1-delta) are respectively 1/24 and 1/12. This is an illustrative
fixed-parameter calculation, not an optimization of the complete argument.
The local loads are actually realizable: use aligned earlier residues and
distinct residues modulo 5 for the two or three congruences.

Increasing delta can reduce the immediate loss, but it also increases the
factor nu=1/(1-delta) in the bounds for future coordinate probabilities.
Those larger probabilities increase subsequent comparison losses. Thus
there is a tradeoff, rather than a free way to erase the extra prime-5 loss.

## What the sharpness experiments do and do not establish

At the tested parameters delta_5=.195 and delta_7=.1425, complete small-stage
enumeration gives the worst combined loss

    409327/1932805 = 0.2117787360856372...

The comparison is approximately 0.211786. Consequently, replacing those two
stage estimates by their exact worst case offers only a very small saving
at these parameters. At prime 11, a tested actual family under the chosen
prefix measure achieves about .05509, against a comparison value about
.05529. The attempted truncation refinement therefore cannot claim the much
larger saving predicted by simply truncating the comparison variable.

These examples establish local or specified-parameter limitations. They do
not establish that one family realizes every worst stage simultaneously, or
that every possible schedule or stronger comparison has total at least 1.
A separate rigorous all-schedule barrier exists only for the *factored
quadratic continuation* specified in
`notes/quadratic-tail-barrier-2026-10-25.md`; it must not be broadened to the
whole stop-loss framework.

## Relation to the current calculations

The present core and group models retain which actual residues can occur
together. They seek to exploit compatibility information discarded by the
stagewise comparison. This is a motivated change of approach, not a logically
forced consequence of an absolute optimality theorem for the bound 10.

No covering with minimum modulus 10, or even greater than 2, has been found.
Failure of an upper estimate to drop below 1 is not evidence of such a cover.
