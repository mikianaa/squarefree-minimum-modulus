# A limited barrier for a quadratic continuation

This is a limitation of one specified estimate, not a proof that all moment
methods fail. It is separate from the seven-prime obstruction.

Suppose the initial prefix has squarefree period P and probability measure mu.
Write c_d=max_a mu(a mod d) and

    M = sum_{d|P} 3^omega(d) c_d.

Since c_d>=1/d, every initial measure, even one not restricted to the residual,
satisfies M>=product_{p|P}(1+3/p).

Consider using only the conventional quadratic stage cost

    M / (4 delta (1-delta) q^2)

and the factored moment update M -> M*(1+3/((1-delta)q)) at each later prime q.
If K is the optimal subsequent cost per unit moment, the optimal current cost
per unit moment is

    F_q(K) = K + (1+sqrt(1+12qK))^2/(4q^2).

Indeed, minimize A/(delta*(1-delta))+B/(1-delta)+K with A=1/(4q^2)
and B=3K/q. The minimizing delta is 1/(1+sqrt(1+12qK)), which is <=1/2.

Start with K=0 at a finite terminal prime and apply F backwards. This is a
lower bound on the value of the infinite-tail estimate, uniformly over all
distortion schedules. Integer square roots and divisions rounded DOWN give
rigorous rational lower bounds.

| Prefix through | Finite tail through | Lower bound on the best possible estimate |
| --- | --- | ---: |
| 7 | 109 | >1.00196 |
| 11 | 241 | >1.00229 |
| 13 | 709 | >1.00022 |
| 17 | 3943 | >1.00005 |

Thus this particular quadratic catalogue estimate cannot certify loss <1
starting at any of these prefixes, regardless of the chosen initial measure.
It does not rule out first-moment bounds, stop-loss bounds, future anchor
normalization, compatible-intersection estimates, or a later starting prefix.

Reproduce the exact fractions with `python code/quadratic_tail_barrier.py`.
They are saved in `data/quadratic-tail-barrier.json`.
