"""Monte Carlo of the independent comparison model for single stages (agent sqaudit).

Stage k with prime p: earlier primes q present independently with probability nu_q/q,
N = number present, S = 2^N / p (valid for p > C0), loss = E (S - delta)_+ / (1 - delta).
Compares with the certificate rows.  Usage: venv python mc_stage.py opt10.json 23 101
"""
import json, sys
import numpy as np

L = 10 ** 6


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


sch = json.load(open(sys.argv[1]))
edges, Ds = sch["edges"], [int(round(float(d) * L)) for d in sch["params"]]
C0 = sch["C0"]


def D_for(k):
    for e, D in zip(edges, Ds):
        if k <= e:
            return D


primes = sieve(200)
rng = np.random.default_rng(20261023)
nsamp = 20_000_000
for p in [int(a) for a in sys.argv[2:]]:
    k = primes.index(p) + 1
    assert p > C0
    delta = D_for(k) / L
    probs = np.array([L / ((L - D_for(j + 1)) * primes[j]) for j in range(k - 1)])
    # sample N in chunks
    tot = 0.0
    tot2 = 0.0
    done = 0
    chunk = 2_000_000
    while done < nsamp:
        U = rng.random((chunk, k - 1))
        N = (U < probs).sum(axis=1)
        S = 2.0 ** N / p
        x = np.maximum(S - delta, 0.0) / (1 - delta)
        tot += x.sum()
        tot2 += (x * x).sum()
        done += chunk
    mean = tot / done
    se = ((tot2 / done - mean * mean) / done) ** 0.5
    ES = np.prod(1 + probs) / p
    print("p=%d k=%d delta=%.6f  MC stop-loss = %.7f +- %.7f (1 s.e.)   exact E S = %.7f" % (p, k, delta, mean, se, ES))
