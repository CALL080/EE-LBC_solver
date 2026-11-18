"""Data and time of the full-round linear attacks on EE-LBC (numbers quoted in Sect. 4 of the paper).

Run:  python code/complexity.py   (output kept in outputs/complexity.log)

Model: Selcuk, J. Cryptol. 21:131-147 (2008), Corollary 1, eq. (12):
    N = ((Phi^-1(P_S) + Phi^-1(1 - 2^(-a-1))) / 2)^2 * |p - 1/2|^-2,   |p - 1/2| = |c| / 2,
for an a-bit advantage on an m-bit subkey guess, assuming wrong keys give correlation 0.
|c| is the single-trail correlation 2^-L of the bit-21 trail over L S-layers; the restricted hull
(hull.c, outputs/hull_eelbc_lin_21.log) gives a slightly larger expected c^2.
"""
from math import log2
from statistics import NormalDist

Z = NormalDist().inv_cdf
ROUNDS, KEY = 31, 80

print('[1] Layers covered by the trail: 31 rounds = 2 rounds of guessed key + L layers + 1 keyless last round')
for guessed, m in ((0, 0), (1, 4), (2, 20)):
    L = ROUNDS - 1 - guessed
    print(f'  guess {guessed} round(s) ({m:2d} key bits): L = {L} layers, single-trail |c| = 2^-{L}, c^2 = 2^-{2 * L}')

print('\n[2] Key recovery (Algorithm 2) with 20 guessed key bits, L = 28, |c| = 2^-28 (Selcuk eq. 12)')
L, m = 28, 20
for a in (10, 16, 20):
    for ps in (0.5, 0.9):
        n = ((Z(ps) + Z(1 - 2 ** (-a - 1))) ** 2) * 2 ** (2 * L)
        t_list = 2 ** (m - a) * 2 ** (KEY - m)          # trial encryptions for the remaining key bits
        print(f'  a = {a:2d}, P_S = {ps}: N = 2^{log2(n):.2f} known plaintexts;'
              f' exhaustive part = 2^{log2(t_list):.0f} encryptions; counters 2^17; evaluation 2^{m + 16} table look-ups')

print('\n[3] Distinguisher on the full 31 rounds without key guess, L = 30, |c| = 2^-30')
# Test: accept "EE-LBC" if |empirical correlation| > t0 = |c|/2 (fixed threshold).  Z = sqrt(N) * empirical correlation,
# d = sqrt(N) |c|, t = d/2, f = N/2^64.  Normal approximation; the two sampling models must be used consistently:
#   distinct plaintexts (without replacement): cipher Z ~ N(+-d, 1-f); random permutation Z ~ N(0, 1), since the
#     codebook correlation of a random permutation has variance about 2^-64 and adds f to the sampling variance 1-f;
#   independent plaintexts (with replacement): cipher Z ~ N(+-d, 1); random permutation Z ~ N(0, 1+f).
# False positive (two-sided) 2 Phi(-t/sqrt(v0)); false negative Phi((t-d)/sqrt(v1)) - Phi((-t-d)/sqrt(v1)).
c = 2.0 ** -30
def errors(d, f):
    """(FP, FN) for distinct and for independent plaintexts, threshold d/2 on the scale of Z."""
    t, out = d / 2, {}
    for name, v0, v1 in (('distinct', 1.0, 1 - f), ('independent', 1 + f, 1.0)):
        fp = 2 * NormalDist().cdf(-t / v0 ** 0.5)
        fn = NormalDist().cdf((t - d) / v1 ** 0.5) - NormalDist().cdf((-t - d) / v1 ** 0.5)
        out[name] = (fp, fn)
    return out


for lgN in (60, 61, 62, 63):
    e = errors(c * (2 ** lgN) ** 0.5, 2.0 ** (lgN - 64))
    print(f'  N = 2^{lgN}: distinct plaintexts FP {e["distinct"][0]:.3f}, FN {e["distinct"][1]:.3f};'
          f' independent plaintexts FP {e["independent"][0]:.3f}, FN {e["independent"][1]:.3f}')

print('\n[4] Rounds needed for the best single trail to fall below 2^-64 (B_R = 2R for R >= 8, last round keyless)')
R = next(r for r in range(8, 64) if 2 * r > 64)
print(f'  smallest layer count with B_R > 64: {R}; plus the keyless last round: {R + 1} rounds (EE-LBC has {ROUNDS})')

print('\n[5] Spread of |c| over keys: normal model of Leander (EUROCRYPT 2011, Sect. 4.2)')
# c = eta * s + X, eta = 2^-L (the trail of Proposition 1), X ~ N(0, sigma^2) for all other trails, and
# sigma^2 = (hull sum) - eta^2, where the hull sum over masks with <= 2 active S-boxes per round is a LOWER
# bound on the true sum, so the model may understate the spread.  |c| / eta ~ |N(1, ratio)|.
import os, re
here = os.path.dirname(os.path.abspath(__file__))
hull = {}
for line in open(os.path.join(here, 'outputs', 'hull_eelbc_lin_21.log')):
    m = re.match(r'R=\s*(\d+)\s+log2\(restricted sum\) =\s*(\S+)', line)
    if m:
        hull[int(m.group(1))] = float(m.group(2))
ND = NormalDist()


def frac_above(t, ratio):            # P(|Z| >= t) for Z ~ N(1, ratio)
    sd = ratio ** 0.5
    return 1 - (ND.cdf((t - 1) / sd) - ND.cdf((-t - 1) / sd))


def quantile(q, ratio):              # t with P(|Z| < t) = q, by bisection
    lo, hi = 0.0, 10.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if 1 - frac_above(mid, ratio) < q:
            lo = mid
        else:
            hi = mid
    return lo


measured = {}
p = os.path.join(here, 'outputs', 'corr_16keys.log')
if os.path.exists(p):
    for line in open(p):
        m = re.match(r'lin\s+r=\s*(\d+).*corr=(\S+)', line)
        if m:
            measured.setdefault(int(m.group(1)) - 1, []).append(log2(abs(float(m.group(2)))))
for L in (8, 10, 12, 28, 30):
    ratio = 2 ** (hull[L] + 2 * L) - 1
    qs = [log2(quantile(q, ratio)) for q in (0.25, 0.5, 0.75)]
    line = (f'  L={L:2d}: sigma^2/eta^2 >= {ratio:.3f}; model log2|c| + L at Q1/median/Q3: '
            f'{qs[0]:+.2f} / {qs[1]:+.2f} / {qs[2]:+.2f}')
    if L in measured:
        v = sorted(measured[L])
        k = len(v)
        line += (f';  measured ({k} keys): {v[round(0.25 * (k - 1))] + L:+.2f} / {v[round(0.5 * (k - 1))] + L:+.2f} / '
                 f'{v[round(0.75 * (k - 1))] + L:+.2f}')
    else:
        line += (f';  fraction of keys with |c| >= 2^-L: {frac_above(1, ratio):.2f}, >= 2^-(L+0.5): '
                 f'{frac_above(2 ** -0.5, ratio):.2f}, >= 2^-(L+1): {frac_above(0.5, ratio):.2f}')
    print(line)

print('\n[6] Median complexities (Leander, EUROCRYPT 2011, Definition 1) from the fixed-key hull over 10000 keys')
# outputs/fixedhull_stats_window.log: distribution of log2|c| + L over the round keys that the attacks use:
# rounds 3..30 (off = 2) for the key recovery, rounds 1..30 (off = 0) for the distinguisher.
# fixedhull.c reproduces the measured correlations at L = 8..13 (summarize_experiments.log [5], check_keyrec_hull.log [1]).
med = {}
for line in open(os.path.join(here, 'outputs', 'fixedhull_stats_window.log')):
    m = re.match(r'stats L=\s*(\d+) off=(\d+) keys=\d+.*median ([+-][\d.]+)', line)
    if m:
        med[(int(m.group(1)), int(m.group(2)))] = float(m.group(3))
lc = -28 + med[(28, 2)]
print(f'  key recovery, rounds 3..30: median log2|c| = {lc:.2f}; N = 2^6/c^2 = 2^{6 - 2 * lc:.2f}, N = 2^7/c^2 = 2^{7 - 2 * lc:.2f}'
      ' (success measured for keys with N c^2 >= 2^6: check_keyrec_hull.log [2])')
ld = -30 + med[(30, 0)]
for lgN in (62, 63):
    e = errors(2 ** ld * 2 ** (lgN / 2), 2.0 ** (lgN - 64))
    print(f'  distinguisher, rounds 1..30: median log2|c| = {ld:.2f}; N = 2^{lgN}, threshold |c_med|/2 (FN at most this for half of the keys):'
          f' distinct plaintexts FP {e["distinct"][0]:.3f}, FN {e["distinct"][1]:.3f};'
          f' independent plaintexts FP {e["independent"][0]:.3f}, FN {e["independent"][1]:.3f}')
