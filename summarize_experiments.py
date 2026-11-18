"""Summaries of outputs/experiments.log, keyrec_*.log and corr_16keys.log (Tables 3 and 5).

Run:  python code/summarize_experiments.py   (output kept in outputs/summarize_experiments.log)
"""
import math
import os
import re
from collections import defaultdict

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')

print('[1] Linear: p[21] xor y[21] over r rounds (r-1 layers); per key log2|c|, and log2 of the mean c^2')
lin = defaultdict(list)
diff = defaultdict(list)
for line in open(os.path.join(D, 'experiments.log')):
    m = re.match(r'lin\s+r=\s*(\d+).*N=2\^(\d+).*(ks|indep)\s+corr=(\S+)', line)
    if m:
        lin[int(m.group(1))].append((m.group(3), float(m.group(4)), int(m.group(2))))
    m = re.match(r'diff\s+r=\s*(\d+).*N=2\^(\d+).*hits=(\d+)', line)
    if m:
        diff[int(m.group(1))].append((int(m.group(3)), int(m.group(2))))
for r in sorted(lin):
    v = lin[r]
    per = ' '.join(f'{"i" if k == "indep" else ""}{math.log2(abs(c)):.2f}' for k, c, _ in v)
    mean = math.log2(sum(c * c for _, c, _ in v) / len(v)) / 2
    lo, hi = min(math.log2(abs(c)) for _, c, _ in v), max(math.log2(abs(c)) for _, c, _ in v)
    print(f'  r={r:2d} layers={r - 1:2d} N=2^{v[0][2]} predicted {-(r - 1):4d}  measured [{lo:.2f}, {hi:.2f}]'
          f'  sqrt(mean c^2)=2^{mean:.2f}  ({len(v)} keys; i = independent round keys): {per}')

print('\n[2] Differential: (p, p xor 1) -> y-difference 1 over r-1 layers; predicted 2^(-3(r-1)) for the 1-bit trail')
for r in sorted(diff):
    v = diff[r]
    tot, n = sum(h for h, _ in v), sum(2 ** lg for _, lg in v)
    print(f'  r={r} layers={r - 1} predicted 2^{-3 * (r - 1)}  measured 2^{math.log2(tot / n):.2f}  ({tot} hits in {len(v)} keys)')

print('\n[3] Key recovery (20-bit guess; r rounds = 2 guessed + L = r-3 trail rounds + 1 keyless): rank of the correct key')
for fn in ('keyrec_r12_32seeds.log', 'keyrec_longer.log'):
    p = os.path.join(D, fn)
    if not os.path.exists(p):
        continue
    by = defaultdict(list)
    for line in open(p):
        m = re.search(r'keyrec r=\s*(\d+) trail_layers=\s*(\d+) N=2\^(\d+).*rank of true 20-bit key = (\d+)', line)
        if m:
            r, L, lg, rank = map(int, m.groups())
            by[(r, L, lg)].append(rank)
    for (r, L, lg) in sorted(by):
        ranks = by[(r, L, lg)]
        print(f'  {fn}: r={r} L={L} N=2^{lg} (= 2^{lg - 2 * L}/c^2): rank 1 in {sum(x == 1 for x in ranks)}/{len(ranks)},'
              f' rank <= 2^10 in {sum(x <= 1024 for x in ranks)}/{len(ranks)}, ranks {sorted(ranks)}')

print('\n[4] Spread of the correlation over 16 keys (key schedule): quartiles of log2|c| against the prediction -L')
p = os.path.join(D, 'corr_16keys.log')
if os.path.exists(p):
    by = defaultdict(list)
    for line in open(p):
        m = re.match(r'lin\s+r=\s*(\d+).*corr=(\S+)', line)
        if m:
            by[int(m.group(1)) - 1].append(math.log2(abs(float(m.group(2)))))
    for L in sorted(by):
        v = sorted(by[L])
        q = lambda f: v[min(len(v) - 1, int(f * (len(v) - 1) + 0.5))]
        print(f'  L={L:2d} ({len(v)} keys): min {v[0]:.2f}  Q1 {q(0.25):.2f}  median {q(0.5):.2f}  Q3 {q(0.75):.2f}'
              f'  max {v[-1]:.2f}   predicted {-L}')

print('\n[5] Fixed-key hull (fixedhull.c, masks with <= 2 active S-boxes) against the measured correlations of the same keys')
pf, pm = os.path.join(D, 'fixedhull_check.log'), os.path.join(D, 'corr_16keys.log')
if os.path.exists(pf) and os.path.exists(pm):
    pred, meas = {}, {}
    for line in open(pf):
        m = re.match(r'check L=\s*(\d+)(?: off=0)? seed=(\d+)\s+predicted corr=(\S+)', line)   # off=0: rounds 1..L
        if m:
            pred[(int(m.group(1)), int(m.group(2)))] = float(m.group(3))
    for line in open(pm):
        m = re.match(r'lin\s+r=\s*(\d+).*seed=(\d+).*corr=(\S+)', line)
        if m:
            meas[(int(m.group(1)) - 1, int(m.group(2)))] = float(m.group(3))
    for L in sorted({k[0] for k in pred}):
        keys = [k for k in pred if k[0] == L and k in meas]
        same = sum((pred[k] > 0) == (meas[k] > 0) for k in keys)
        dev = sorted(abs(math.log2(abs(meas[k]) / abs(pred[k]))) for k in keys)
        print(f'  L={L:2d}: {len(keys)} keys, same sign {same}/{len(keys)}, |log2(measured/predicted)| median {dev[len(dev) // 2]:.2f}, max {dev[-1]:.2f}')
p = os.path.join(D, 'fixedhull_stats.log')
if os.path.exists(p):
    print('  distribution over 2000 keys (log2|c| + L):')
    for line in open(p):
        print('   ', line.strip())
