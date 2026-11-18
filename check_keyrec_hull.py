"""Key-recovery runs against the fixed-key hull of the same keys (Sect. 4.3).

Run:  python3 check_keyrec_hull.py   (output kept in outputs/check_keyrec_hull.log)
Inputs (outputs/):
  keyrec_r12_32seeds.log, keyrec_longer.log   rank of the correct 20-bit guess and the measured |c| of the correct guess
  fixedhull_check_keyrec.log                  fixedhull.c check with off = 2 (rounds 3 .. r-1) for the same seeds
[1] measured |c| against the hull prediction (the measurement has a relative noise of about 1/sqrt(N c^2));
[2] success restricted to keys whose predicted |c| is at least the nominal 2^-L, i.e. N c^2 >= N c_0^2.
"""
import math
import os
import re
from collections import defaultdict

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
pred = {}
for line in open(os.path.join(D, 'fixedhull_check_keyrec.log')):
    m = re.match(r'check L=\s*(\d+) off=2 seed=(\d+)\s+predicted corr=(\S+)', line)
    if m:
        pred[(int(m.group(1)), int(m.group(2)))] = abs(float(m.group(3)))
runs = []
for fn in ('keyrec_r12_32seeds.log', 'keyrec_longer.log'):
    for line in open(os.path.join(D, fn)):
        m = re.search(r'keyrec r=\s*(\d+) trail_layers=\s*(\d+) N=2\^(\d+) seed=(\d+)\s+\|corr\(true key\)\|=2\^(\S+) .*'
                      r'rank of true 20-bit key = (\d+)', line)
        if m:
            r, L, lg, seed = map(int, m.groups()[:4])
            runs.append(dict(r=r, L=L, lg=lg, seed=seed, meas=float(m.group(5)), rank=int(m.group(6))))

print('[1] measured log2|c| of the correct guess against the fixed-key hull (rounds 3 .. r-1), per key and data size')
by = defaultdict(list)
for x in runs:
    p = pred.get((x['L'], x['seed']))
    if p:
        by[(x['L'], x['lg'])].append(x['meas'] - math.log2(p))
for k in sorted(by):
    v = sorted(abs(d) for d in by[k])
    print(f'  L={k[0]:2d} N=2^{k[1]}: {len(v)} keys, |log2(measured/predicted)| median {v[len(v) // 2]:.2f}, max {v[-1]:.2f}')

print('\n[2] success for keys with predicted |c| >= 2^-L (N c^2 >= N c_0^2), against all keys')
grp = defaultdict(lambda: [0, 0, 0, 0, 0, 0])
for x in runs:
    p = pred.get((x['L'], x['seed']))
    nc = x['lg'] - 2 * x['L']
    g = grp[nc]
    g[0] += 1; g[1] += x['rank'] == 1; g[2] += x['rank'] <= 1024
    if p and math.log2(p) >= -x['L']:
        g[3] += 1; g[4] += x['rank'] == 1; g[5] += x['rank'] <= 1024
for nc in sorted(grp):
    g = grp[nc]
    print(f'  N c_0^2 = 2^{nc}: all keys: rank 1 {g[1]}/{g[0]}, top 2^10 {g[2]}/{g[0]};'
          f'  keys with |c| >= c_0: rank 1 {g[4]}/{g[3]}, top 2^10 {g[5]}/{g[3]}')
