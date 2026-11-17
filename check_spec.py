"""Facts about the EE-LBC specification that need no search (S-box tables, permutation, the paper's examples).

Run:  python code/check_spec.py   (output kept in outputs/check_spec.log)
Each block prints a heading [n] that notes/spec_ja.md and README.md refer to.
"""
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eelbc import SBOX, SBOX_INV, PRESENT_SBOX, PERM, encrypt, decrypt, sub16, perm64, M80


def ddt(s):
    t = [[0] * 16 for _ in range(16)]
    for x in range(16):
        for a in range(16):
            t[a][s[x] ^ s[x ^ a]] += 1
    return t


def lat(s):
    """LAT[a][b] = #{x : a.x = b.S(x)} - 8  (bias times 16)."""
    t = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for b in range(16):
            n = sum(1 for x in range(16) if (bin(a & x).count('1') + bin(b & s[x]).count('1')) % 2 == 0)
            t[a][b] = n - 8
    return t


def show(t, name):
    print(f'{name} (rows: input mask/difference, columns: output)')
    print('     ' + ' '.join(f'{b:>2x}' for b in range(16)))
    for a in range(16):
        print(f'  {a:x}: ' + ' '.join(f'{v:>2d}' for v in t[a]))


def degree(s):
    """Algebraic degree of each output bit (ANF via Moebius transform)."""
    degs = []
    for j in range(4):
        f = [(s[x] >> j) & 1 for x in range(16)]
        for i in range(4):
            for x in range(16):
                if x >> i & 1:
                    f[x] ^= f[x ^ (1 << i)]
        degs.append(max((bin(x).count('1') for x in range(16) if f[x]), default=0))
    return degs


W1 = (1, 2, 4, 8)

print('[1] S-box')
print('  EE-LBC  S =', ' '.join(f'{v:X}' for v in SBOX))
print('  PRESENT S =', ' '.join(f'{v:X}' for v in PRESENT_SBOX))
diff = [x for x in range(16) if SBOX[x] != PRESENT_SBOX[x]]
print('  entries that differ from PRESENT:', ', '.join(f'S({x:X})={SBOX[x]:X} (PRESENT {PRESENT_SBOX[x]:X})' for x in diff))
print('  bijective:', sorted(SBOX) == list(range(16)), ' involution:', all(SBOX[SBOX[x]] == x for x in range(16)))
print('  fixed points:', [x for x in range(16) if SBOX[x] == x], ' algebraic degree per output bit (LSB first):', degree(SBOX))

for name, s in (('EE-LBC', SBOX), ('PRESENT', PRESENT_SBOX)):
    D, L = ddt(s), lat(s)
    mx = max(D[a][b] for a in range(1, 16) for b in range(16))
    ml = max(abs(L[a][b]) for a in range(16) for b in range(1, 16))
    nmax = sum(1 for a in range(1, 16) for b in range(16) if D[a][b] == mx)
    nlin = sum(1 for a in range(16) for b in range(1, 16) if abs(L[a][b]) == ml)
    print(f'\n[2] {name}: differential uniformity {mx} (max probability {mx}/16), entries equal to it: {nmax}')
    print(f'    {name}: max |LAT| {ml} (max bias {ml}/16), entries equal to it: {nlin}')
    one = [(a, b, D[a][b]) for a in W1 for b in W1 if D[a][b]]
    print(f'    {name}: single-bit -> single-bit differential transitions (a, b, DDT):', one or 'none')
    onel = [(a, b, L[a][b]) for a in W1 for b in W1 if L[a][b]]
    print(f'    {name}: single-bit -> single-bit linear approximations (a, b, LAT):', onel or 'none')
    if name == 'EE-LBC':
        show(D, '    DDT')
        show(L, '    LAT')

print('\n[3] Permutation (Table 1)')
print('  P(i) = 16*(i mod 4) + floor(i/4):', PERM == [16 * (i % 4) + i // 4 for i in range(64)])
print('  same as PRESENT pLayer P(i) = 16 i mod 63 (i<63), P(63)=63:',
      PERM == [(16 * i) % 63 if i < 63 else 63 for i in range(64)])
print('  fixed points:', [i for i in range(64) if PERM[i] == i])
print('  order of P:', next(k for k in range(1, 10) if all(_ == i for i, _ in enumerate(
    __import__('functools').reduce(lambda acc, _: [PERM[v] for v in acc], range(k), list(range(64)))))))

print('\n[4] The example in Sect. 4 (Differential Cryptanalysis)')
P1, P2 = 0x4172756e61204b47, 0x4172756e61204b48
C1, C2 = 0xfff2f964d2012604, 0x5a3c5405043f0d45
print(f'  P1 xor P2 = {P1 ^ P2:016x}  ({bin(P1 ^ P2).count("1")} bits; the text says "just one-bit difference")')
print(f'  C1 xor C2 = {C1 ^ C2:016x}  ({bin(C1 ^ C2).count("1")} bits)')
print('  P1 as ASCII:', repr(P1.to_bytes(8, 'big').decode('ascii')), ' P2:', repr(P2.to_bytes(8, 'big').decode('ascii')))
print('  no key is given, so these are not test vectors; trying a few obvious keys under every reading:')


def asc(s):
    b = s.encode('ascii')[:10].ljust(10, b'\0')
    return int.from_bytes(b, 'big')


keys = {'0': 0, '1..1': M80, '0123456789abcdef0123': 0x0123456789abcdef0123,
        'ascii "Aruna Gupt"': asc('Aruna Gupta'), 'ascii "Aruna KG"': asc('Aruna KG'),
        'ascii "1234567890"': asc('1234567890'), 'ascii "abcdefghij"': asc('abcdefghij'),
        'ascii "0000000000"': asc('0000000000'), 'ascii "EE-LBC"': asc('EE-LBC'),
        'ascii "Sathyabama"': asc('Sathyabama'), 'ascii "secretkey!"': asc('secretkey!')}
options = [dict(perm=p, ks_sbox=k, counter=c, ks_shift=sh, final_key=f)
           for p, k, c, sh, f in itertools.product(('table', 'inverse'), ('all', 'top'), ('i', 'i-1'),
                                                    ('rotate', 'shift'), (False, True))]


def swap8(x):
    return int.from_bytes(x.to_bytes(8, 'big'), 'little')


hits = 0
for kname, key in keys.items():
    for opt in options:
        for endian in ('big', 'little'):
            p1 = P1 if endian == 'big' else swap8(P1)
            c = encrypt(p1, key, **opt)
            c = c if endian == 'big' else swap8(c)
            if c in (C1, C2):
                hits += 1
                print('   MATCH', kname, opt, endian)
print(f'  tried {len(keys)} keys x {len(options)} readings x 2 byte orders: {hits} matches')

print('\n[5] Decryption as written in Sect. 3B (same substitution layer, inverse permutation, reversed keys)')
rng = random.Random(5)
from eelbc import round_keys, PERM_INV
bad = 0
for _ in range(100):
    pt, key = rng.getrandbits(64), rng.getrandbits(80)
    ct = encrypt(pt, key)
    rk = round_keys(key)
    x = ct
    for r in reversed(range(31)):  # literal: inverse permutation, the same S-box, add the round key
        x = sub16(perm64(x, PERM_INV)) ^ rk[r]
    bad += (x != pt)
print(f'  literal reading fails to decrypt in {bad}/100 random cases; with the inverse S-box it succeeds in',
      sum(decrypt(encrypt(p, k), k) == p for p, k in ((rng.getrandbits(64), rng.getrandbits(80)) for _ in range(100))), '/100')

print('\n[6] Last round: no key follows the last S-layer and permutation')
key, pt = rng.getrandbits(80), rng.getrandbits(64)
rk = round_keys(key)
y = sub16(perm64(encrypt(pt, key), PERM_INV), SBOX_INV)  # anyone can strip P and S of round 31
print('  S^-1(P^-1(C)) == state after 30 rounds xor K_31:', y == encrypt(pt, key, rounds=30) ^ rk[30])

print('\n[7] One-active-S-box-per-round transitions at the fixed points of P (checked on the implementation)')
# Differential: S-box 0, input 1 -> output 1 (DDT 2), and P(0) = 0.  Linear: S-box 5, mask 2 -> 2 (LAT 4), P(21) = 21.
R, N = 3, 1 << 16
key = rng.getrandbits(80)
hit = sum(encrypt(x, key, rounds=R) ^ encrypt(x ^ 1, key, rounds=R) == 1 for x in (rng.getrandbits(64) for _ in range(N)))
print(f'  differential e_0 -> e_0 over {R} rounds: {hit}/{N} pairs, expected (2/16)^{R} * N = {N * (2 / 16) ** R:.0f}')
agree = sum(((x >> 21) & 1) == ((encrypt(x, key, rounds=R) >> 21) & 1) for x in (rng.getrandbits(64) for _ in range(N)))
print(f'  linear  e_21 -> e_21 over {R} rounds: |correlation| = {abs(2 * agree / N - 1):.4f}, '
      f'expected (4/8)^{R} = {(4 / 8) ** R:.4f} (sign depends on the key)')

print('\n[8] The four S-box conditions of PRESENT (Bogdanov et al., CHES 2007, Sect. 4), Walsh value = 2 * LAT')
for name, s in (('PRESENT', PRESENT_SBOX), ('EE-LBC', SBOX)):
    D, L = ddt(s), lat(s)
    c1 = max(D[a][b] for a in range(1, 16) for b in range(1, 16)) <= 4
    v2 = [(a, b) for a in W1 for b in W1 if D[a][b]]
    c3 = max(abs(2 * L[a][b]) for a in range(1, 16) for b in range(1, 16)) <= 8
    # Condition 4 is printed as 'S_b^W(a) = +-4', which PRESENT's own S-box does not meet (some values are 0);
    # its App. IV uses 'bias of any single-bit approximation <= 2^-3', i.e. |Walsh| <= 4.  We test that form.
    v4 = [(a, b, 2 * L[a][b]) for a in W1 for b in W1 if abs(2 * L[a][b]) > 4]
    print(f'  {name}: (1) DDT <= 4: {c1};  (2) no 1-bit -> 1-bit differential: {not v2} {v2 or ""};'
          f'  (3) |Walsh| <= 8: {c3};  (4*) 1-bit |Walsh| <= 4: {not v4} {v4 or ""}')

print('\n[9] With the PRESENT S-box, the PRESENT readings of (ii)-(iii) and a final key, EE-LBC is PRESENT-80')
present = dict(sbox=PRESENT_SBOX, ks_sbox='top', counter='i-1', final_key=True)
vectors = [(0, 0, 0x5579C1387B228445), (0, M80, 0xE72C46C0F5945049),
           ((1 << 64) - 1, 0, 0xA112FFC72F68417B), ((1 << 64) - 1, M80, 0x3333DCD3213210D2)]  # PRESENT, App. I
ok = all(encrypt(p, k, **present) == c for p, k, c in vectors)
print('  all four test vectors of PRESENT-80 (Bogdanov et al., App. I) reproduced:', ok)
print('  EE-LBC itself unchanged: E_0(0) =', f'{encrypt(0, 0):016x}', '(expected 1b7b51a05ea03b5c)')
