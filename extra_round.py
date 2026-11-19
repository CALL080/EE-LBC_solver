"""The variant of Sect. 6 (Scope of the Attacks): guess also the four bits of K_31 after S-box 5 of round 30.

Run:  python3 extra_round.py   (output kept in outputs/extra_round.log)
[1] the outputs of S-box 5 in round 30 go to bits 5, 21, 37, 53 of x_30, so y = x_30 xor K_31 and four bits of K_31 give
    the S-box output nibble; S^-1 then gives bit 21 of x_29 xor K_30 (the input bit of the trail one round earlier).
[2] wrong guesses of these four bits: correlation between the predicted bit and the right one (no tie if never +-1).
[3] data factor: the trail loses one round (correlation x2, data / 4), and the guess grows from 20 to 24 bits, which
    raises N c^2 for the same success in the model of Selcuk (J. Cryptol. 2008, Eq. (12)) by the ratio printed here.
This variant is not evaluated experimentally; the paper only states the data factor.
"""
from math import log2
from statistics import NormalDist

S = [0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2]
SI = [S.index(v) for v in range(16)]


def P(j):
    return 16 * (j % 4) + j // 4


print('[1] bits of x_30 fed by the outputs of S-box 5 in round 30:', [P(j) for j in range(20, 24)])


def b1(v):
    return (SI[v] >> 1) & 1


rho = {d: sum(1 - 2 * (b1(v) ^ b1(v ^ d)) for v in range(16)) / 16 for d in range(1, 16)}
print('[2] correlation of the predicted bit for the 15 wrong guesses:', sorted(set(abs(r) for r in rho.values())),
      '; equal to +-1 for', sum(abs(r) == 1 for r in rho.values()), 'guesses')
z = NormalDist().inv_cdf
f = {a: (z(0.9) + z(1 - 2 ** -(a + 1))) ** 2 for a in (20, 24)}
ratio = f[24] / f[20]
print(f'[3] Selcuk factor for rank 1 with probability 0.9: a = 20: {f[20]:.2f}, a = 24: {f[24]:.2f}, ratio 2^{log2(ratio):.2f};'
      f' data factor 4 / {ratio:.3f} = {4 / ratio:.2f} (2^{2 - log2(ratio):.2f})')
