"""EE-LBC (Gupta and Sasikala, Int. J. Electr. Comput. Eng. Syst. 16(7):507-515, 2025), reference implementation.

Every line below cites the proposal (source/EE-LBC_proposal.pdf; page numbers are the journal's, 507-515).
Where the proposal is ambiguous the choice is a named option; the default is our reading of the text.

  block 64 bits, key 80 bits, 31 rounds                               (p. 508 bullets, p. 509 Sect. 3)
  round i = 1..31:  B <- P( S16( B xor K_i ) )                       (p. 510 "three phases 31 times")
  no key is added after the last round                               (p. 509 "One key K_i ... 1<=i<=31"; Fig. 1)
  S = {C,5,6,B,7,0,A,D,1,E,F,8,4,9,3,2} on every nibble              (p. 510 Substitution Layer)
  P : bit i -> 16*(i mod 4) + floor(i/4), bit 0 = LSB                (p. 510 Table 1, "ith bit is moved to the
                                                                      position indicated by the value inside the cell")
  key schedule on the 80-bit register (K79..K0):                     (p. 509-510 Key Scheduler)
      K_1 = Key64 = (K79..K16)
      for i = 2..31:  rotate left by 61; Key64 through the substitution layer;
                      (K19..K15) ^= counter;  K_i = (K79..K16)

Options (the points the proposal leaves open, see notes/spec_ja.md):
  perm      'table'  : bit i -> PERM[i] as the text says (default)
            'inverse': the cell index is the new position (the figure's label "NewBit Position")
  ks_sbox   'all'    : "Key64 goes through the substitution layer" = all 16 nibbles of K79..K16 (default, literal)
            'top'    : only the nibble K79..K76 (as in PRESENT-80)
  counter   'i'      : the index i of the round key being generated, i = 2..31 (default, "counter of the current round")
            'i-1'    : 1..30 (as in PRESENT-80's round_counter)
  ks_shift  'rotate' : 80-bit rotation by 61 (default)
            'shift'  : logical shift of the 80-bit register (literal "<< 61"; destroys the key, kept only for checks)
  final_key False    : no key after round 31 (default, as stated)
            True     : a 32nd key is generated and XORed after round 31 (as in PRESENT; not in the proposal)
"""

SBOX = [0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2]
SBOX_INV = [SBOX.index(y) for y in range(16)]
PRESENT_SBOX = [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2]

# Table 1 read row-wise: cell number i (0-based) holds the new position of input bit i.
PERM = [16 * (i % 4) + i // 4 for i in range(64)]
PERM_INV = [PERM.index(j) for j in range(64)]
assert sorted(SBOX) == list(range(16)) and sorted(PERM) == list(range(64))

ROUNDS = 31
M64 = (1 << 64) - 1
M80 = (1 << 80) - 1

DEFAULT = dict(perm='table', ks_sbox='all', counter='i', ks_shift='rotate', final_key=False, sbox=SBOX)


def sub16(x, box=SBOX):
    """Apply the 4-bit S-box to each of the 16 nibbles of a 64-bit word."""
    y = 0
    for j in range(16):
        y |= box[(x >> (4 * j)) & 0xF] << (4 * j)
    return y


def perm64(x, table):
    """Bit i of x goes to bit table[i] of the result (bit 0 = LSB)."""
    y = 0
    for i in range(64):
        y |= ((x >> i) & 1) << table[i]
    return y


def _tables(opt):
    if opt['perm'] == 'table':
        return PERM, PERM_INV
    if opt['perm'] == 'inverse':
        return PERM_INV, PERM
    raise ValueError(opt['perm'])


def round_keys(key, **kw):
    """Round keys K_1..K_31 (and K_32 if final_key) of an 80-bit key."""
    opt = {**DEFAULT, **kw}
    assert 0 <= key <= M80
    k = key
    keys = [k >> 16]
    n = ROUNDS + (1 if opt['final_key'] else 0)
    for i in range(2, n + 1):
        if opt['ks_shift'] == 'rotate':
            k = ((k << 61) | (k >> 19)) & M80
        elif opt['ks_shift'] == 'shift':
            k = (k << 61) & M80
        else:
            raise ValueError(opt['ks_shift'])
        if opt['ks_sbox'] == 'all':
            k = (sub16(k >> 16, opt['sbox']) << 16) | (k & 0xFFFF)
        elif opt['ks_sbox'] == 'top':
            k = (opt['sbox'][k >> 76] << 76) | (k & ((1 << 76) - 1))
        else:
            raise ValueError(opt['ks_sbox'])
        c = i if opt['counter'] == 'i' else i - 1
        k ^= (c & 0x1F) << 15
        keys.append(k >> 16)
    return keys


def encrypt(pt, key, rounds=ROUNDS, **kw):
    opt = {**DEFAULT, **kw}
    p, _ = _tables(opt)
    rk = round_keys(key, **opt)
    x = pt
    for r in range(rounds):
        x = perm64(sub16(x ^ rk[r], opt['sbox']), p)
    if opt['final_key'] and rounds == ROUNDS:
        x ^= rk[ROUNDS]
    return x


def decrypt(ct, key, rounds=ROUNDS, **kw):
    """The inverse of encrypt.  The proposal's Sect. 3B keeps the same substitution layer for decryption;
    the S-box is not an involution, so that literal reading does not invert encryption (see check_spec.py)."""
    opt = {**DEFAULT, **kw}
    _, pinv = _tables(opt)
    rk = round_keys(key, **opt)
    x = ct
    if opt['final_key'] and rounds == ROUNDS:
        x ^= rk[ROUNDS]
    for r in reversed(range(rounds)):
        x = sub16(perm64(x, pinv), [opt['sbox'].index(y) for y in range(16)]) ^ rk[r]
    return x


if __name__ == '__main__':
    import random
    rng = random.Random(1)
    for opt in [DEFAULT, dict(perm='inverse', ks_sbox='top', counter='i-1', final_key=True)]:
        for _ in range(200):
            pt, key = rng.getrandbits(64), rng.getrandbits(80)
            assert decrypt(encrypt(pt, key, **opt), key, **opt) == pt
    print('decrypt(encrypt(x)) == x for 400 random (plaintext, key) pairs')
    for key in (0, M80):
        print(f'key={key:020x}  E(0)={encrypt(0, key):016x}  E(ffff..)={encrypt(M64, key):016x}')
