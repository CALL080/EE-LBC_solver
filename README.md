# Code for "One Fixed Point Suffices: Full-Round Linear Cryptanalysis of the Lightweight Block Cipher EE-LBC"

This repository reproduces every number, table, and figure of the paper.
EE-LBC is the 64-bit block cipher with an 80-bit key and 31 rounds of A. Gupta and T. Sasikala,
"Lightweight Block Cipher for Security in Resource-Constrained Network",
Int. J. Electr. Comput. Eng. Syst. 16(7):507-515, 2025, DOI 10.32985/ijeces.16.7.2 (called *the proposal* below;
the PDF of the proposal is not included).

## What the paper shows

1. The bit permutation of EE-LBC is that of PRESENT, and its S-box differs from the PRESENT S-box in four entries.
   These entries violate all four S-box conditions of PRESENT.
2. For every round key, bit 21 (a fixed point of the permutation) gives a one-round linear approximation with one active
   S-box and correlation +-1/2. By branch-and-bound, the minimum linear weight is exactly `2R` for 8 <= R <= 34.
3. No key is added after the last round, so the trail has to cover only 30 of the 31 rounds. For half of the keys,
   the full key is recovered with 2^62.3 known plaintexts and 2^70 encryptions (success 0.97), or with 2^63.3 known
   plaintexts and 2^60 encryptions (success 0.94). The proposal estimates 2^124 plaintexts.
4. The claims of the proposal (maximum differential probability 2^-2, at least 10 active S-boxes over 5 rounds and
   62 over 25 rounds) do not hold: the correct values are 2^-1, 5, and 25.

## Quick start

Requirements: Python 3.8 or later (standard library only; we used 3.12) and gcc with OpenMP (we used gcc 13.3 on
Ubuntu 24.04 under WSL).

```bash
python3 eelbc.py          # reference implementation, self-test
python3 check_spec.py     # facts that need no search (S-box, permutation, PRESENT conditions, Proposition 1)
./run_all.sh              # builds the C programs and regenerates outputs/ (about 2 minutes)
./run_all.sh --full       # also reruns the reduced-round experiments (about two hours more on 20 threads)
```

`outputs/` contains the logs used for the paper, and `outputs/figs/` the TeX inputs of Fig. 2 and Tables 3 and 5
generated from them. A rerun overwrites them; the searches print their running time and node count, which depend on
the machine, while the weights and all other numbers do not. The random seeds of the experiments are fixed, so the
measurements are the same on every run.

## Conventions

- Bit 0 is the least significant bit, as in the proposal. A 64-bit value is written in hexadecimal with bit 63 first.
- Round `i` (1 <= i <= 31) maps `x` to `P(S(x xor K_i))`; no key is added after the last round.
- `P(j) = 16 * (j mod 4) + floor(j / 4)` (Table 1 of the proposal). The S-box is `C 5 6 B 7 0 A D 1 E F 8 4 9 3 2`.
- The points that the proposal leaves open are options of `eelbc.py` (see its header):
  `perm` (direction of Table 1 of the proposal), `ks_sbox` (S-box on all of Key64 or on its top nibble), `counter`
  (`i` or `i-1`), `ks_shift` (rotation or logical shift), `final_key` (no key after round 31, or a 32nd key as in
  PRESENT), `sbox` (the S-box; `PRESENT_SBOX` gives PRESENT). The defaults are the readings of Sect. 2 of the paper.
  The results of Sects. 3 and 5 depend only on the S-box and `P`, and the attacks hold for every reading (Sect. 4.4).
- The weight of a trail is `-log2` of the product of the DDT probabilities (differential) or of the squared
  correlations (linear) of its active S-boxes. An `R`-round trail passes through `R` S-layers.
- In the experiments, `r` rounds of the cipher leave `L = r - 1` trail rounds (the last round is removed without a key),
  and the key recovery on `r` rounds leaves `L = r - 3` (two rounds of guessed key and the keyless last round).

## Claims of the paper and how to check them

| Paper | Claim | Command | Expected output (file in `outputs/`) |
|---|---|---|---|
| Sect. 2 | Encryption and decryption are inverse (400 random pairs); `E_0(0) = 1b7b51a05ea03b5c` | `python3 eelbc.py` | `eelbc_selftest.log` |
| Sect. 2.2 | The S-box differs from the PRESENT S-box at inputs 4, 8, D, E | `python3 check_spec.py` | `check_spec.log`, `[1]` |
| Sect. 2.2 | `P` is the PRESENT permutation; fixed points 0, 21, 42, 63 | `python3 check_spec.py` | `[3]` |
| Sect. 2.3 (vi) | Decryption with the same S-box layer fails (100 of 100) | `python3 check_spec.py` | `[5]` |
| Sect. 2.3 (vii) | The example plaintexts differ in four bits; no ciphertext reproduced for 11 keys x 32 readings x 2 byte orders | `python3 check_spec.py` | `[4]` |
| Sect. 3, Table 1 | DDT maximum 8 (`F -> F`), largest absolute correlation 3/4 | `python3 check_spec.py` | `[2]` |
| Sect. 3, Table 1 | PRESENT meets the four S-box conditions of PRESENT; EE-LBC violates all four | `python3 check_spec.py` | `[8]` |
| Prop. 1 | Over three rounds: 133 of 2^16 pairs follow the bit-0 characteristic (128 expected); correlation of bit 21 is 0.1252 (0.125 expected) | `python3 check_spec.py` | `[7]` |
| Theorem 1, Table 2, Fig. 2 | Minimum linear weight `2R` for 8 <= R <= 34; minimum differential weight `3R - log2 6` for 6 <= R <= 34; `R` active S-boxes | `./run_all.sh` | `matsui_eelbc_*.log`, lines `R= .. B= ..` |
| Sect. 3.1 | The search reproduces PRESENT: differential 2, 4, 8, 12, 20, 24, 28; linear 2, 4, 8, 12, 16, 20, 24; 10 active S-boxes over 5 rounds | `./run_all.sh` | `matsui_present_*.log` |
| Sect. 3.1 | Linear hull over masks with at most two active S-boxes: 2^-59.29 for 30 rounds (EE-LBC), 2^-81.64 (PRESENT) | `./run_all.sh` | `hull_eelbc_lin_21.log`, `hull_present_lin_21.log`, line `R=30` |
| Fig. 2 | Plot data | `python3 make_figs.py` | `figs/fig_weights.tex` |
| Sect. 3.2, Table 3 | Measured correlation of the bit-21 trail (within 0.25 bits of 2^-L in RMS) and probability of the bit-0 characteristic (within 0.1 bits of 2^-3L) for 2 to 13 trail rounds | `./run_all.sh --full` | `experiments.log`, `summarize_experiments.log` `[1]` `[2]`, `figs/tab_exp_rows.tex` |
| Sect. 4.1 | Distinguisher with 2^63 known plaintexts: both error probabilities at most 0.10 for half of the keys (about 0.08 for the single-trail correlation 2^-30) | `python3 complexity.py` | `complexity.log`, `[6]` (`[3]`) |
| Sect. 4.2 | Selcuk's model: 2^61.26 known plaintexts for the 20-bit guess with success 0.9 (not used for Table 4) | `python3 complexity.py` | `complexity.log`, `[2]` |
| Sect. 4.3, Table 5 | Rank of the correct 20-bit guess on 12, 14 and 16 rounds: with `N = 2^6/c^2`, among the best 2^10 in 62 of 64 runs and first in 48; with `2^7/c^2` on 12 rounds, first in 30 of 32 | `./run_all.sh --full` | `keyrec_r12_32seeds.log`, `keyrec_longer.log`, `summarize_experiments.log` `[3]`, `figs/tab_keyrec_rows.tex` |
| Sect. 4.3 | The fixed-key hull has the sign of the measured correlation for all 48 keys of `corr_16keys.log`; median deviation 0.18 to 0.31 bits | `./run_all.sh` (measurements: `--full`) | `fixedhull_check.log`, `summarize_experiments.log` `[5]` |
| Sect. 4.3 | Over 2000 keys: median of `|c|` is 2^-28.13 at 28 rounds and 2^-30.17 at 30; 45%, 62%, 75% of the keys have `|c|` >= 2^-28, 2^-28.5, 2^-29 | `./run_all.sh` | `fixedhull_stats.log` |
| Table 4 | For half of the keys: 2^62.26 and 2^63.26 known plaintexts (= 2^6/c^2 and 2^7/c^2 with the median `|c|`) | `python3 complexity.py` | `complexity.log`, `[6]` |
| Sect. 6 | With the PRESENT S-box, a final key, and the PRESENT readings of the open points, the code reproduces the four PRESENT-80 test vectors | `python3 check_spec.py` | `[9]` |

Run times on one PC (20 threads): `check_spec.py` about 15 s, `./run_all.sh` about 2 min,
`./run_all.sh --full` about two hours more (the longest single runs are 2^34 encryptions of 14 rounds and the key
recovery on 16 rounds with 2^32 known plaintexts).

## Files

- `eelbc.py` — reference implementation with the open points as options; self-test.
- `check_spec.py` — facts that need no search: S-box tables, PRESENT conditions, permutation, the example of the proposal, decryption, Proposition 1, PRESENT-80.
- `matsui.c` — exact minimum weights and numbers of active S-boxes by branch-and-bound (EE-LBC and PRESENT).
- `hull.c` — sum of squared correlations (or probabilities) over all trails with at most two active S-boxes per round.
- `experiments.c` — reduced-round measurements of Table 3 (`lin`, `diff`) and the key recovery of Table 5 (`keyrec`: 20 guessed key bits of rounds 1 and 2, `r - 3` trail rounds, keyless last round).
- `fixedhull.c` — correlation for a fixed key as the signed sum over trails with at most two active S-boxes per round (`check`: the keys of the measurements; `stats`: distribution over many keys).
- `complexity.py` — data and time of the full-round attacks: Selcuk's model, the distinguisher, and the complexities for half of the keys.
- `summarize_experiments.py`, `make_figs.py` — Tables 3 and 5 and Fig. 2 from the logs.
- `run_all.sh` — builds and runs everything; each C file also gives its build command in its header.

## License

MIT (see `LICENSE`).
