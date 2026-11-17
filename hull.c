/* Linear hull (ELP) and differential (EDP) of EE-LBC / PRESENT restricted to masks (differences) with at most
 * two active S-boxes at every S-layer input.  Every term is a nonnegative trail contribution, so the result is a
 * lower bound on ELP = sum over all trails of c^2 (Nyberg's linear-hull theorem, independent round keys) and on EDP.
 *
 *   gcc -O3 -o hull hull.c -lm
 *   ./hull <eelbc|present> <lin|diff> <in_bit> <out_bit> <Rmax>
 * Prints log2 of the restricted ELP/EDP from e_in to e_out over R S-layers, and the single-trail value for comparison.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static const int S_EE[16] = {0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2};
static const int S_PR[16] = {0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2};
static double T[16][16];   /* c^2 or probability */
#define NST (240 + 120 * 225)
static int pairidx[16][16];

static int idx_of(uint64_t m) {   /* -1 if more than two active S-boxes */
    int act[3], n = 0;
    for (int j = 0; j < 16; j++) if ((m >> (4 * j)) & 15) { if (n == 2) return -1; act[n++] = j; }
    if (n == 1) return act[0] * 15 + (int)((m >> (4 * act[0])) & 15) - 1;
    int v1 = (m >> (4 * act[0])) & 15, v2 = (m >> (4 * act[1])) & 15;
    return 240 + pairidx[act[0]][act[1]] * 225 + (v1 - 1) * 15 + (v2 - 1);
}
static uint64_t mask_of(int s) {
    if (s < 240) return (uint64_t)(s % 15 + 1) << (4 * (s / 15));
    s -= 240;
    int p = s / 225, r = s % 225, j1 = 0, j2 = 0;
    for (int a = 0; a < 16; a++) for (int b = a + 1; b < 16; b++) if (pairidx[a][b] == p) { j1 = a; j2 = b; }
    return (uint64_t)(r / 15 + 1) << (4 * j1) | (uint64_t)(r % 15 + 1) << (4 * j2);
}
static uint64_t P(uint64_t x) {
    uint64_t y = 0;
    for (int i = 0; i < 64; i++) if (x >> i & 1) y |= 1ULL << (16 * (i % 4) + i / 4);
    return y;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s <eelbc|present> <lin|diff> <in_bit> <out_bit> <Rmax>\n", argv[0]); return 1; }
    const int *S = strcmp(argv[1], "present") == 0 ? S_PR : S_EE;
    int lin = strcmp(argv[2], "lin") == 0, ib = atoi(argv[3]), ob = atoi(argv[4]), Rmax = atoi(argv[5]);
    for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++) {
            int v = 0;
            if (lin) { for (int x = 0; x < 16; x++) v += (__builtin_popcount(a & x) + __builtin_popcount(b & S[x])) % 2 == 0; v -= 8; T[a][b] = (v / 8.0) * (v / 8.0); }
            else { for (int x = 0; x < 16; x++) v += (S[x] ^ S[x ^ a]) == b; T[a][b] = v / 16.0; }
        }
    int k = 0;
    for (int a = 0; a < 16; a++) for (int b = a + 1; b < 16; b++) pairidx[a][b] = k++;
    /* sparse transition lists: state -> (next state, weight) */
    int *cnt = calloc(NST, sizeof(int)), **to = malloc(NST * sizeof(int *));
    double **w = malloc(NST * sizeof(double *));
    for (int s = 0; s < NST; s++) {
        uint64_t m = mask_of(s);
        int act[2], n = 0;
        for (int j = 0; j < 16; j++) if ((m >> (4 * j)) & 15) act[n++] = j;
        to[s] = malloc(225 * sizeof(int));
        w[s] = malloc(225 * sizeof(double));
        int a0 = (m >> (4 * act[0])) & 15, a1 = n > 1 ? (m >> (4 * act[1])) & 15 : 0;
        for (int b0 = 1; b0 < 16; b0++) {
            if (T[a0][b0] == 0) continue;
            for (int b1 = (n > 1 ? 1 : 0); b1 < (n > 1 ? 16 : 1); b1++) {
                double x = T[a0][b0] * (n > 1 ? T[a1][b1] : 1.0);
                if (x == 0) continue;
                uint64_t out = (uint64_t)b0 << (4 * act[0]) | (n > 1 ? (uint64_t)b1 << (4 * act[1]) : 0);
                int t = idx_of(P(out));
                if (t < 0) continue;
                to[s][cnt[s]] = t; w[s][cnt[s]] = x; cnt[s]++;
            }
        }
    }
    double *v = calloc(NST, sizeof(double)), *nv = calloc(NST, sizeof(double));
    int sin = idx_of(1ULL << ib), sout = idx_of(1ULL << ob);
    double single = lin ? T[2][2] : T[1][1];
    v[sin] = 1.0;
    printf("# %s %s, e_%d -> e_%d, masks with <= 2 active S-boxes at every S-layer input\n", argv[1], argv[2], ib, ob);
    for (int R = 1; R <= Rmax; R++) {
        memset(nv, 0, NST * sizeof(double));
        for (int s = 0; s < NST; s++) if (v[s] != 0) for (int i = 0; i < cnt[s]; i++) nv[to[s][i]] += v[s] * w[s][i];
        double *tmp = v; v = nv; nv = tmp;
        /* the last S-layer's output must be e_out after P: v holds masks at the next S-layer input */
        printf("R=%2d  log2(restricted sum) = %8.3f   log2(single 1-bit trail) = %8.3f\n", R, log2(v[sout]), R * log2(single));
    }
    return 0;
}
