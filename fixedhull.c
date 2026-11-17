/* Fixed-key linear hull of EE-LBC from bit 21 to bit 21, restricted to masks with at most two active S-boxes
 * at every S-layer input, with the signs given by the actual round keys:
 *   c(K) ~ sum over trails  prod_r (-1)^{<u_r, K_r>} C_S(u_r -> v_r),   u_{r+1} = P(v_r).
 * The keys are those of experiments.c (same splitmix seeds and key schedule), so the values can be compared
 * with the measured correlations, which also carry the sign (-1)^{K_{L+1}[21]} of the last key addition.
 *
 *   gcc -O3 -march=native -fopenmp -o fixedhull fixedhull.c -lm
 *   ./fixedhull check <L> <seed_from> <seed_to>     signed prediction of the measured correlation for these seeds
 *   ./fixedhull stats <L> <nkeys>                    quantiles of log2|c| + L over nkeys keys (seeds 1000000..)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static const int SB[16] = {0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2};
#define NST (240 + 120 * 225)
static int pairidx[16][16];
static double CS[16][16];
static int *cnt, **to;
static double **w;
static uint64_t smask[NST];

static int idx_of(uint64_t m) {
    int act[3], n = 0;
    for (int j = 0; j < 16; j++) if ((m >> (4 * j)) & 15) { if (n == 2) return -1; act[n++] = j; }
    if (n == 1) return act[0] * 15 + (int)((m >> (4 * act[0])) & 15) - 1;
    int v1 = (m >> (4 * act[0])) & 15, v2 = (m >> (4 * act[1])) & 15;
    return 240 + pairidx[act[0]][act[1]] * 225 + (v1 - 1) * 15 + (v2 - 1);
}
static uint64_t P(uint64_t x) {
    uint64_t y = 0;
    for (int i = 0; i < 64; i++) if (x >> i & 1) y |= 1ULL << (16 * (i % 4) + i / 4);
    return y;
}

/* key schedule and seeds exactly as in experiments.c */
static uint16_t S16[65536];
static uint64_t sub(uint64_t x) {
    return (uint64_t)S16[x & 0xFFFF] | (uint64_t)S16[(x >> 16) & 0xFFFF] << 16 | (uint64_t)S16[(x >> 32) & 0xFFFF] << 32 | (uint64_t)S16[x >> 48] << 48;
}
static void round_keys(uint64_t khi16, uint64_t klo64, int n, uint64_t *rk) {
    typedef unsigned __int128 u128;
    u128 M = ((u128)1 << 80) - 1, k = ((u128)khi16 << 64 | klo64) & M;
    rk[0] = (uint64_t)(k >> 16);
    for (int i = 2; i <= n; i++) {
        k = ((k << 61) | (k >> 19)) & M;
        uint64_t top = sub((uint64_t)(k >> 16));
        k = ((u128)top << 16) | (k & 0xFFFF);
        k ^= (u128)(i & 0x1F) << 15;
        rk[i - 1] = (uint64_t)(k >> 16);
    }
}
static uint64_t splitmix(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static double hull(int L, const uint64_t *rk, double *v, double *nv) {
    memset(v, 0, NST * sizeof(double));
    v[idx_of(1ULL << 21)] = 1.0;
    for (int r = 0; r < L; r++) {
        memset(nv, 0, NST * sizeof(double));
        for (int s = 0; s < NST; s++) {
            if (v[s] == 0) continue;
            double x = (__builtin_popcountll(smask[s] & rk[r]) & 1) ? -v[s] : v[s];
            for (int i = 0; i < cnt[s]; i++) nv[to[s][i]] += x * w[s][i];
        }
        double *t = v; v = nv; nv = t;
    }
    return v[idx_of(1ULL << 21)];
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: see header\n"); return 1; }
    for (int x = 0; x < 65536; x++) { S16[x] = 0; for (int j = 0; j < 4; j++) S16[x] |= SB[(x >> (4 * j)) & 15] << (4 * j); }
    for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++) {
            int agree = 0;
            for (int x = 0; x < 16; x++) agree += (__builtin_popcount(a & x) + __builtin_popcount(b & SB[x])) % 2 == 0;
            CS[a][b] = (2.0 * agree - 16.0) / 16.0;
        }
    int k = 0;
    for (int a = 0; a < 16; a++) for (int b = a + 1; b < 16; b++) pairidx[a][b] = k++;
    cnt = calloc(NST, sizeof(int)); to = malloc(NST * sizeof(int *)); w = malloc(NST * sizeof(double *));
    for (int s = 0; s < NST; s++) {
        uint64_t m = 0;
        if (s < 240) m = (uint64_t)(s % 15 + 1) << (4 * (s / 15));
        else {
            int q = s - 240, p = q / 225, rr = q % 225, j1 = 0, j2 = 0;
            for (int a = 0; a < 16; a++) for (int b = a + 1; b < 16; b++) if (pairidx[a][b] == p) { j1 = a; j2 = b; }
            m = (uint64_t)(rr / 15 + 1) << (4 * j1) | (uint64_t)(rr % 15 + 1) << (4 * j2);
        }
        smask[s] = m;
        int act[2], n = 0;
        for (int j = 0; j < 16; j++) if ((m >> (4 * j)) & 15) act[n++] = j;
        to[s] = malloc(225 * sizeof(int)); w[s] = malloc(225 * sizeof(double));
        int a0 = (m >> (4 * act[0])) & 15, a1 = n > 1 ? (m >> (4 * act[1])) & 15 : 0;
        for (int b0 = 1; b0 < 16; b0++) {
            if (CS[a0][b0] == 0) continue;
            for (int b1 = (n > 1 ? 1 : 0); b1 < (n > 1 ? 16 : 1); b1++) {
                double x = CS[a0][b0] * (n > 1 ? CS[a1][b1] : 1.0);
                if (x == 0) continue;
                uint64_t out = (uint64_t)b0 << (4 * act[0]) | (n > 1 ? (uint64_t)b1 << (4 * act[1]) : 0);
                int t = idx_of(P(out));
                if (t < 0) continue;
                to[s][cnt[s]] = t; w[s][cnt[s]] = x; cnt[s]++;
            }
        }
    }
    int L = atoi(argv[2]);
    if (strcmp(argv[1], "check") == 0) {
        int from = atoi(argv[3]), upto = argc > 4 ? atoi(argv[4]) : from;
        double *v = malloc(NST * sizeof(double)), *nv = malloc(NST * sizeof(double));
        for (int seed = from; seed <= upto; seed++) {
            uint64_t s = (uint64_t)seed, rk[32];
            uint64_t khi = splitmix(&s) & 0xFFFF, klo = splitmix(&s);
            round_keys(khi, klo, 31, rk);
            double c = hull(L, rk, v, nv);
            double meas = ((rk[L] >> 21) & 1) ? -c : c;   /* sign of the key bit added before the stripped round */
            printf("check L=%2d seed=%d  predicted corr=%+.6e  log2|c|+L=%+.3f\n", L, seed, meas, log2(fabs(c)) + L);
        }
    } else {
        int nk = atoi(argv[3]);
        double *res = malloc(nk * sizeof(double));
#pragma omp parallel
        {
            double *v = malloc(NST * sizeof(double)), *nv = malloc(NST * sizeof(double));
#pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < nk; i++) {
                uint64_t s = 1000000ULL + i, rk[32];
                uint64_t khi = splitmix(&s) & 0xFFFF, klo = splitmix(&s);
                round_keys(khi, klo, 31, rk);
                res[i] = log2(fabs(hull(L, rk, v, nv))) + L;
            }
        }
        int cmp(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
        qsort(res, nk, sizeof(double), cmp);
        int ge0 = 0, gem05 = 0, gem1 = 0;
        for (int i = 0; i < nk; i++) { ge0 += res[i] >= 0; gem05 += res[i] >= -0.5; gem1 += res[i] >= -1; }
        printf("stats L=%2d keys=%d  log2|c|+L: min %+.2f  Q1 %+.2f  median %+.2f  Q3 %+.2f  max %+.2f;"
               "  fraction >= 0: %.3f, >= -0.5: %.3f, >= -1: %.3f\n", L, nk, res[0], res[nk / 4], res[nk / 2], res[3 * nk / 4],
               res[nk - 1], (double)ge0 / nk, (double)gem05 / nk, (double)gem1 / nk);
    }
    return 0;
}
