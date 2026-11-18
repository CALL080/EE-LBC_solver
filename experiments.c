/* Reduced-round experiments on EE-LBC (default reading of code/eelbc.py: 80-bit rotation by 61,
 * Key64 through the S-box layer, counter i, no key after the last round).
 *
 *   gcc -O3 -march=native -fopenmp -o experiments experiments.c
 *   ./experiments check            -> E(0) under key 0 (must equal eelbc.py: 1b7b51a05ea03b5c)
 *   ./experiments lin  <r> <log2N> <seed> [indep]   -> correlation of p[21] xor y[21],  y = S^-1(P^-1(C))
 *   ./experiments diff <r> <log2N> <seed> [indep]   -> fraction of pairs (p, p^1) with y-difference 1
 * r = number of rounds; the last round is stripped without key, so the trail covers r-1 S-layers.
 * 'indep' replaces the key schedule by independent uniform round keys.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static const int SB[16] = {0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2};
static int SI[16];
static uint16_t S16[65536], SI16[65536];
static uint64_t PT[8][256], PIT[8][256];

static uint64_t sub(uint64_t x, const uint16_t *t) {
    return (uint64_t)t[x & 0xFFFF] | (uint64_t)t[(x >> 16) & 0xFFFF] << 16 | (uint64_t)t[(x >> 32) & 0xFFFF] << 32 | (uint64_t)t[x >> 48] << 48;
}
static uint64_t perm(uint64_t x, uint64_t T[8][256]) {
    uint64_t y = 0;
    for (int k = 0; k < 8; k++) y |= T[k][(x >> (8 * k)) & 0xFF];
    return y;
}
static void init(void) {
    for (int x = 0; x < 16; x++) SI[SB[x]] = x;
    for (int x = 0; x < 65536; x++) {
        S16[x] = SI16[x] = 0;
        for (int j = 0; j < 4; j++) { S16[x] |= SB[(x >> (4 * j)) & 15] << (4 * j); SI16[x] |= SI[(x >> (4 * j)) & 15] << (4 * j); }
    }
    for (int k = 0; k < 8; k++)
        for (int v = 0; v < 256; v++) {
            PT[k][v] = PIT[k][v] = 0;
            for (int i = 0; i < 8; i++) if (v >> i & 1) {
                int b = 8 * k + i, pb = 16 * (b % 4) + b / 4, qb = 4 * (b % 16) + b / 16;  /* P and P^-1 */
                PT[k][v] |= 1ULL << pb;
                PIT[k][v] |= 1ULL << qb;
            }
        }
}
static void round_keys(uint64_t khi16, uint64_t klo64, int n, uint64_t *rk) {
    typedef unsigned __int128 u128;
    u128 M = ((u128)1 << 80) - 1, k = ((u128)khi16 << 64 | klo64) & M;
    rk[0] = (uint64_t)(k >> 16);
    for (int i = 2; i <= n; i++) {
        k = ((k << 61) | (k >> 19)) & M;
        uint64_t top = sub((uint64_t)(k >> 16), S16);
        k = ((u128)top << 16) | (k & 0xFFFF);
        k ^= (u128)(i & 0x1F) << 15;
        rk[i - 1] = (uint64_t)(k >> 16);
    }
}
static uint64_t enc(uint64_t x, const uint64_t *rk, int r) {
    for (int i = 0; i < r; i++) x = perm(sub(x ^ rk[i], S16), PT);
    return x;
}
static uint64_t strip(uint64_t c) { return sub(perm(c, PIT), SI16); }   /* = state after r-1 rounds xor K_r */

static uint64_t splitmix(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

int main(int argc, char **argv) {
    init();
    uint64_t rk[32];
    if (argc >= 2 && strcmp(argv[1], "check") == 0) {
        round_keys(0, 0, 31, rk);
        printf("key 0: E(0) = %016llx (eelbc.py: 1b7b51a05ea03b5c)\n", (unsigned long long)enc(0, rk, 31));
        round_keys(0xFFFF, ~0ULL, 31, rk);
        printf("key 1..1: E(1..1) = %016llx (eelbc.py: bc921e297eac883f)\n", (unsigned long long)enc(~0ULL, rk, 31));
        return 0;
    }
    if (argc >= 5 && strcmp(argv[1], "keyrec") == 0) {
        /* Algorithm 2 with 20 guessed round-key bits: K_1 nibbles 4..7 and K_2 nibble 5.
         * r rounds = 2 guessed + (r-3) trail layers + 1 stripped round.  Statistic: x_2[21] xor y[21]. */
        int r = atoi(argv[2]), lg = atoi(argv[3]);
        uint64_t seed = strtoull(argv[4], 0, 0), s = seed;
        uint64_t khi = splitmix(&s) & 0xFFFF, klo = splitmix(&s);
        round_keys(khi, klo, 31, rk);
        static long long cnt[1 << 16];                      /* sum of (-1)^{y[21]} per 16 plaintext bits 16..31 */
        memset(cnt, 0, sizeof cnt);
        long long N = 1LL << lg;
#pragma omp parallel
        {
            static long long loc[64][1 << 16];
            int id = 0;
#ifdef _OPENMP
            id = __builtin_omp_get_thread_num();
#endif
            memset(loc[id], 0, sizeof loc[id]);
#pragma omp for schedule(static)
            for (long long t = 0; t < N; t++) {
                uint64_t st = seed * 0x100000001B3ULL ^ (uint64_t)t * 0x9E3779B97F4A7C15ULL;
                uint64_t p = splitmix(&st), y = strip(enc(p, rk, r));
                loc[id][(p >> 16) & 0xFFFF] += ((y >> 21) & 1) ? -1 : 1;
            }
#pragma omp critical
            for (int u = 0; u < (1 << 16); u++) cnt[u] += loc[id][u];
        }
        int true_k1 = (int)((rk[0] >> 16) & 0xFFFF), true_k2 = (int)((rk[1] >> 20) & 15);
        static double score[1 << 20];
#pragma omp parallel for schedule(dynamic, 16)
        for (int g1 = 0; g1 < (1 << 16); g1++) {
            /* x_1 nibble 5 = bit 1 of S(p_j xor k_j) for j = 4..7, placed at bits 0..3 of the nibble */
            static __thread int nib5[1 << 16];
            for (int u = 0; u < (1 << 16); u++) {
                int v = u ^ g1, z = 0;
                for (int j = 0; j < 4; j++) z |= ((SB[(v >> (4 * j)) & 15] >> 1) & 1) << j;
                nib5[u] = z;
            }
            for (int g2 = 0; g2 < 16; g2++) {
                long long acc = 0;
                for (int u = 0; u < (1 << 16); u++) acc += ((SB[nib5[u] ^ g2] >> 1) & 1) ? -cnt[u] : cnt[u];
                score[(g1 << 4) | g2] = fabs((double)acc) / N;
            }
        }
        int truek = (true_k1 << 4) | true_k2;
        /* rank = 1 + number of guesses with a strictly larger score; ties (equal score) are counted separately,
         * and best_wrong is the largest score of a wrong guess, so that rank 1 with ties = 0 is a strict first place */
        long long rank = 1, ties = 0;
        double best_wrong = 0;
        for (int g = 0; g < (1 << 20); g++) {
            if (g == truek) continue;
            if (score[g] > score[truek]) rank++;
            else if (score[g] == score[truek]) ties++;
            if (score[g] > best_wrong) best_wrong = score[g];
        }
        printf("keyrec r=%2d trail_layers=%2d N=2^%d seed=%llu  |corr(true key)|=2^%.3f (expected 2^%d)  rank of true 20-bit key = %lld / 2^20"
               "  ties=%lld  best wrong 2^%.3f\n",
               r, r - 3, lg, (unsigned long long)seed, log2(score[truek]), -(r - 3), rank, ties, log2(best_wrong));
        return 0;
    }
    if (argc < 5) { fprintf(stderr, "usage: see header\n"); return 1; }
    int lin = strcmp(argv[1], "lin") == 0, r = atoi(argv[2]), lg = atoi(argv[3]);
    uint64_t seed = strtoull(argv[4], 0, 0), s = seed;
    int indep = argc > 5 && strcmp(argv[5], "indep") == 0;
    uint64_t khi = splitmix(&s) & 0xFFFF, klo = splitmix(&s);
    if (indep) for (int i = 0; i < 31; i++) rk[i] = splitmix(&s);
    else round_keys(khi, klo, 31, rk);
    long long N = 1LL << lg, cnt = 0;
#pragma omp parallel for reduction(+ : cnt) schedule(static)
    for (long long t = 0; t < N; t++) {
        uint64_t st = seed * 0x100000001B3ULL ^ (uint64_t)t * 0x9E3779B97F4A7C15ULL;
        uint64_t p = splitmix(&st);
        if (lin) {
            uint64_t y = strip(enc(p, rk, r));
            cnt += (((p ^ y) >> 21) & 1) == 0;
        } else {
            uint64_t y0 = strip(enc(p, rk, r)), y1 = strip(enc(p ^ 1, rk, r));
            cnt += (y0 ^ y1) == 1;
        }
    }
    if (lin) {
        double c = 2.0 * cnt / N - 1.0;
        printf("lin  r=%2d layers=%2d N=2^%d seed=%llu %s  corr=%+.6e  log2|c|=%7.3f  expected log2|c|=%d  (sigma=2^%.1f)\n",
               r, r - 1, lg, (unsigned long long)seed, indep ? "indep" : "ks   ", c, log2(fabs(c)), -(r - 1), -lg / 2.0);
    } else {
        printf("diff r=%2d layers=%2d N=2^%d seed=%llu %s  hits=%lld  log2(p)=%7.3f  expected 1-bit trail log2(p)=%d\n",
               r, r - 1, lg, (unsigned long long)seed, indep ? "indep" : "ks   ", cnt, cnt ? log2((double)cnt / N) : -INFINITY, -3 * (r - 1));
    }
    return 0;
}
