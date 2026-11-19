/* Exact minimum weights of differential and linear trails for 64-bit SPNs with the PRESENT bit permutation
 * (EE-LBC and PRESENT), by Matsui's branch-and-bound algorithm (Matsui, EUROCRYPT 1994, Sect. 4).
 *
 *   gcc -O3 -march=native -o matsui matsui.c -lm
 *   ./matsui <eelbc|present> <lin|diff|aslin|asdiff> <Rmax> [time_limit_seconds_per_round]
 *
 * Weight of a trail = sum over active S-boxes of -log2(DDT/16) (diff) or -log2((LAT/8)^2) (lin, LAT = #agree-8),
 * i.e. the trail probability / squared correlation of the usual Markov estimate with independent round keys.
 * Modes aslin/asdiff count active S-boxes instead (weight 1 per active S-box).
 * An R-round trail has R S-layers: a_1 -S-> b_1 -P-> a_2 -S-> ... -S-> b_R.
 * B[R] is computed from B[1..R-1] (exact), by iterative deepening on the bound; the printed trail attains B[R].
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

static const int S_EE[16] = {0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2};
static const int S_PR[16] = {0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2};

#define EPS 1e-9
#define MAXR 64
static double W[16][16];                 /* weight of a->b, INFINITY if impossible */
static int NOUT[16], OUTB[16][16];       /* outputs of input a sorted by weight */
static double OUTW[16][16];
static int ORD_OUT[15];                  /* nonzero outputs sorted by WMIN_OUT (first round) */
static double WMIN_OUT[16], WMIN_IN[16], WMIN;
static uint64_t PT[16][16];              /* P applied to value v in S-box j */
static double B[MAXR + 1];
static double Bound;
static int R;
static uint64_t cur_a[MAXR + 2], cur_b[MAXR + 2], best_a[MAXR + 2], best_b[MAXR + 2];
static int found;
static long long nodes;
static clock_t t0;
static double tlimit;
static int timed_out;

static inline int nib(uint64_t x, int j) { return (int)((x >> (4 * j)) & 0xF); }

static void record(double w) {
    Bound = w;
    found = 1;
    memcpy(best_a, cur_a, sizeof cur_a);
    memcpy(best_b, cur_b, sizeof cur_b);
}

static void start_round(int r, uint64_t a, double acc);

/* enumerate outputs of the k-th active S-box of round r */
static void enum_sbox(int r, uint64_t a, const int *act, int nact, int k, uint64_t outb, uint64_t outP, double acc) {
    if (timed_out) return;
    if (k == nact) {
        cur_a[r] = a;
        cur_b[r] = outb;
        start_round(r + 1, outP, acc);
        return;
    }
    if ((++nodes & 0xFFFFF) == 0 && tlimit > 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > tlimit) { timed_out = 1; return; }
    int j = act[k], v = nib(a, j);
    double rest = (nact - k - 1) * WMIN + B[R - r];
    for (int t = 0; t < NOUT[v]; t++) {
        double w = OUTW[v][t];
        if (acc + w + rest >= Bound - EPS) break;
        int b = OUTB[v][t];
        enum_sbox(r, a, act, nact, k + 1, outb | ((uint64_t)b << (4 * j)), outP | PT[j][b], acc + w);
    }
}

static void start_round(int r, uint64_t a, double acc) {
    int act[16], nact = 0;
    for (int j = 0; j < 16; j++) if (nib(a, j)) act[nact++] = j;
    if (r == R) {               /* last round: best output for each active S-box */
        double w = acc;
        uint64_t outb = 0;
        for (int k = 0; k < nact; k++) {
            int v = nib(a, act[k]);
            w += WMIN_IN[v];
            outb |= (uint64_t)OUTB[v][0] << (4 * act[k]);
        }
        if (w < Bound - EPS) { cur_a[r] = a; cur_b[r] = outb; record(w); }
        return;
    }
    if (acc + nact * WMIN + B[R - r] >= Bound - EPS) return;
    enum_sbox(r, a, act, nact, 0, 0, 0, acc);
}

/* first round: choose output values b_1 S-box by S-box, each costing its cheapest input */
static void first_round(int j, uint64_t outb, uint64_t outP, double acc, int nact) {
    if (timed_out) return;
    if (j == 16) {
        if (!nact) return;
        uint64_t a = 0;
        for (int i = 0; i < 16; i++) {
            int b = nib(outb, i);
            if (!b) continue;
            int best = -1;
            for (int x = 1; x < 16; x++) if (W[x][b] <= WMIN_OUT[b] + EPS) { best = x; break; }
            a |= (uint64_t)best << (4 * i);
        }
        cur_a[1] = a;
        cur_b[1] = outb;
        if (R == 1) { if (acc < Bound - EPS) record(acc); return; }
        start_round(2, outP, acc);
        return;
    }
    first_round(j + 1, outb, outP, acc, nact);            /* S-box j inactive */
    for (int t = 0; t < 15; t++) {
        int b = ORD_OUT[t];
        double w = WMIN_OUT[b];
        if (acc + w + B[R - 1] >= Bound - EPS) break;
        first_round(j + 1, outb | ((uint64_t)b << (4 * j)), outP | PT[j][b], acc + w, nact + 1);
    }
}

static void print_trail(void) {
    for (int r = 1; r <= R; r++) {
        printf("    r%-2d", r);
        for (int j = 15; j >= 0; j--) {
            int a = nib(best_a[r], j), b = nib(best_b[r], j);
            if (a) printf(" S%d:%X->%X(%.3f)", j, a, b, W[a][b]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s <eelbc|present> <lin|diff|aslin|asdiff> <Rmax> [tlimit]\n", argv[0]); return 1; }
    const int *S = strcmp(argv[1], "present") == 0 ? S_PR : S_EE;
    const char *mode = argv[2];
    int Rmax = atoi(argv[3]);
    tlimit = argc > 4 ? atof(argv[4]) : 0;
    int lin = strstr(mode, "lin") != NULL, as = strncmp(mode, "as", 2) == 0;

    for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++) {
            int v = 0;
            if (lin) { for (int x = 0; x < 16; x++) v += (__builtin_popcount(a & x) + __builtin_popcount(b & S[x])) % 2 == 0; v -= 8; v = abs(v); }
            else { for (int x = 0; x < 16; x++) v += (S[x] ^ S[x ^ a]) == b; }
            if (a == 0 || b == 0 || v == 0) W[a][b] = INFINITY;
            else if (as) W[a][b] = 1.0;
            else W[a][b] = lin ? -2.0 * log2(v / 8.0) : -log2(v / 16.0);
        }
    WMIN = INFINITY;
    for (int a = 1; a < 16; a++) {
        NOUT[a] = 0;
        for (int b = 1; b < 16; b++) if (isfinite(W[a][b])) { OUTB[a][NOUT[a]] = b; OUTW[a][NOUT[a]] = W[a][b]; NOUT[a]++; }
        for (int i = 1; i < NOUT[a]; i++)            /* insertion sort by weight */
            for (int k = i; k > 0 && OUTW[a][k] < OUTW[a][k - 1]; k--) {
                double tw = OUTW[a][k]; OUTW[a][k] = OUTW[a][k - 1]; OUTW[a][k - 1] = tw;
                int tb = OUTB[a][k]; OUTB[a][k] = OUTB[a][k - 1]; OUTB[a][k - 1] = tb;
            }
        WMIN_IN[a] = OUTW[a][0];
        if (WMIN_IN[a] < WMIN) WMIN = WMIN_IN[a];
    }
    for (int b = 1; b < 16; b++) { WMIN_OUT[b] = INFINITY; for (int a = 1; a < 16; a++) if (W[a][b] < WMIN_OUT[b]) WMIN_OUT[b] = W[a][b]; }
    for (int t = 0; t < 15; t++) ORD_OUT[t] = t + 1;
    for (int i = 1; i < 15; i++)
        for (int k = i; k > 0 && WMIN_OUT[ORD_OUT[k]] < WMIN_OUT[ORD_OUT[k - 1]]; k--) { int tb = ORD_OUT[k]; ORD_OUT[k] = ORD_OUT[k - 1]; ORD_OUT[k - 1] = tb; }
    for (int j = 0; j < 16; j++)
        for (int v = 0; v < 16; v++) {
            uint64_t y = 0;
            for (int i = 0; i < 4; i++) if (v >> i & 1) { int bit = 4 * j + i; y |= 1ULL << (16 * (bit % 4) + bit / 4); }
            PT[j][v] = y;
        }

    printf("# cipher=%s mode=%s  (weight = -log2 of %s)\n", argv[1], mode, as ? "1 per active S-box" : lin ? "squared correlation" : "probability");
    B[0] = 0;
    double step = as ? 1.0 : (lin ? 2.0 : 3.0);
    for (R = 1; R <= Rmax; R++) {
        t0 = clock();
        nodes = 0;
        timed_out = 0;
        double lo = B[R - 1] + WMIN;               /* B[R] >= B[R-1] + B[1] */
        double hi = B[R - 1] + step;
        found = 0;
        while (!found && !timed_out) {
            Bound = hi + EPS;
            first_round(0, 0, 0, 0.0, 0);
            if (!found) hi += 1.0;
        }
        if (timed_out) { printf("R=%2d  TIMEOUT after %.0fs; B[R] >= %.4f\n", R, tlimit, lo); fflush(stdout); break; }
        B[R] = Bound;
        printf("R=%2d  B=%9.4f  (%.1fs, %lld nodes)\n", R, B[R], (double)(clock() - t0) / CLOCKS_PER_SEC, nodes);
        print_trail();
        fflush(stdout);
    }
    return 0;
}
