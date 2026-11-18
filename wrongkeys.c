/* Correlation between the bit x_{2,21} predicted by a wrong 20-bit guess and by the correct one (Sect. 4.2).
 * x_{2,21} = bit 1 of S(z xor g2), z_j = bit 1 of S(u_j xor g1_j) for the four plaintext nibbles u_4..u_7.
 * Over uniform plaintext bits u, the correlation of two guesses factorizes over the nibbles, so all 2^20 guesses
 * are evaluated exactly.  rho does not depend on the trail rounds; the statistic of a wrong guess has mean rho * c.
 *
 *   gcc -O3 -o wrongkeys wrongkeys.c -lm
 *   ./wrongkeys <ntrue> <seed>      distribution of |rho| over all wrong guesses, for ntrue random correct guesses
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

static const int SB[16] = {0xC, 0x5, 0x6, 0xB, 0x7, 0x0, 0xA, 0xD, 0x1, 0xE, 0xF, 0x8, 0x4, 0x9, 0x3, 0x2};
static int b1(int x) { return (SB[x] >> 1) & 1; }

int main(int argc, char **argv) {
    int ntrue = argc > 1 ? atoi(argv[1]) : 4;
    uint64_t s = argc > 2 ? strtoull(argv[2], 0, 0) : 1;
    /* P[g][h][a][b] = #{u : b1(u^g) = a, b1(u^h) = b} / 16 */
    static double P[16][16][2][2];
    for (int g = 0; g < 16; g++) for (int h = 0; h < 16; h++) for (int u = 0; u < 16; u++)
        P[g][h][b1(u ^ g)][b1(u ^ h)] += 1.0 / 16;
    for (int t = 0; t < ntrue; t++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        int true1 = (int)(s >> 40) & 0xFFFF, true2 = (int)(s >> 20) & 15;
        int hist[129] = {0};            /* |rho| in steps of 1/128 */
        int ties = 0;
        double maxr = 0;
        for (int g = 0; g < (1 << 20); g++) {
            int g1 = g >> 4, g2 = g & 15;
            if (g1 == true1 && g2 == true2) continue;
            double rho = 0;
            for (int z = 0; z < 16; z++) for (int zt = 0; zt < 16; zt++) {
                double p = 1;
                for (int j = 0; j < 4; j++) p *= P[(g1 >> (4 * j)) & 15][(true1 >> (4 * j)) & 15][(z >> j) & 1][(zt >> j) & 1];
                if (p == 0) continue;
                rho += ((b1(z ^ g2) ^ b1(zt ^ true2)) ? -p : p);
            }
            double a = fabs(rho);
            if (a > 1 - 1e-12) ties++;
            if (a > maxr) maxr = a;
            hist[(int)lround(a * 128)]++;
        }
        printf("correct guess K1[4..7]=%04x K2[5]=%x: max |rho| over wrong guesses = %.4f; |rho| = 1 for %d guesses;", true1, true2, maxr, ties);
        printf(" counts with |rho| >= 0.75: %d, >= 0.5: %d\n",
               hist[96] + hist[97] + hist[98] + hist[99] + hist[100] + hist[101] + hist[102] + hist[103] + hist[104] + hist[105] +
               hist[106] + hist[107] + hist[108] + hist[109] + hist[110] + hist[111] + hist[112] + hist[113] + hist[114] + hist[115] +
               hist[116] + hist[117] + hist[118] + hist[119] + hist[120] + hist[121] + hist[122] + hist[123] + hist[124] + hist[125] +
               hist[126] + hist[127] + hist[128],
               ({ int c = 0; for (int k = 64; k <= 128; k++) c += hist[k]; c; }));
    }
    return 0;
}
