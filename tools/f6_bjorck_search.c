/*
 * f6_bjorck_search.c
 *
 * Exhaustive search for the 48 Bjorck biunimodular sequences of length 6.
 *
 * A sequence z is biunimodular when |z_a| = 1 and |(F z)_a| = 1, where F
 * is the unnormalised 6-point Fourier matrix.  The known algebraic
 * alphabet is
 *
 *     A = { d^m * zeta^k : -2 <= m <= 2, 0 <= k < 12 },
 *     d = (1 - sqrt(3))/2 + i*sqrt(sqrt(3)/2),
 *     zeta = exp(2*pi*i/12).
 *
 * We fix z_0 = 1 (global phase) and enumerate the remaining five
 * coordinates.  The search returns exactly 48 sequences.  The codes
 * printed here are the ones embedded in python/generate_facts.py and
 * lean4/MUB6/Certificates.lean.
 *
 * Build and run:
 *   cc -O3 -o f6_search f6_bjorck_search.c -lm
 *   ./f6_search > f6_solutions.txt
 */

#include <stdio.h>
#include <complex.h>
#include <math.h>

#define M 5
#define K 12
#define N (M * K)   /* 60 */

int main(void) {
    const double pi = acos(-1.0);
    const double complex zeta = cexp(I * 2 * pi / 12.0);
    const double complex d =
        (1.0 - sqrt(3.0)) / 2.0 + I * sqrt(sqrt(3.0) / 2.0);

    double complex A[N];
    int t = 0;
    for (int m = -2; m <= 2; ++m) {
        double complex dm = cpow(d, (double)m);
        for (int k = 0; k < 12; ++k) {
            A[t++] = dm * cpow(zeta, k);
        }
    }

    const double complex q = cexp(I * 2 * pi / 6.0);
    double complex C[6][6][N];
    for (int a = 0; a < 6; ++a) {
        for (int b = 0; b < 6; ++b) {
            double complex qq = cpow(q, (double)(a * b));
            for (int j = 0; j < N; ++j) {
                C[a][b][j] = qq * A[j];
            }
        }
    }

    /* m = 0, k = 0 is the scalar 1: index (2 * 12) = 24. */
    const int one_idx = 24;
    int idx[6];
    idx[0] = one_idx;
    long long solutions = 0;

    for (idx[1] = 0; idx[1] < N; ++idx[1])
    for (idx[2] = 0; idx[2] < N; ++idx[2])
    for (idx[3] = 0; idx[3] < N; ++idx[3])
    for (idx[4] = 0; idx[4] < N; ++idx[4])
    for (idx[5] = 0; idx[5] < N; ++idx[5]) {
        double complex s[6];
        for (int a = 0; a < 6; ++a) {
            s[a] = A[one_idx]
                 + C[a][1][idx[1]] + C[a][2][idx[2]]
                 + C[a][3][idx[3]] + C[a][4][idx[4]]
                 + C[a][5][idx[5]];
        }
        int ok = 1;
        for (int a = 0; a < 6; ++a) {
            double m = creal(s[a]) * creal(s[a])
                     + cimag(s[a]) * cimag(s[a]);
            if (fabs(m - 6.0) > 1e-8) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            ++solutions;
            printf("%d %d %d %d %d %d\n",
                   idx[0], idx[1], idx[2], idx[3], idx[4], idx[5]);
        }
    }

    fprintf(stderr, "solutions: %lld\n", solutions);
    return 0;
}
