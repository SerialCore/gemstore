/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2026, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/eigen.h>
#include <gemstore/thread.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

void eigen_tridiagonal(double **a, int n, double *d, double *e, double *et, int lt)
{
    int i, j, k;
    double sigma, beta, uu, eps, vmax;

    double r, c, s, t, ks;
    double d0, d1, e0, e1;
    double vc, vs;
    double *dd, *ee, *ds, *es, *as;
    int *od;
    double temp, g;

    eps = 1E-17;
    e[0] = 0;
    for (i = n - 1; i >= 1; i--) {
        sigma = 0;
        for (j = 0; j <= i - 1; j++) {
            sigma += a[i][j] * a[i][j];
        }
        sigma = sqrt(sigma);
        if (a[i][i - 1] < 0) {
            sigma *= -1;
        }
        d[i] = a[i][i];
        e[i] = -sigma;
        beta = sigma * (sigma + a[i][i - 1]);
        a[i][i - 1] += sigma;

        if (fabs(beta) > eps * fabs(a[i][i])) {
            for (j = 0; j < i; j++) {
                a[j][i] = 0;
                for (k = 0; k < i; k++)
                {
                    a[j][i] += a[j][k] * a[i][k];
                }
            }

            uu = 0;
            for (j = 0; j < i; j++) {
                uu += a[i][j] * a[j][i];
            }
            uu = uu / beta / beta;

            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    a[j][k] = a[j][k] - a[i][j] * a[k][i] / beta - a[i][k] * a[j][i] / beta + a[i][j] * a[i][k] * uu;
                }
            }
            a[i][i] = beta;
        }
        else {
            a[i][i] = 0;
        }
    }
    d[0] = a[0][0];

    if (lt > 0) {
        a[0][0] = 1;
        for (i = 1; i < n; i++) {
            if (0 != a[i][i]) {
                for (j = 0; j < i; j++) {
                    a[j][i] = 0;
                    for (k = 0; k < i; k++) {
                        a[j][i] += a[j][k] * a[i][k];
                    }
                    a[j][i] /= a[i][i];
                }
                for (j = 0; j < i; j++) {
                    for (k = 0; k < i; k++) {
                        a[j][k] -= a[j][i] * a[i][k];
                    }
                }
                for (j = 0; j < i; j++) {
                    a[i][j] = 0;
                    a[j][i] = 0;
                }
            }
            a[i][i] = 1;
        }
    }

    dd = (double *)malloc(sizeof(double) * n);
    ee = (double *)malloc(sizeof(double) * n);
    ds = (double *)malloc(sizeof(double) * n);
    es = (double *)malloc(sizeof(double) * n);
    as = (double *)malloc(sizeof(double) * n);
    od = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
        ds[i] = d[i];
        es[i] = e[i];
        od[i] = i;
    }

    t = 0;
    e[0] = 0;
    for (i = 0; i < n - 1; i++) {
        while (fabs(e[i + 1]) > eps * (fabs(d[i]) + fabs(d[i + 1]))) {
            g = (d[i + 1] - d[i]) / (2 * e[i + 1]);
            if (g >= 0) {
                ks = d[i] - e[i + 1] / (g + sqrt(1 + g * g));
            }
            else {
                ks = d[i] - e[i + 1] / (g - sqrt(1 + g * g));
            }
            for (j = n - 1; j > i; j--) {
                r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                if (0 == r) {
                    continue;
                }
                c = (d[j] - ks) / r;
                s = -e[j] / r;

                d0 = d[j - 1];
                d1 = d[j];

                e0 = e[j - 1];
                e1 = e[j];

                d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                e[j - 1] = c * e0;
                e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                if (j < n - 1) {
                    e[j + 1] = -s * t + c * e[j + 1];
                }

                t = -s * e0;

                if (lt >= n - 1) {
                    for (k = 0; k < n; k++) {
                        vc = a[j - 1][k];
                        vs = a[j][k];
                        a[j - 1][k] = c * vc + s * vs;
                        a[j][k] = -s * vc + c * vs;
                    }
                }
            }
        }
    }

    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (d[j] > d[j + 1]) {
                temp = d[j];
                d[j] = d[j + 1];
                d[j + 1] = temp;

                k = od[j];
                od[j] = od[j + 1];
                od[j + 1] = k;
            }
        }
    }

    for (i = 0; i < n; i++) {
        dd[i] = d[i];
        ee[i] = e[i];
    }

    if (lt >= n - 1) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                a[i][j] = as[od[i]];
            }
        }
    }
    else {
        for (i = 0; i < n; i++) {
            od[i] = -1;
            d[i] = ds[i];
            e[i] = es[i];
        }

        t = 0;
        e[0] = 0;
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                d[j] -= dd[i];
            }
            while (1) {
                int bk = 0;
                for (j = 0; j < n - 1; j++) {
                    if (fabs(e[j + 1]) <= eps * (fabs(d[j]) + fabs(d[j + 1])) && fabs(d[j]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[j] < 0) {
                        od[j] = i;
                        bk = 1;
                        break;
                    }
                }
                if (fabs(d[n - 1]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[n - 1] < 0 && 0 == bk) {
                    od[n - 1] = i;
                    bk = 1;
                }
                if (1 == bk) {
                    break;
                }

                ks = 0;
                for (j = n - 1; j > 0; j--) {
                    if (od[j - 1] < 0) {
                        r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                        if (0 == r) {
                            continue;
                        }
                        c = (d[j] - ks) / r;
                        s = -e[j] / r;

                        d0 = d[j - 1];
                        d1 = d[j];

                        e0 = e[j - 1];
                        e1 = e[j];

                        d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                        d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                        e[j - 1] = c * e0;
                        e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                        if (j < n - 1) {
                            e[j + 1] = -s * t + c * e[j + 1];
                        }

                        t = -s * e0;

                        for (k = 0; k < n; k++) {
                            vc = a[j - 1][k];
                            vs = a[j][k];
                            a[j - 1][k] = c * vc + s * vs;
                            a[j][k] = -s * vc + c * vs;
                        }
                    }
                }
            }
            for (j = 0; j < n; j++) {
                d[j] += dd[i];
            }
        }

        for (i = 0; i < n; i++) {
            as[i] = d[i];
        }

        for (i = 0; i < n; i++) {
            if (od[i] >= 0) {
                d[od[i]] = as[i];
            }
        }
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                if (od[i] >= 0) {
                    a[od[i]][j] = as[i];
                }
            }
        }
    }

    for (i = lt; i < n; i++){
        for (j = 0; j < n; j++) {
            a[i][j] = 0;
        }
    }
    for (i = 0; i < lt; i++) {
        et[i] = e[i];
    }
    for (i = 0; i < n; i++) {
        d[i] = dd[i];
        e[i] = ee[i];
    }
    for (i = 0; i < lt; i++) {
        vmax = 0;
        for (j = 0; j < n; j++) {
            if (fabs(a[i][j]) > fabs(vmax)) {
                vmax = a[i][j];
            }
        }
        if (vmax < 0) {
            for (j = 0; j < n; j++) {
                a[i][j] *= -1;
            }
        }
    }

    free(dd);
    free(ee);
    free(ds);
    free(es);
    free(as);
    free(od);
}

void eigen_standard(double **a, int n, double *d, double **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }

    double **aa, *e, *et;
    int i, j;
    aa = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    for (i = 0; i < n; i++) {
        aa[i] = (double *)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
    }

    eigen_tridiagonal(aa, n, d, e, et, lt);

    if (lt > 0 && vt != NULL) {
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                vt[i][j] = aa[i][j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        free(aa[i]);
    }
    free(aa);
    free(e);
    free(et);
}

void eigen_general(double **a, double **b, int n, double *d, double **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }

    double **G, **IG, **IGA, **S, *e, *et;
    int i, j, k, ii;
    double s, ds;

    G = (double **)malloc(sizeof(double *) * n);
    IG = (double **)malloc(sizeof(double *) * n);
    IGA = (double **)malloc(sizeof(double *) * n);
    S = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);

    for (i = 0; i < n; i++) {
        G[i] = (double *)malloc(sizeof(double) * n);
        IG[i] = (double *)malloc(sizeof(double) * n);
        IGA[i] = (double *)malloc(sizeof(double) * n);
        S[i] = (double *)malloc(sizeof(double) * n);
    }

    for (j = 0; j < n; j++) {
        s = 0;
        for (k = 0; k <= j - 1; k++) {
            s = s + G[j][k] * G[j][k];
        }
        ds = b[j][j] - s;
        if (ds <= 0) {
            printf("error_cholesky\n");
            exit(1);
            return;
        }
        ds = fabs(ds);
        G[j][j] = sqrt(ds);

        for (i = j + 1; i < n; i++) {
            s = 0;
            for (k = 0; k <= j - 1; k++) {
                s = s + G[i][k] * G[j][k];
            }
            G[i][j] = (b[i][j] - s) / G[j][j];
        }
    }

    for (ii = 0; ii < n; ii++) {
        for (i = 0; i < n; i++) {
            IG[i][ii] = 0;
        }
        IG[ii][ii] = 1;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                IG[j][ii] -= G[j][i] / G[i][i] * IG[i][ii];
            }
        }
        for (i = ii; i < n; i++) {
            IG[i][ii] /= G[i][i];
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            IGA[i][j] = 0;
            for (k = 0; k <= i; k++) {
                IGA[i][j] += IG[i][k] * a[k][j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            S[i][j] = 0;
            for (k = 0; k <= j; k++) {
                S[i][j] += IGA[i][k] * IG[j][k];
            }
            S[j][i] = S[i][j];
        }
    }

    eigen_tridiagonal(S, n, d, e, et, lt);

    if (lt > 0 && vt != NULL) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < lt; i++) {
                vt[i][j] = 0;
                for (k = j; k < n; k++) {
                    vt[i][j] += S[i][k] * IG[k][j];
                }
            }
        }
    }

    for (i = 0; i < n; i++) {
        free(G[i]);
        free(IG[i]);
        free(IGA[i]);
        free(S[i]);
    }
    free(G);
    free(IG);
    free(IGA);
    free(S);
    free(e);
    free(et);
}

#ifdef LAPACKE

#include <stdio.h>
#include <lapacke.h>

void lapack_general(double **a, double **b, int n, double *e, double **vt, int lt)
{
    lapack_int N = n;          /* dimension of the matrices */
    double *A = (double *)malloc(sizeof(double) * N * N);
    double *B = (double *)malloc(sizeof(double) * N * N);
    lapack_int info;

    /* copy 2D arrays to 1D arrays in column-major order */
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            A[j * n + i] = a[i][j];
            B[j * n + i] = b[i][j];
        }
    }

    /* calculate eigenvalues and right eigenvectors */
    info = LAPACKE_dsygv(LAPACK_COL_MAJOR,     /* storage： col major */
                         1,                    /* ITYPE=1： A x = λ B x */
                         'V',                  /* calculate eigenvector 'V' / or just eigenvalue 'N' */
                         'U',                  /* use up triangular 'U' or 'L' */
                         N,
                         A, N,
                         B, N,
                         e);

    if (info == 0) {
        /* copy the eigenvectors from the 1D array to the 2D array */
        if (vt != NULL) {
             for (int i = 0; i < lt; i++) {
                for (int j = 0; j < N; j++) {
                    vt[i][j] = A[i * N + j];
                }
            }
        }
    } else if (info > 0) {
        printf("DSYGV fault，info = %d\n", info);
    } else {
        printf("Illegal argument，info = %d\n", info);
    }
}

#endif