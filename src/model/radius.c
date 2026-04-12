/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/radius.h>

#include <gemstore/basis/orbit.h>
#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void radius_meson_rms(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len)
{
    int nmax = input->nmax;
    double rmax = input->rmax;
    double rmin = input->rmin;
    int L = (int)input->L;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
    }

    /* construct matrices */
    matrix_t mR2 = matrix_init(nmax, nmax);
    matrix_t mOver = matrix_init(nmax, nmax);

    /* prepare variables */
    argsOrbit_t args_bra;
    argsOrbit_t args_ket;
    double factor;
    double r2sum;
    double oversum;
    double fm = 5.06773093854369882649;

    /* precompute matrices */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            args_bra = basis[i];
            args_ket = basis[j];

            factor = 1.0 / sqrt(args_bra.scale + args_ket.scale);
            mR2.value[i][j] = integral_rms_radius(GRnlr, factor, &args_bra, &args_ket);
            mOver.value[i][j] = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_ket);
        }
    }

    /* Cholesky decomposition, construct matrix R' of r^2 under orthogonal basis */
    matrix_t Lm = matrix_init(nmax, nmax);
    matrix_t Linv = matrix_init(nmax, nmax);
    matrix_cholesky_decomp(&mOver, &Lm);
    matrix_inverse_lowertri(&Lm, &Linv);

    matrix_t temp = matrix_init(nmax, nmax);
    matrix_t LinvT = matrix_init(nmax, nmax);
    matrix_t Rprime = matrix_init(nmax, nmax);
    matrix_product(&Linv, &mR2, &temp);             /* Linv * R */
    matrix_transpose(&Linv, &LinvT);                /* Linv^T = L^{-T} */
    matrix_product(&temp, &LinvT, &Rprime);         /* R' = Linv * R * Linv^T */

    /* Free temporary matrices after transformation (before loop) */
    matrix_free(&Linv);
    matrix_free(&temp);
    matrix_free(&LinvT);

    /* calculate rms radius with orthogonal coefficients d */
    for (int n = 0; n < len; n++) {
        double *d_vec = (double *)calloc(nmax, sizeof(double));
        r2sum = 0.0;
        oversum = 0.0;

        /* d = L^T c */
        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                d_vec[i] += Lm.value[j][i] * vector->value[n][j];   /* L^T[i][j] = L[j][i] */
            }
            oversum += d_vec[i] * d_vec[i];
        }

        /* r2sum = d^T R' d */
        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                r2sum += d_vec[i] * Rprime.value[i][j] * d_vec[j];
            }
        }

        radius->value[n] = (oversum > 1e-12) ? sqrt(r2sum / oversum) / fm : 0.0;

        double maxd = 0.0;
        for (int i = 0; i < nmax; i++) if (fabs(d_vec[i]) > maxd) maxd = fabs(d_vec[i]);
        printf("State %2d: RMS=%.6f  max|d|=%.3f  oversum=%.12f\n", n+1, radius->value[n], maxd, oversum);

        free(d_vec);
    }

    /* Free matrices after loop completes */
    matrix_free(&Lm);
    matrix_free(&Rprime);

    free(basis);
    matrix_free(&mR2);
    matrix_free(&mOver);
}