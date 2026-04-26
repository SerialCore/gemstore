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
    double factor;
    double coef;
    double r2sum;
    double oversum;
    double norm;
    double rms2;
    double fm = 5.06773093854369882649;

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);
            mR2.value[i][j] = integral_wfn_radius(GRnlr, factor, &basis[i], &basis[j]);
            mOver.value[i][j] = integral_wfn_overlap(GRnlr, factor, &basis[i], &basis[j]);
        }
    }

    /* calculate rms radius with orthogonalized coefficients */
    for (int n = 0; n < len; n++) {
        r2sum = 0.0;
        oversum = 0.0;
        norm = 0.0;

        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                coef = vector->value[n][i] * vector->value[n][j];
                r2sum += coef * mR2.value[i][j];
                oversum += coef * mOver.value[i][j];
            }
            norm += vector->value[n][i] * vector->value[n][i];
        }

        rms2 = (norm > 1e-12) ? r2sum / norm : 0.0;
        radius->value[n] = sqrt(rms2) / fm;
    }

    free(basis);
    matrix_free(&mR2);
    matrix_free(&mOver);
}