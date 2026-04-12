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

    /* prepare variables */
    argsOrbit_t args_bra;
    argsOrbit_t args_ket;
    double factor;
    double coef;
    double r2sum;
    double oversum;
    double fm = 5.06773093854369882649;

    /* construct matrices */
    matrix_t mR2 = matrix_init(nmax, nmax);
    matrix_t mOver = matrix_init(nmax, nmax);

    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            args_bra = basis[i];
            args_ket = basis[j];

            factor = 1.0 / sqrt(args_bra.scale + args_ket.scale);
            mR2.value[i][j] = integral_rms_radius(GRnlr, factor, &args_bra, &args_ket);
            mOver.value[i][j] = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_ket);
        }
    }

    /* calculate RMS radius */
    for (int n = 0; n < len; n++) {
        r2sum = 0.0;
        oversum = 0.0;
        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                coef = vector->value[n][i] * vector->value[n][j];
                r2sum += coef * mR2.value[i][j];
                oversum += coef * mOver.value[i][j];
            }
        }
        radius->value[n] = sqrt(r2sum / oversum) / fm;
    }

    free(basis);
    matrix_free(&mR2);
    matrix_free(&mOver);
}
