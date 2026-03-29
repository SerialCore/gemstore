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

#include <math.h>
#include <stdlib.h>

void radius_meson_rms(const argsInput_t *args_input, const matrix_t *vector, array_t *r_out, int v_len)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    int L = (int)args_input->L;

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
    double overlap;
    double coef;
    double r2;
    double sum;
    double fm = 5.06773093854369882649;

    for (int n = 0; n < v_len; n++) {
        sum = 0.0;
        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                args_bra = basis[i];
                args_ket = basis[j];

                coef = vector->value[n][i] * vector->value[n][j];
                factor = 1.0 / sqrt(args_bra.scale + args_ket.scale);

                overlap = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_ket);
                r2 = integral_rms_radius(GRnlr, factor, &args_bra, &args_ket);
                sum += coef * r2 / overlap;
            }
        }
        r_out->value[n] = sqrt(sum) / fm;
    }

    free(basis);
}