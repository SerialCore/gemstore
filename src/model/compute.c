/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/spectra.h>
#include <gemstore/model/radius.h>
#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

#include <gemstore/types.h>
#include <gemstore/print.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void compute_spectra_meson(const argsInput_t *input)
{
    int nmax = input->nmax;

    array_t eigenvalue = array_init(nmax);
    array_t rmsradius = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    argsGIModelDy_t args_dynmc = {0};

    if (input->model == MODEL_GI_STRING) {
        args_dynmc.model = MODEL_GI_STRING;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(input, &input->params, &args_dynmc, &eigenvalue, &eigenvector, nmax);
    }
    else if (input->model == MODEL_GI_SCREEN) {
        args_dynmc.model = MODEL_GI_SCREEN;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(input, &input->params, &args_dynmc, &eigenvalue, &eigenvector, nmax);
    }
    else if (input->model == MODEL_GI_QUADRA) {
        args_dynmc.model = MODEL_GI_QUADRA;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(input, &input->params, &args_dynmc, &eigenvalue, &eigenvector, nmax);
    }
    else {
        array_free(&eigenvalue);
        array_free(&rmsradius);
        matrix_free(&eigenvector);
        return;
    }

    radius_meson_rms(input, &eigenvector, &rmsradius, nmax);

    /* debug the results */
    for (int n = 0; n < nmax; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        for (int i = 0; i < nmax; i++) {
            double c = fabs(eigenvector.value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }
        printf("State %2d:  mass=%2.6f  RMS=%2.3f  max|c|=%.3f  ||c||^2=%.10f\n", 
            n+1, eigenvalue.value[n], rmsradius.value[n], maxc, norm);
    }

    write_meson_spectra(input, &eigenvalue, &rmsradius, &eigenvector, nmax);

    array_free(&eigenvalue);
    array_free(&rmsradius);
    matrix_free(&eigenvector);
}