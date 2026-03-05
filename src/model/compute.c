/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/model/spectra.h>
#include <gemstore/math/matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compute_spectra_meson(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, const char *model)
{
    array_t eigenvalue = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    argsModelDy_t args_dynmc = {0};

    if (strcmp(model, "GIString") == 0) {
        args_dynmc.model = MODEL_GI_STRING;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(f1, f2, S, L, J, nmax, rmax, rmin, &argsGIString_meson, &args_dynmc, &eigenvalue, &eigenvector, nmax);

    } else if (strcmp(model, "GIScreen") == 0) {
        args_dynmc.model = MODEL_GI_SCREEN;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(f1, f2, S, L, J, nmax, rmax, rmin, &argsGIScreen_meson, &args_dynmc, &eigenvalue, &eigenvector, nmax);

    } else if (strcmp(model, "GIQuadra") == 0) {
        args_dynmc.model = MODEL_GI_QUADRA;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_GI(f1, f2, S, L, J, nmax, rmax, rmin, &argsGIQuadra_meson, &args_dynmc, &eigenvalue, &eigenvector, nmax);

    } else {
        printf("Unknown model: %s\n", model);
        array_free(&eigenvalue);
        matrix_free(&eigenvector);
        return;
    }

    array_print(&eigenvalue);
    matrix_print(&eigenvector);

    array_free(&eigenvalue);
    matrix_free(&eigenvector);
}