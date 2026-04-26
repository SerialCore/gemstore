/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>
#include <gemstore/param/argset.h>
#include <gemstore/param/meson.h>
#include <gemstore/param/bbbar.h>
#include <gemstore/param/ccbar.h>

#include <gemstore/math/matrix.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/model/spectra.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double call_meson_GIScreen(int f1, int f2, int N, double S, double L, double J, int nmax, double rmax, double rmin, const double *params)
{
    array_t eigenvalue = array_init(nmax);

    argsInput_t args_input = {
        .f1 = f1,
        .f2 = f2,
        .S = S,
        .L = L,
        .J = J,
        .nmax = nmax,
        .rmax = rmax,
        .rmin = rmin
    };
    argsGIModel_t args_model = {
        .mn = params[0],
        .ms = params[1],
        .mc = params[2],
        .mb = params[3],
        .b = params[4],
        .mu = params[5],
        .c = params[6],
        .sigma_0 = params[7],
        .s = params[8],
        .epsilon_cont = params[9],
        .epsilon_sov = params[10],
        .epsilon_sos = params[11],
        .epsilon_tens = params[12]
    };
    argsGIModelDy_t args_dynmc = {
        .model = MODEL_GISCREEN,
        .system = SYSTEM_MESON
    };

    spectra_meson_GEM(&args_input, &args_model, &args_dynmc, &eigenvalue, NULL, 0);
    double e_out = eigenvalue.value[N - 1];

    array_free(&eigenvalue);
    return e_out;
}

void call_minuit2_GIScreen(const char* system)
{
    double *params = (double *)malloc(13 * sizeof(double));
    if (strcmp(system, "meson") == 0) minuit2_meson_GIScreen(params);
    else if (strcmp(system, "bbbar") == 0) minuit2_bbbar_GIScreen(params);
    else if (strcmp(system, "ccbar") == 0) minuit2_ccbar_GIScreen(params);
    else {fprintf(stderr, "Unknown fitting system: %s\n", system); exit(1);}

    free(params);
}