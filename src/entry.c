/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/fileio.h>
#include <gemstore/fitting.h>

#include <gemstore/numerical/spectra.h>
#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/model.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void call_spectra_meson(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, const char *model)
{
    array_t eigenvalue = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    argsModelDy_t args_dynmc = {0};

    if (strcmp(model, "NRScreen") == 0) {
        args_dynmc.model = MODEL_NR_SCREEN;
        args_dynmc.system = SYSTEM_MESON;
        spectra_meson_NR(f1, f2, S, L, J, nmax, rmax, rmin, &argsNRScreen_meson, &args_dynmc, &eigenvalue, &eigenvector, nmax);

    } else if (strcmp(model, "GIString") == 0) {
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

double call_fitting_meson_GIScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params)
{
    array_t eigenvalue = array_init(nmax);

    argsModel_t args_model = {
        .mn = params[0],
        .ms = params[1],
        .mc = params[2],
        .mb = params[3],
        .mt = 172.57,
        .b1 = params[4],
        .mu = params[5],
        .c = params[6],
        .sigma_0 = params[7],
        .s = params[8],
        .epsilon_Coul = 0.0,
        .epsilon_cont = params[9],
        .epsilon_sov = params[10],
        .epsilon_sos = params[11],
        .epsilon_tens = params[12],
    };
    argsModelDy_t args_dynmc = {
        .model = MODEL_GI_SCREEN,
        .system = SYSTEM_MESON
    };

    spectra_meson_GI(f1, f2, S, L, J, nmax, rmax, rmin, &args_model, &args_dynmc, &eigenvalue, NULL, 0);
    double e_out = eigenvalue.value[N - 1];

    array_free(&eigenvalue);
    return e_out;
}

double call_fitting_meson_GIQuadra(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params)
{
    array_t eigenvalue = array_init(nmax);

    argsModel_t args_model = {
        .mn = params[0],
        .ms = params[1],
        .mc = params[2],
        .mb = params[3],
        .mt = 172.57,
        .b1 = params[4],
        .b2 = params[5],
        .mu = params[6],
        .c = params[7],
        .sigma_0 = params[8],
        .s = params[9],
        .epsilon_Coul = 0.0,
        .epsilon_cont = params[10],
        .epsilon_sov = params[11],
        .epsilon_sos = params[12],
        .epsilon_tens = params[13],
    };
    argsModelDy_t args_dynmc = {
        .model = MODEL_GI_QUADRA,
        .system = SYSTEM_MESON
    };

    spectra_meson_GI(f1, f2, S, L, J, nmax, rmax, rmin, &args_model, &args_dynmc, &eigenvalue, NULL, 0);
    double e_out = eigenvalue.value[N - 1];

    array_free(&eigenvalue);
    return e_out;
}

void call_minuit2_chi2()
{
    double *params = (double *)malloc(20 * sizeof(double));
    perform_fit(params);
}