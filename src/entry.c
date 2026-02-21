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

#include <stdio.h>
#include <stdlib.h>

void call_spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin)
{
    array_t eigenvalue = array_init(nmax);

    spectra_meson_NRScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, NULL, 0, NULL);
    array_print(&eigenvalue);

    array_free(&eigenvalue);
}

void call_spectra_meson_GIScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin)
{
    array_t eigenvalue = array_init(nmax);

    spectra_meson_GIScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, NULL, 0, NULL);
    array_print(&eigenvalue);

    array_free(&eigenvalue);
}

double call_fitting_meson_NRScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params)
{
    array_t eigenvalue = array_init(nmax);

    argsModel_t args_model = {
        .mn = params[0],
        .ms = params[1],
        .mc = params[2],
        .mb = params[3],
        .alpha_s = params[4],
        .b1 = params[5],
        .mu = params[6],
        .c = params[7],
        .sigma = params[8]
    };

    spectra_meson_NRScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, NULL, 0, &args_model);
    double e_out = eigenvalue.value[N - 1];

    array_free(&eigenvalue);
    return e_out;
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

    spectra_meson_GIScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, NULL, 0, &args_model);
    double e_out = eigenvalue.value[N - 1];

    array_free(&eigenvalue);
    return e_out;
}

void call_minuit2_chi2()
{
    double *params = (double *)malloc(20 * sizeof(double));
    perform_fit(params);
}