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
    double *e = (double *)malloc(nmax * sizeof(double));
    array_t eigenvalue = {
        .con = nmax,
        .value = e
    };

    spectra_meson_NRScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, NULL, 0, NULL);
    array_print(&eigenvalue);

    array_free(&eigenvalue);
}

double call_fitting_meson_NRScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params)
{
    double *e = (double *)malloc(nmax * sizeof(double));
    array_t eigenvalue = {
        .con = nmax,
        .value = e
    };

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
    double e_out = e[N - 1];

    array_free(&eigenvalue);
    return e_out;
}

void call_minuit2_chi2()
{
    double *params = (double *)malloc(9 * sizeof(double));
    perform_fit(params);
}