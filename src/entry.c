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
    double *e1 = (double *)malloc(nmax * sizeof(double));
    double **v = (double **)malloc(nmax * sizeof(double *));
    for (int i = 0; i < nmax; i++) {
        v[i] = (double *)malloc(nmax * sizeof(double));
    }

    array_t eigenvalue = {
        .con = nmax,
        .value = e1
    };

    matrix_t eigenvector = {
        .row = nmax,
        .col = nmax,
        .value = v
    };

    spectra_meson_NRScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, &eigenvector, 0);

    array_print(&eigenvalue);

    array_free(&eigenvalue);
    matrix_free(&eigenvector);
}

void call_fitting_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, double *e_out, double **v_out)
{
    array_t eigenvalue = {
        .con = nmax,
        .value = e_out
    };

    matrix_t eigenvector = {
        .row = nmax,
        .col = nmax,
        .value = v_out
    };

    spectra_meson_NRScreen(f1, f2, S, L, J, nmax, rmax, rmin, &eigenvalue, &eigenvector, 0);
}

void call_minuit2()
{
    double a[] = {1.1, 1.2, 1.3};
    test(a);
}