/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_ENTRY
#define GEMSTORE_ENTRY

/* C entry function to be called by c++ */
#ifdef __cplusplus
extern "C" {
#endif

double call_fitting_meson_NRScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params);

double call_fitting_meson_GIScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params);

#ifdef __cplusplus
}
#endif

void call_spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin);

void call_spectra_meson_GIScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin);

void call_minuit2_chi2();

#endif