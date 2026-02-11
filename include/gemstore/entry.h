/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_ENTRY
#define GEMSTORE_ENTRY

void call_minuit2();

void call_spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin);

/* C entry function to be called by c++ */
#ifdef __cplusplus
extern "C" {
#endif

void call_fitting_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, double *e_out, double **v_out);

#ifdef __cplusplus
}
#endif

#endif