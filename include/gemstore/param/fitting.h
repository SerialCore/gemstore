/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_FITTING
#define GEMSTORE_PARAM_FITTING

/* C entry function to be called by c++ */
#ifdef __cplusplus
extern "C" {
#endif

double call_meson_GIScreen(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params);

double call_meson_GIQuadra(int f1, int f2, int N, int S, int L, int J, int nmax, double rmax, double rmin, const double *params);

#ifdef __cplusplus
}
#endif

void call_minuit2_GIScreen(const char* system);

void call_minuit2_GIQuadra(const char* system);

#endif