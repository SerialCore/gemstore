/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_FITCCBAR
#define GEMSTORE_PARAM_FITCCBAR

/* Make sure c program can only see this c++ entry function */
#ifdef __cplusplus
extern "C" {
#endif

void minuit2_ccbar_GIScreen(double *params_out);

void minuit2_ccbar_GIQuadra(double *params_out);

#ifdef __cplusplus
}
#endif

#endif