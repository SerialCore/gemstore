/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_MESON
#define GEMSTORE_PARAM_MESON

/* Make sure c program can only see this c++ entry function */
#ifdef __cplusplus
extern "C" {
#endif

void minuit2_meson_GIScreen(double *params_out);

void minuit2_meson_GIQuadra(double *params_out);

#ifdef __cplusplus
}
#endif

#endif