/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_FITTING
#define GEMSTORE_FITTING

/* Make sure c program can only see this c++ entry function */
#ifdef __cplusplus
extern "C" {
#endif

void perform_fit(double *params_out);

#ifdef __cplusplus
}
#endif

#endif