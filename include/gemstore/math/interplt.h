/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_INTERPLT
#define GEMSTORE_MATH_INTERPLT

/* fix divergent anomalies with cubic interpolation */
void interpolate_divergence(double *data, int n);

/* fix divergent anomalies with slop-based extrapolation */
void extrapolate_divergence(double *data, int n);

#endif
