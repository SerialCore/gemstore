/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_INTERPLT
#define GEMSTORE_MATH_INTERPLT

/* use quadratic interpolation to fix anomalies */
void interpolate_quadratic(double *data, int n);

#endif