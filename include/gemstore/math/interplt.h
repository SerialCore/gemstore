/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_INTERPLT
#define GEMSTORE_MATH_INTERPLT

#include <gemstore/math/matrix.h>

/* use quadratic interpolation to fix anomalies */
void interpolate_quadratic(array_t *data);

#endif