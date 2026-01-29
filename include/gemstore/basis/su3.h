/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_SU3
#define GEMSTORE_BASIS_SU3

/* Calculate dimension of SU(3) representation */
int su3_dimension(int upper, int lower);

/* Product SU(3) representations and find sums */
void su3_product(int a_upper, int a_lower, int b_upper, int b_lower);

#endif