/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef DEBUG_H
#define DEBUG_H

void getq(double mn, double ms, double mc, double mb, int q, double *mi);

void debug01(int q1, int q2, int q3, int f12, double J, int P, int Lmax);

void debug02(int q1, int q2, int q3, int f12);

#endif