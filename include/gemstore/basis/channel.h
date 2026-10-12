/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Three-body channel labels: which quarks form the pair, whether that
 * pair is identical, and the 1↔2 exchange phase. These are index
 * permutations, not the linear Jacobi maps in jacobi.h.
 * Four-body channels would be channel_4b_*.
 */

#ifndef GEMSTORE_BASIS_CHANNEL
#define GEMSTORE_BASIS_CHANNEL

/* Channel c: 1 = pair(1,2) spectator 3,
 *            2 = pair(3,1) spectator 2,
 *            3 = pair(2,3) spectator 1.
 */

int channel_3b_pair_identical(int id1, int id2, int id3, int c);

/* 1↔2 phase: f12 * (-1)^{1+sij+lρ}. */
double channel_3b_exchange_eta(int f12, double sij, int lrho);

/* Permute a triple onto the pair of channel c, and back.
 * c=1 (x1,x2,x3), c=2 (x3,x1,x2), c=3 (x2,x3,x1). */
void channel_3b_ijk(double x1, double x2, double x3, double *xi, double *xj, double *xk, int c);
void channel_3b_123(double *x1, double *x2, double *x3, double xi, double xj, double xk, int c);

#endif
