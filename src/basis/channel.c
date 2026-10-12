/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/channel.h>

#include <stdio.h>

int channel_3b_pair_identical(int id1, int id2, int id3, int c)
{
    switch (c) {
        case 1: return id1 == id2;
        case 2: return id3 == id1;
        case 3: return id2 == id3;
        default: return 0;
    }
}

double channel_3b_exchange_eta(int f12, double sij, int lrho)
{
    int phase = ((1 + (int)sij + lrho) % 2 == 0) ? 1 : -1;
    return (double)(f12 * phase);
}

void channel_3b_ijk(double x1, double x2, double x3, double *xi, double *xj, double *xk, int c)
{
    switch (c) {
        case 1: *xi = x1; *xj = x2; *xk = x3; break;
        case 2: *xi = x3; *xj = x1; *xk = x2; break;
        case 3: *xi = x2; *xj = x3; *xk = x1; break;
        default:
            printf("error_channel_3b_ijk\n");
            *xi = 0;
            *xj = 0;
            *xk = 0;
            break;
    }
}

void channel_3b_123(double *x1, double *x2, double *x3, double xi, double xj, double xk, int c)
{
    switch (c) {
        case 1: *x1 = xi; *x2 = xj; *x3 = xk; break;
        case 2: *x1 = xj; *x2 = xk; *x3 = xi; break;
        case 3: *x1 = xk; *x2 = xi; *x3 = xj; break;
        default:
            printf("error_channel_3b_123\n");
            *x1 = 0;
            *x2 = 0;
            *x3 = 0;
            break;
    }
}
