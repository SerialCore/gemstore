/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/su3.h>

#include <stdio.h>

int su3_dimension(int upper, int lower)
{
    return (upper + 1) * (lower + 1) * (upper + lower + 2) / 2;
}

/* Print representation details */
static void print_rep(int upper, int lower);
static void print_rep(int upper, int lower)
{
    printf("{%d, %d} dim=%d\n", upper, lower, su3_dimension(upper, lower));
}

/* Find direct sums */
static void indice(int n, int nprime, int m, int mprime);
static void indice(int n, int nprime, int m, int mprime)
{
    print_rep(n + nprime, m + mprime);
    
    int min_p = (n < nprime) ? n : nprime;
    for (int p = 1; p <= min_p; p++) {
        print_rep(n + nprime - 2 * p, m + mprime + p);
    }
    
    int min_q = (m < mprime) ? m : mprime;
    for (int q = 1; q <= min_q; q++) {
        print_rep(n + nprime + q, m + mprime - 2 * q);
    }
}

void su3_product(int a_upper, int a_lower, int b_upper, int b_lower)
{
    printf("\"direct product\"\n");
    print_rep(a_upper, a_lower);
    print_rep(b_upper, b_lower);

    printf("\"direct sum\"\n");
    int min_i = (a_upper < b_lower) ? a_upper : b_lower;
    for (int i = 0; i <= min_i; i++) {
        int min_j = (a_lower < b_upper) ? a_lower : b_upper;
        for (int j = 0; j <= min_j; j++) {
            indice(a_upper - i, b_upper - j, a_lower - j, b_lower - i);
        }
    }
}