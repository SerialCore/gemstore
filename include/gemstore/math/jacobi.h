/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_JACOBI
#define GEMSTORE_MATH_JACOBI

/* Jacobi channel c: 1 = pair(1,2) spectator 3,
 *                   2 = pair(3,1) spectator 2,
 *                   3 = pair(2,3) spectator 1.
 *
 * Map basis channel a onto operator channel c:
 *   ρ_a = α ρ_c + β λ_c
 *   λ_a = γ ρ_c + δ λ_c
 * jacobi_3b_r_1 maps along ρ. jacobi_3b_r_2 swaps the roles of ρ and λ.
 * jacobi_3b_p_1, jacobi_3b_p_2 and jacobi_3b_p_i are the momentum maps for ρ, λ and one quark.
 */

void jacobi_3b_r_1(double m1, double m2, double m3, int a, int c,
    double *alpha, double *beta, double *gamma, double *delta);

void jacobi_3b_r_2(double m1, double m2, double m3, int a, int c,
    double *alpha, double *beta, double *gamma, double *delta);

void jacobi_3b_p_1(double m1, double m2, double m3, int a, int c,
    double *alpha, double *beta, double *gamma, double *delta);

void jacobi_3b_p_2(double m1, double m2, double m3, int a, int c,
    double *alpha, double *beta, double *gamma, double *delta);

void jacobi_3b_p_i(double m1, double m2, double m3, int a, int c,
    double *alpha, double *beta, double *gamma, double *delta);

#endif
