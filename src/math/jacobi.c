/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/jacobi.h>

void jacobi_3b_r_1(double m1, double m2, double m3, int a, int c,
    double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac)
{
    double alpha[3][3] = {{1, -(m3 / (m1 + m3)), -(m3 / (m2 + m3))}, {-(m2 / (m1 + m2)), 1, -(m2 / (m2 + m3))}, {-(m1 / (m1 + m2)), -(m1 / (m1 + m3)), 1}};
    double beta[3][3] = {{0, 1, -1}, {-1, 0, 1}, {1, -1, 0}};
    double gamma[3][3] = {{0, -1 + (m2 * m3) / ((m1 + m2) * (m1 + m3)), 1 - (m1 * m3) / ((m1 + m2) * (m2 + m3))}, {1 - (m2 * m3) / ((m1 + m2) * (m1 + m3)), 0, -1 + (m1 * m2) / ((m1 + m3) * (m2 + m3))}, {-1 + (m1 * m3) / ((m1 + m2) * (m2 + m3)), 1 - (m1 * m2) / ((m1 + m3) * (m2 + m3)), 0}};
    double delta[3][3] = {{1, -(m2 / (m1 + m2)), -(m1 / (m1 + m2))}, {-(m3 / (m1 + m3)), 1, -(m1 / (m1 + m3))}, {-(m3 / (m2 + m3)), -(m2 / (m2 + m3)), 1}};
    *alpha_ac = alpha[a - 1][c - 1];
    *beta_ac = beta[a - 1][c - 1];
    *gamma_ac = gamma[a - 1][c - 1];
    *delta_ac = delta[a - 1][c - 1];
}

void jacobi_3b_r_2(double m1, double m2, double m3, int a, int c,
    double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac)
{
    double alpha[3][3] = {{0, 1, -1}, {-1, 0, 1}, {1, -1, 0}};
    double beta[3][3] = {{1, -(m3 / (m1 + m3)), -(m3 / (m2 + m3))}, {-(m2 / (m1 + m2)), 1, -(m2 / (m2 + m3))}, {-(m1 / (m1 + m2)), -(m1 / (m1 + m3)), 1}};
    double gamma[3][3] = {{1, -(m2 / (m1 + m2)), -(m1 / (m1 + m2))}, {-(m3 / (m1 + m3)), 1, -(m1 / (m1 + m3))}, {-(m3 / (m2 + m3)), -(m2 / (m2 + m3)), 1}};
    double delta[3][3] = {{0, -1 + (m2 * m3) / ((m1 + m2) * (m1 + m3)), 1 - (m1 * m3) / ((m1 + m2) * (m2 + m3))}, {1 - (m2 * m3) / ((m1 + m2) * (m1 + m3)), 0, -1 + (m1 * m2) / ((m1 + m3) * (m2 + m3))}, {-1 + (m1 * m3) / ((m1 + m2) * (m2 + m3)), 1 - (m1 * m2) / ((m1 + m3) * (m2 + m3)), 0}};
    *alpha_ac = alpha[a - 1][c - 1];
    *beta_ac = beta[a - 1][c - 1];
    *gamma_ac = gamma[a - 1][c - 1];
    *delta_ac = delta[a - 1][c - 1];
}

void jacobi_3b_p_1(double m1, double m2, double m3, int a, int c,
    double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac)
{
    double alpha[3][3] = {{1, -(m2 / (m1 + m2)), -(m1 / (m1 + m2))}, {-(m3 / (m1 + m3)), 1, -(m1 / (m1 + m3))}, {-(m3 / (m2 + m3)), -(m2 / (m2 + m3)), 1}};
    double beta[3][3] = {{0, (m1 * (m1 + m2 + m3)) / ((m1 + m2) * (m1 + m3)), -((m2 * (m1 + m2 + m3)) / ((m1 + m2) * (m2 + m3)))}, {-((m1 * (m1 + m2 + m3)) / ((m1 + m2) * (m1 + m3))), 0, (m3 * (m1 + m2 + m3)) / ((m1 + m3) * (m2 + m3))}, {(m2 * (m1 + m2 + m3)) / ((m1 + m2) * (m2 + m3)), -((m3 * (m1 + m2 + m3)) / ((m1 + m3) * (m2 + m3))), 0}};
    double gamma[3][3] = {{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}};
    double delta[3][3] = {{1, -(m3 / (m1 + m3)), -(m3 / (m2 + m3))}, {-(m2 / (m1 + m2)), 1, -(m2 / (m2 + m3))}, {-(m1 / (m1 + m2)), -(m1 / (m1 + m3)), 1}};
    *alpha_ac = alpha[a - 1][c - 1];
    *beta_ac = beta[a - 1][c - 1];
    *gamma_ac = gamma[a - 1][c - 1];
    *delta_ac = delta[a - 1][c - 1];
}

void jacobi_3b_p_2(double m1, double m2, double m3, int a, int c,
    double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac)
{
    double alpha[3][3] = {{0, (m1 * (m1 + m2 + m3)) / ((m1 + m2) * (m1 + m3)), -((m2 * (m1 + m2 + m3)) / ((m1 + m2) * (m2 + m3)))}, {-((m1 * (m1 + m2 + m3)) / ((m1 + m2) * (m1 + m3))), 0, (m3 * (m1 + m2 + m3)) / ((m1 + m3) * (m2 + m3))}, {(m2 * (m1 + m2 + m3)) / ((m1 + m2) * (m2 + m3)), -((m3 * (m1 + m2 + m3)) / ((m1 + m3) * (m2 + m3))), 0}};
    double beta[3][3] = {{1, -(m2 / (m1 + m2)), -(m1 / (m1 + m2))}, {-(m3 / (m1 + m3)), 1, -(m1 / (m1 + m3))}, {-(m3 / (m2 + m3)), -(m2 / (m2 + m3)), 1}};
    double gamma[3][3] = {{1, -(m3 / (m1 + m3)), -(m3 / (m2 + m3))}, {-(m2 / (m1 + m2)), 1, -(m2 / (m2 + m3))}, {-(m1 / (m1 + m2)), -(m1 / (m1 + m3)), 1}};
    double delta[3][3] = {{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}};
    *alpha_ac = alpha[a - 1][c - 1];
    *beta_ac = beta[a - 1][c - 1];
    *gamma_ac = gamma[a - 1][c - 1];
    *delta_ac = delta[a - 1][c - 1];
}

void jacobi_3b_p_i(double m1, double m2, double m3, int a, int c,
    double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac)
{
    double alpha[3][3] = {{m2 / (m1 + m2), -1, m1 / (m1 + m2)}, {-1, m3 / (m1 + m3), m1 / (m1 + m3)}, {m2 / (m2 + m3), m3 / (m2 + m3), -1}};
    double beta[3][3] = {{-(m1 / (m1 + m2)), -(m2 / (m1 + m2)), m1 / (m1 + m2) + m2 / (m1 + m2)}, {-(m1 / (m1 + m3)), m1 / (m1 + m3) + m3 / (m1 + m3), -(m3 / (m1 + m3))}, {m2 / (m2 + m3) + m3 / (m2 + m3), -(m2 / (m2 + m3)), -(m3 / (m2 + m3))}};
    double gamma[3][3] = {{1, 0, -1}, {0, -1, 1}, {-1, 1, 0}};
    double delta[3][3] = {{1, -1, 0}, {-1, 0, 1}, {0, 1, -1}};
    *alpha_ac = alpha[a - 1][c - 1];
    *beta_ac = beta[a - 1][c - 1];
    *gamma_ac = gamma[a - 1][c - 1];
    *delta_ac = delta[a - 1][c - 1];
}
