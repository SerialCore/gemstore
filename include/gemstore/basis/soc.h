/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_SOC
#define GEMSTORE_BASIS_SOC

/* return Clebsch-Gordan coefficient for given angular momenta:
 * <j1m1j2m2|jm>
 */
double clebsch_gordan(double j1, double m1, double j2, double m2, double j, double m);

/* return 3J symbol for given angular momenta:
 * j1, j2, j
 * m1, m2, m
 */
double threeJ_symbol(double j1, double m1, double j2, double m2, double j, double m);

/* return 6J symbol for given angular momenta:
 * j1, j2, j12
 * j3,  j, j23
 */
double sixJ_symbol(double j1, double j2, double j12, double j3, double j, double j23);

/* return 9J symbol for given angular momenta:
 *  j1,  j2, j12
 *  j3,  j4, j34
 * j13, j24,   j
*/
double nineJ_symbol(double j1, double j2, double j12, double j3, double j4, double j34, double j13, double j24, double j);

/* return Kronecker delta for given arrays of indices, n is length of array */
double kronecker_delta(double *i, double *j, int n);

/* return matrix element of spin-spin operator */
double operator_sdots(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp);

/* return matrix element of spin1-orbit operator */
double operator_ldots1(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of spin2-orbit operator */
double operator_ldots2(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of tensor operator */
double operator_tensor(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

#endif