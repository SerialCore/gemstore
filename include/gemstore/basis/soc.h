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

/* define operator function for sl coupling */
typedef double (*operator_sl)(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of center potential in sl couping*/
double operator_center_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of spin-spin operator in sl couping */
double operator_sdots_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of spin1-orbit operator in sl couping */
double operator_ldots1_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of spin2-orbit operator in sl couping */
double operator_ldots2_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* return matrix element of tensor operator in sl couping */
double operator_tensor_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j);

/* define operator function for sl coupling */
typedef double (*operator_jj)(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

/* return matrix element of center potential in jj couping*/
double operator_center_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

/* return matrix element of spin-spin operator in jj couping */
double operator_sdots_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

/* return matrix element of spin1-orbit operator in jj couping */
double operator_ldots1_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

/* return matrix element of spin2-orbit operator in jj couping */
double operator_ldots2_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

/* return matrix element of tensor operator in jj couping */
double operator_tensor_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j);

#endif