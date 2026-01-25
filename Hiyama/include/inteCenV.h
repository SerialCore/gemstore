/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INTECENV_H
#define INTECENV_H

#include <basis.h>
#include <sumckdk.h>

typedef struct
{
	double mq;
	double ms;
	double mc;
	double mb;
	double b;
	double K;
	double alpha_ss;
	double alpha_so;
	double alpha_ten;
	double C;
	double Lambda;
} vargs;

typedef double (*vCenPart)(double, basis_qnum, basis_qnum, vargs, int);

double vNfi(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vConLine(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vOgeCoul(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vOgeCont(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vSorp(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vSorn(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vTens(double r, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vTi(double p, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vT1(double p, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double vT2(double p, basis_qnum qf, basis_qnum qi, vargs varg, int c);

double InteCenV(double b11, int n, basis_qnum qf, basis_qnum qi, vCenPart v, vargs varg, int c);

double Nnl(int l, double nu);

void coordinatesTransformation_r_1(double m1, double m2, double m3, int a, int c, double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac);

void coordinatesTransformation_r_2(double m1, double m2, double m3, int a, int c, double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac);

void coordinatesTransformation_p_1(double m1, double m2, double m3, int a, int c, double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac);

void coordinatesTransformation_p_2(double m1, double m2, double m3, int a, int c, double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac);

void coordinatesTransformation_p_i(double m1, double m2, double m3, int a, int c, double *alpha_ac, double *beta_ac, double *gamma_ac, double *delta_ac);

void getT1r(basis_qnum qf, basis_qnum qi, double *txrp_cent, int c, int vt);

void getT2r(basis_qnum qf, basis_qnum qi, double *txrp_cent, int c, int vt);

void getT1p(basis_qnum qf, basis_qnum qi, double *txrp_cent, int c, int vt);

void getT2p(basis_qnum qf, basis_qnum qi, double *txrp_cent, int c, int vt);

void getTpi(basis_qnum qf, basis_qnum qi, double *txrp_cent, int c, int vt);

typedef void (*txrp_vtype)(basis_qnum, basis_qnum, double *, int *, int);

void t1r_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t2r_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1p_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t2p_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void tpi_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_tens(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_soii(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_soij(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_soji(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_sojj(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_sorp(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

void t1r_sorn(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

double inteVcenPart(txrp_vtype get_trxp_vtype, sumckdk_scdk scdk, basis_qnum qf, basis_qnum qi, vargs varg, vCenPart v, int c);

#endif