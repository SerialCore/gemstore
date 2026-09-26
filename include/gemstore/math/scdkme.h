/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * SCDK kinematics, operator types, and radial integrals (recycle/inteCenV.h, vtype.h).
 */

#ifndef GEMSTORE_MATH_SCDKME
#define GEMSTORE_MATH_SCDKME

#include <gemstore/basis/basis.h>
#include <gemstore/math/sumckdk.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>
#include <gemstore/types.h>

#include <pthread.h>

typedef struct {
    double mn, ms, mc, mb;
    double alpha[3];
    double gamma[3];
    double b, c, sigma0, s, f, mu;
    double econt, etens, esov, esos, eCoul;
    model_type_t model;
} scdk_vargs_t;

void scdk_vargs_from_model(scdk_vargs_t *varg, const argsGIModel_t *model, model_type_t mtype);

void getijk(double x1, double x2, double x3, double *xi, double *xj, double *xk, int c);
void get123(double *x1, double *x2, double *x3, double xi, double xj, double xk, int c);

typedef void (*txrp_vtype)(basis_qnum, basis_qnum, double *, int *, int);
typedef double (*scdk_inte_fn)(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);

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
void t1r_sorr(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void t1r_sorl(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tir_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tir_tens(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tir_soii(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tir_soij(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tjr_cent(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tjr_tens(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tjr_sojj(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);
void tjr_soji(basis_qnum qf, basis_qnum qi, double *t, int *lent, int c);

double inteVcenPartA(txrp_vtype get_trxp_vtype, sumckdk_scdk scdk, basis_qnum qf, basis_qnum qi,
    scdk_vargs_t varg, scdk_inte_fn iv, int c);

double inteNfi(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVogeG(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVcont(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVtens(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsovii(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsovjj(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsovij(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsovji(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVstring(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsosii(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteVsosjj(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepogeG(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepcont(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteptens(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsovii(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsovjj(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsovji(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsovij(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsosii(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double intepsosjj(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteTi(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);
double inteRMS(double b11, int n, basis_qnum qf, basis_qnum qi, scdk_vargs_t varg, int c);

typedef void (*vtype)(sumckdk_scdk *, double, double, double, double, double, double, double, double,
    int, int, int, int, int, int, int, int, basis_qnum, basis_qnum, pthread_mutex_t *, int, int);

void sumckdk_scdk_vtype(sumckdk_scdk *scdk, vtype vsodt, int fg, basis_list qnlist, matrix_t **mlsj,
    int nf, int nfp, int ni, int nip, pthread_mutex_t *mutex, int lock, int c);

void vcent(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vcont(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vtens(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vsoii(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vsoji(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vsojj(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);
void vsoij(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

#endif
