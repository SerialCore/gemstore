/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef VTYPE_H
#define VTYPE_H

#include <basis.h>
#include <matrix.h>
#include <sumckdk.h>

#include <unistd.h>
#include <pthread.h>

typedef void (*vtype)(sumckdk_scdk *, double, double, double, double, double, double, double, double, int, int, int, int, int, int, int, int, basis_qnum, basis_qnum, pthread_mutex_t *, int, int);

void getijk(double x1, double x2, double x3, double *xi, double *xj, double *xk, int c);

void get123(double *x1, double *x2, double *x3, double xi, double xj, double xk, int c);

void sumckdk_scdk_vtype(sumckdk_scdk *scdk, vtype vsodt, int fg, basis_list qnlist, matrix_t **mlsj, int nf, int nfp, int ni, int nip, pthread_mutex_t *mutex, int lock, int c);

void vcent(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

void vcont(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

void vtens(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

void vsorp(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

void vsorn(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c);

#endif