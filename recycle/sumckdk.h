/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SUMCKDK_H
#define SUMCKDK_H

#include <unistd.h>
#include <pthread.h>

typedef struct {
	int len;
	double *c;
	double **d;
}sumckdk_cdlm;

typedef struct {
	int len;
	double *coe;
	int **nfgh;
}sumckdk_scdk;

void sumckdk_cdlm_init(sumckdk_cdlm *cdlm);

void sumckdk_cdlm_logs(sumckdk_cdlm cdlm);

void sumckdk_cdlm_push(sumckdk_cdlm *cdlm, double c, double dx, double dy, double dz);

void sumckdk_cdlm_cali(sumckdk_cdlm *cdlm, int l, int m);

void sumckdk_cdlm_calf(sumckdk_cdlm *cdlm, int l, int m);

void sumckdk_cdlm_free(sumckdk_cdlm *cdlm);

void sumckdk_scdk_init(sumckdk_scdk *scdk);

void sumckdk_scdk_free(sumckdk_scdk *scdk);

void sumckdk_scdk_logs(sumckdk_scdk scdk);

void sumckdk_scdk_push(sumckdk_scdk *scdk, double coe, int nf, int nf1, int nf2, int nf3, int nf4, int ng, int ng1, int ng2, int ng3, int ng4, int nh11, int nh12, int nh13, int nh14, int nh21, int nh22, int nh23, int nh24);

void sumckdk_scdk_spfy(sumckdk_scdk *scdk);

void sumckdk_scdk_calc(sumckdk_scdk *scdk, double cit, int l1, int l2, int l3, int l4, int lh1, int lh2, int m1, int m2, int m3, int m4, int mh1, int mh2, pthread_mutex_t *mutex, int lock);

double sumckdk_scdk_vals(sumckdk_scdk scdk, double coef, double f1, double f2, double f3, double f4, double coeg, double g1, double g2, double g3, double g4, double h11, double h12, double h13, double h14, double h21, double h22, double h23, double h24);

#endif