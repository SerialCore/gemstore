/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef MFI_H
#define MFI_H

#include <basis.h>
#include <matrix.h>
#include <sumckdk.h>
#include <inteCenV.h>
#include <thread.h>

typedef struct{
        basis_list qnlist_spfy;
        basis_list qnlist_full;
        matrix_t **mlsj;
        sumckdk_scdk ****scdk_cent_1;
        sumckdk_scdk ****scdk_cent_2;
        sumckdk_scdk ****scdk_cent_3;
        sumckdk_scdk ****scdk_cont_1;
        sumckdk_scdk ****scdk_cont_2;
        sumckdk_scdk ****scdk_cont_3;
        sumckdk_scdk ****scdk_tens_1;
        sumckdk_scdk ****scdk_tens_2;
        sumckdk_scdk ****scdk_tens_3;
        sumckdk_scdk ****scdk_sorp_1;
        sumckdk_scdk ****scdk_sorp_2;
        sumckdk_scdk ****scdk_sorp_3;
        sumckdk_scdk ****scdk_sorn_1;
        sumckdk_scdk ****scdk_sorn_2;
        sumckdk_scdk ****scdk_sorn_3;
	vargs varg;
    matrix_t Nfi;
	matrix_t Hfi;
	matrix_t v;
	double *e1;
	double *e2;
} margs;

void malloc_scdk(int *len_part, int len_list, sumckdk_scdk *****scdk);

void free_scdk(int *len_part, int len_list, sumckdk_scdk *****scdk);

void *calc_scdk(void *args);

void *calc_scdk_mt(void *args);

void basis_mlsj_sl(margs *marg);

void getlsj_sl(margs *marg, double m1, double m2, double m3, double rmin, double rmax, int nmax, int f12, double J, int P, int Lmax);

void *getmfi(void *args);

#endif