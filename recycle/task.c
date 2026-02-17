/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <task.h>
#include <mfi.h>
#include <eigen.h>
#include <sumckdk.h>
#include <inteCenV.h>
#include <thread.h>

#include <stdio.h>
#include <stdlib.h>

void getq(double mn, double ms, double mc, double mb, int q, double *mi);
void getq(double mn, double ms, double mc, double mb, int q, double *mi)
{
    switch (q) {
    case 1:
        *mi = mn;
        break;
    case 2:
        *mi = ms;
        break;
    case 3:
        *mi = mc;
        break;
    case 4:
        *mi = mb;
        break;
    default:
        printf("error_getq\n");
        break;
    }
}

void task_baryon(int q1, int q2, int q3, int f12, double J, int P, int Lmax)
{
    double GeV = 1;
    double MeV = 0.001 * GeV;
    double fm = 5.0676896 / GeV;

    double mq = 300 * MeV;
    double ms = 510 * MeV;
    double mc = 1750 * MeV;
    double mb = 5112 * MeV;

    double b = 0.165 * GeV * GeV;
    double K = 90 * MeV;
    double alpha_ss = 1.2;
    double alpha_so = 0.077;
    double alpha_ten = 0.077;
    double C = -1139 * MeV;
    double Lambda = 3.5 / fm;

    double rmin = 0.2 * fm;
    double rmax = 2.0 * fm;
    int nmax = 6;

    int i, j;

    double m1, m2, m3;
    char cP;

    margs marg;
    vargs varg;

    varg.mq = mq;
    varg.ms = ms;
    varg.mc = mc;
    varg.mb = mb;

    varg.b = b;
    varg.K = K;
    varg.alpha_ss = alpha_ss;
    varg.alpha_so = alpha_so;
    varg.alpha_ten = alpha_ten;
    varg.C = C;
    varg.Lambda = Lambda;

    marg.varg = varg;

    getq(mq, ms, mc, mb, q1, &m1);
    getq(mq, ms, mc, mb, q2, &m2);
    getq(mq, ms, mc, mb, q3, &m3);

    if (P > 0)
    {
        cP = '+';
    }
    else if (P < 0)
    {
        cP = '-';
    }
    else
    {
        cP = '?';
    }
    printf("%d/2^%c:", (int)(2 * J), cP);

    getlsj_sl(&marg, m1, m2, m3, rmin, rmax, nmax, f12, J, P, Lmax);

    /*
        basis_list_logs(marg.qnlist_spfy);
        basis_list_logs(marg.qnlist_full);
    */

    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_1));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_2));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_3));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_1));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_2));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_3));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_1));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_2));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_3));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_1));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_2));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_3));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_1));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_2));
    malloc_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_3));

    thread_load(15 * marg.qnlist_spfy.len_list, calc_scdk_mt, &marg, getNumCores());
    matrix_init(&(marg.Nfi), marg.qnlist_full.len_list, marg.qnlist_full.len_list);
    matrix_init(&(marg.Hfi), marg.qnlist_full.len_list, marg.qnlist_full.len_list);
    matrix_init(&(marg.v), marg.qnlist_full.len_list, marg.qnlist_full.len_list);
    marg.e1 = (double *)malloc(sizeof(double) * marg.qnlist_full.len_list);
    marg.e2 = (double *)malloc(sizeof(double) * marg.qnlist_full.len_list);
    thread_load(marg.qnlist_full.len_list, getmfi, &marg, getNumCores());

    /**/

    int n = marg.qnlist_full.len_list;
    double **Nfi = marg.Nfi.p;
    double **Hfi = marg.Hfi.p;
    double **v = marg.v.p;
    double *e1 = marg.e1;
    double *e2 = marg.e2;
    int info;

    /*printmatrix(marg->Nfi);
    printmatrix(marg->Hfi);*/
    eigv2Mul(Hfi, Nfi, n, e1, e2, v, n, &info);
    //array_print(marg->v);
    /*printf("info=%d\n",info);*/
    printArrayD1(e1, 3);

    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_1));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_2));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cent_3));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_1));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_2));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_cont_3));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_1));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_2));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_tens_3));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_1));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_2));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorp_3));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_1));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_2));
    free_scdk(marg.qnlist_spfy.len_part, marg.qnlist_spfy.len_list, &(marg.scdk_sorn_3));

    for (i = 0; i < marg.qnlist_spfy.len_list; i++)
    {
        for (j = 0; j < marg.qnlist_spfy.len_part[i]; j++)
        {
            matrix_free(&marg.mlsj[i][j]);
        }
        free(marg.mlsj[i]);
    }
    free(marg.mlsj);

    for (i = 0; i < marg.qnlist_spfy.len_list; i++)
    {
        free(marg.qnlist_spfy.qnum[i]);
    }
    free(marg.qnlist_spfy.qnum);
    free(marg.qnlist_spfy.len_part);

    for (i = 0; i < marg.qnlist_full.len_list; i++)
    {
        free(marg.qnlist_full.qnum[i]);
    }
    free(marg.qnlist_full.qnum);
    free(marg.qnlist_full.len_part);

    matrix_free(&marg.Nfi);
    matrix_free(&marg.Hfi);
    matrix_free(&marg.v);
    free(marg.e1);
    free(marg.e2);
}
