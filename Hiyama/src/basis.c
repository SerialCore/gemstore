/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <basis.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void basis_list_init(basis_list *qnlist)
{
    qnlist->qnum = (basis_qnum **)malloc(sizeof(basis_qnum *) * 0);
    qnlist->len_part = (int *)malloc(sizeof(int) * 0);
    qnlist->len_list = 0;
}

void basis_list_logs(basis_list qnlist)
{
    int i, j;
    printf("-------------------basis_list start--------------------\n");
    for (i = 0; i < qnlist.len_list; i++)
    {
        printf("%3d: |", i);
        for (j = 0; j < qnlist.len_part[i]; j++)
        {
            printf("j=%d map1=%3d map2=%3d coe=%10.6f c=%d nrho=%2d nlam=%2d |", j, qnlist.qnum[i][j].map1, qnlist.qnum[i][j].map2, qnlist.qnum[i][j].coe, qnlist.qnum[i][j].c, qnlist.qnum[i][j].nrho, qnlist.qnum[i][j].nlam);
        }
        printf("\n");
    }
    printf("-------------------basis_list end---------------------\n\n");
}

void basis_list_push(basis_list *qnlist, int newqnQ, int map1, int map2, double coe, double m1, double m2, double m3, double s1, double s2, double s3, double t1, double t2, double t3, double tij, double T, int c, int lrho, int llam, int L, double sij, double jl, double J, int nrho, int nlam, double nurho, double nulam)
{
    int i, j;
    if (1 == newqnQ)
    {
        i = qnlist->len_list;
        j = 0;
        qnlist->qnum = (basis_qnum **)realloc(qnlist->qnum, sizeof(basis_qnum *) * (i + 1));
        qnlist->qnum[i] = (basis_qnum *)malloc(sizeof(basis_qnum) * 1);
        qnlist->len_part = (int *)realloc(qnlist->len_part, sizeof(int) * (i + 1));
        qnlist->len_part[i] = 1;
        qnlist->len_list++;
    }
    else
    {
        i = qnlist->len_list - 1;
        j = qnlist->len_part[i];
        qnlist->qnum[i] = (basis_qnum *)realloc(qnlist->qnum[i], sizeof(basis_qnum) * (j + 1));
        qnlist->len_part[i]++;
    }

    if (map1 < 0 || map2 < 0)
    {
        qnlist->qnum[i][j].map1 = i;
        qnlist->qnum[i][j].map2 = j;
    }
    else
    {
        qnlist->qnum[i][j].map1 = map1;
        qnlist->qnum[i][j].map2 = map2;
    }

    qnlist->qnum[i][j].coe = coe;
    qnlist->qnum[i][j].m1 = m1;
    qnlist->qnum[i][j].m2 = m2;
    qnlist->qnum[i][j].m3 = m3;
    qnlist->qnum[i][j].s1 = s1;
    qnlist->qnum[i][j].s2 = s2;
    qnlist->qnum[i][j].s3 = s3;
    qnlist->qnum[i][j].t1 = t1;
    qnlist->qnum[i][j].t2 = t2;
    qnlist->qnum[i][j].t3 = t3;
    qnlist->qnum[i][j].tij = tij;
    qnlist->qnum[i][j].T = T;
    qnlist->qnum[i][j].c = c;
    qnlist->qnum[i][j].lrho = lrho;
    qnlist->qnum[i][j].llam = llam;
    qnlist->qnum[i][j].L = L;
    qnlist->qnum[i][j].sij = sij;
    qnlist->qnum[i][j].jl = jl;
    qnlist->qnum[i][j].J = J;
    qnlist->qnum[i][j].nrho = nrho;
    qnlist->qnum[i][j].nlam = nlam;
    qnlist->qnum[i][j].nurho = nurho;
    qnlist->qnum[i][j].nulam = nulam;
}

double getnu(double xmin, double xmax, int N, int n)
{
    double a = pow(xmax / xmin, 1. / (N - 1.));
    double xn = xmin * pow(a, n - 1);
    return 1. / xn / xn;
}

void basis_list_push_full(basis_list *qnlist_spfy, basis_list *qnlist_full, double rmin, double rmax, int nmax)
{
    int i, j, nrho, nlam, newqnQ;
    for (i = 0; i < qnlist_spfy->len_list; i++)
    {
        for (nrho = 1; nrho <= nmax; nrho++)
        {
            for (nlam = 1; nlam <= nmax; nlam++)
            {
                for (j = 0; j < qnlist_spfy->len_part[i]; j++)
                {
                    if (0 == j)
                    {
                        newqnQ = 1;
                    }
                    else
                    {
                        newqnQ = 0;
                    }

                    basis_list_push(qnlist_full, newqnQ,
                                    qnlist_spfy->qnum[i][j].map1,
                                    qnlist_spfy->qnum[i][j].map2,
                                    qnlist_spfy->qnum[i][j].coe,
                                    qnlist_spfy->qnum[i][j].m1,
                                    qnlist_spfy->qnum[i][j].m2,
                                    qnlist_spfy->qnum[i][j].m3,
                                    qnlist_spfy->qnum[i][j].s1,
                                    qnlist_spfy->qnum[i][j].s2,
                                    qnlist_spfy->qnum[i][j].s3,
                                    qnlist_spfy->qnum[i][j].t1,
                                    qnlist_spfy->qnum[i][j].t2,
                                    qnlist_spfy->qnum[i][j].t3,
                                    qnlist_spfy->qnum[i][j].tij,
                                    qnlist_spfy->qnum[i][j].T,
                                    qnlist_spfy->qnum[i][j].c,
                                    qnlist_spfy->qnum[i][j].lrho,
                                    qnlist_spfy->qnum[i][j].llam,
                                    qnlist_spfy->qnum[i][j].L,
                                    qnlist_spfy->qnum[i][j].sij,
                                    qnlist_spfy->qnum[i][j].jl,
                                    qnlist_spfy->qnum[i][j].J,
                                    nrho, nlam, getnu(rmin, rmax, nmax, nrho), getnu(rmin, rmax, nmax, nlam));
                }
            }
        }
    }
}