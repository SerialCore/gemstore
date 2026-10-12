/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/basis.h>
#include <gemstore/basis/orbit.h>

#include <stdio.h>
#include <stdlib.h>

void bases_3b_init(bases_3b *qnlist)
{
    qnlist->qnum = (basis_3b **)malloc(sizeof(basis_3b *) * 0);
    qnlist->len_part = (int *)malloc(sizeof(int) * 0);
    qnlist->len_list = 0;
}

void bases_3b_push(bases_3b *qnlist, int newqnQ, int map1, int map2,
    double coe, double m1, double m2, double m3,
    double s1, double s2, double s3,
    double t1, double t2, double t3, double tij, double T,
    int c, int lrho, int llam, int L, double sij, double jl, double J,
    int nrho, int nlam, double nurho, double nulam)
{
    int i;
    int j;

    if (newqnQ == 1) {
        i = qnlist->len_list;
        j = 0;
        qnlist->qnum = (basis_3b **)realloc(qnlist->qnum, sizeof(basis_3b *) * (i + 1));
        qnlist->qnum[i] = (basis_3b *)malloc(sizeof(basis_3b) * 1);
        qnlist->len_part = (int *)realloc(qnlist->len_part, sizeof(int) * (i + 1));
        qnlist->len_part[i] = 1;
        qnlist->len_list++;
    }
    else {
        i = qnlist->len_list - 1;
        j = qnlist->len_part[i];
        qnlist->qnum[i] = (basis_3b *)realloc(qnlist->qnum[i], sizeof(basis_3b) * (j + 1));
        qnlist->len_part[i]++;
    }

    if (map1 < 0 || map2 < 0) {
        qnlist->qnum[i][j].map1 = i;
        qnlist->qnum[i][j].map2 = j;
    }
    else {
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

void bases_3b_push_full(bases_3b *qnlist_spfy, bases_3b *qnlist_full,
    double rmin, double rmax, int nmax)
{
    int newqnQ;

    for (int i = 0; i < qnlist_spfy->len_list; i++) {
        for (int nrho = 1; nrho <= nmax; nrho++) {
            for (int nlam = 1; nlam <= nmax; nlam++) {
                for (int j = 0; j < qnlist_spfy->len_part[i]; j++) {
                    newqnQ = (j == 0) ? 1 : 0;
                    bases_3b_push(qnlist_full, newqnQ,
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
                        nrho,
                        nlam,
                        getnu(nrho, nmax, rmax, rmin),
                        getnu(nlam, nmax, rmax, rmin));
                }
            }
        }
    }
}

void bases_3b_free(bases_3b *qnlist)
{
    for (int i = 0; i < qnlist->len_list; i++) {
        free(qnlist->qnum[i]);
    }
    free(qnlist->qnum);
    free(qnlist->len_part);
    qnlist->qnum = NULL;
    qnlist->len_part = NULL;
    qnlist->len_list = 0;
}
