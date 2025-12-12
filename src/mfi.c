/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <mfi.h>
#include <spin.h>
#include <basis.h>
#include <vtype.h>
#include <matrix.h>
#include <sumckdk.h>
#include <inteCenV.h>
#include <parallel.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

void malloc_scdk(int *len_part, int len_list, sumckdk_scdk *****scdk)
{
    int i, j, k;
    *scdk = (sumckdk_scdk ****)malloc(sizeof(sumckdk_scdk ***) * len_list);
    for (i = 0; i < len_list; i++)
    {
        (*scdk)[i] = (sumckdk_scdk ***)malloc(sizeof(sumckdk_scdk **) * len_part[i]);
        for (j = 0; j < len_part[i]; j++)
        {
            (*scdk)[i][j] = (sumckdk_scdk **)malloc(sizeof(sumckdk_scdk *) * len_list);
            for (k = 0; k < len_list; k++)
            {
                (*scdk)[i][j][k] = (sumckdk_scdk *)malloc(sizeof(sumckdk_scdk) * len_part[k]);
            }
        }
    }
}

void free_scdk(int *len_part, int len_list, sumckdk_scdk *****scdk)
{
    int i, j, k, l;
    for (i = 0; i < len_list; i++)
    {
        for (j = 0; j < len_part[i]; j++)
        {
            for (k = 0; k < len_list; k++)
            {
                for (l = 0; l < len_part[k]; l++)
                {
                    sumckdk_scdk_free(&((*scdk)[i][j][k][l]));
                }
                free((*scdk)[i][j][k]);
            }
            free((*scdk)[i][j]);
        }
        free((*scdk)[i]);
    }
    free(*scdk);
}

void *calc_scdk(void *args)
{
    mt_args *mtargs = (mt_args *)args;
    int ith = mtargs->ith;
    pthread_mutex_t *mutex = mtargs->mutex;
    int lock = mtargs->lock;
    margs *arg = (margs *)mtargs->p;
    int len_list = arg->qnlist_spfy.len_list;
    int *len_part = arg->qnlist_spfy.len_part;
    int nf = ith;
    int nfp, ni, nip;

    for (nfp = 0; nfp < len_part[nf]; nfp++)
    {
        for (ni = 0; ni < len_list; ni++)
        {
            for (nip = 0; nip < len_part[ni]; nip++)
            {
                sumckdk_scdk_vtype(&(arg->scdk_cent_1[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                sumckdk_scdk_vtype(&(arg->scdk_cent_2[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                sumckdk_scdk_vtype(&(arg->scdk_cent_3[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);

                sumckdk_scdk_vtype(&(arg->scdk_cont_1[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                sumckdk_scdk_vtype(&(arg->scdk_cont_2[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                sumckdk_scdk_vtype(&(arg->scdk_cont_3[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);

                sumckdk_scdk_vtype(&(arg->scdk_tens_1[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                sumckdk_scdk_vtype(&(arg->scdk_tens_2[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                sumckdk_scdk_vtype(&(arg->scdk_tens_3[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);

                sumckdk_scdk_vtype(&(arg->scdk_sorp_1[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                sumckdk_scdk_vtype(&(arg->scdk_sorp_2[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                sumckdk_scdk_vtype(&(arg->scdk_sorp_3[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);

                sumckdk_scdk_vtype(&(arg->scdk_sorn_1[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                sumckdk_scdk_vtype(&(arg->scdk_sorn_2[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                sumckdk_scdk_vtype(&(arg->scdk_sorn_3[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
            }
        }
    }
    printf("ith=%d\n", ith);
    return NULL;
}

void *calc_scdk_mt(void *args)
{
    mt_args *mtargs = (mt_args *)args;
    int ith = mtargs->ith;
    pthread_mutex_t *mutex = mtargs->mutex;
    int lock = mtargs->lock;
    margs *arg = (margs *)mtargs->p;
    int len_list = arg->qnlist_spfy.len_list;
    int *len_part = arg->qnlist_spfy.len_part;
    int nf = ith / 15;
    int nfp, ni, nip;

    for (nfp = 0; nfp < len_part[nf]; nfp++)
    {
        for (ni = 0; ni < len_list; ni++)
        {
            for (nip = 0; nip < len_part[ni]; nip++)
            {
                switch (ith % 15)
                {
                case 0:
                    sumckdk_scdk_vtype(&(arg->scdk_cent_1[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                    break;
                case 1:
                    sumckdk_scdk_vtype(&(arg->scdk_cent_2[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                    break;
                case 2:
                    sumckdk_scdk_vtype(&(arg->scdk_cent_3[nf][nfp][ni][nip]), vcent, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
                    break;

                case 3:
                    sumckdk_scdk_vtype(&(arg->scdk_cont_1[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                    break;
                case 4:
                    sumckdk_scdk_vtype(&(arg->scdk_cont_2[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                    break;
                case 5:
                    sumckdk_scdk_vtype(&(arg->scdk_cont_3[nf][nfp][ni][nip]), vcont, 1, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
                    break;

                case 6:
                    sumckdk_scdk_vtype(&(arg->scdk_tens_1[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                    break;
                case 7:
                    sumckdk_scdk_vtype(&(arg->scdk_tens_2[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                    break;
                case 8:
                    sumckdk_scdk_vtype(&(arg->scdk_tens_3[nf][nfp][ni][nip]), vtens, 2, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
                    break;

                case 9:
                    sumckdk_scdk_vtype(&(arg->scdk_sorp_1[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                    break;
                case 10:
                    sumckdk_scdk_vtype(&(arg->scdk_sorp_2[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                    break;
                case 11:
                    sumckdk_scdk_vtype(&(arg->scdk_sorp_3[nf][nfp][ni][nip]), vsorp, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
                    break;

                case 12:
                    sumckdk_scdk_vtype(&(arg->scdk_sorn_1[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 1);
                    break;
                case 13:
                    sumckdk_scdk_vtype(&(arg->scdk_sorn_2[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 2);
                    break;
                case 14:
                    sumckdk_scdk_vtype(&(arg->scdk_sorn_3[nf][nfp][ni][nip]), vsorn, 3, arg->qnlist_spfy, arg->mlsj, nf, nfp, ni, nip, mutex, lock, 3);
                    break;

                default:
                    printf("error_calc_scdk_mt\n");
                    break;
                }
            }
        }
    }

    return NULL;
}

void basis_mlsj_sl(margs *marg)
{
    int i, j;
    int lrho, llam, L, c;
    double sij, S, J, MJ;
    double s1, s2, s3, si, sj, sk, ms1, ms2, ms3, msi, msj, msk, mrho, mlam;
    double coe, cgf;
    basis_list *qnlist = &(marg->qnlist_spfy);
    for (i = 0; i < qnlist->len_list; i++)
    {
        for (j = 0; j < qnlist->len_part[i]; j++)
        {
            matrix_init(&(marg->mlsj[i][j]), 0, 6);
            coe = qnlist->qnum[i][j].coe;
            s1 = qnlist->qnum[i][j].s1;
            s2 = qnlist->qnum[i][j].s2;
            s3 = qnlist->qnum[i][j].s3;
            c = qnlist->qnum[i][j].c;
            lrho = qnlist->qnum[i][j].lrho;
            llam = qnlist->qnum[i][j].llam;
            L = qnlist->qnum[i][j].L;
            sij = qnlist->qnum[i][j].sij;
            S = qnlist->qnum[i][j].jl;
            J = qnlist->qnum[i][j].J;
            MJ = J;

            getijk(s1, s2, s3, &si, &sj, &sk, c);
            for (msi = -si; msi <= si; msi++)
            {
                for (msj = -sj; msj <= sj; msj++)
                {
                    for (msk = -sk; msk <= sk; msk++)
                    {
                        for (mrho = -lrho; mrho <= lrho; mrho++)
                        {
                            for (mlam = -llam; mlam <= llam; mlam++)
                            {
                                cgf = coe * clebsch_gordan(si, msi, sj, msj, sij, msi + msj) * clebsch_gordan(sij, msi + msj, sk, msk, S, msi + msj + msk) * clebsch_gordan(lrho, mrho, llam, mlam, L, mrho + mlam) * clebsch_gordan(S, msi + msj + msk, L, mrho + mlam, J, MJ);
                                if (0 != cgf)
                                {
                                    get123(&ms1, &ms2, &ms3, msi, msj, msk, c);
                                    matrix_push(&(marg->mlsj[i][j]));
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][0] = cgf;
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][1] = ms1;
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][2] = ms2;
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][3] = ms3;
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][4] = mrho;
                                    marg->mlsj[i][j].p[marg->mlsj[i][j].n - 1][5] = mlam;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void getlsj_sl(margs *marg, double m1, double m2, double m3, double rmin, double rmax, int nmax, int f12, double J, int P, int Lmax)
{
    int i, lrho, llam, L;
    double si, sj, sk, sij, S, t1, t2, t3, tij, T;
    int c;

    basis_list_init(&(marg->qnlist_spfy));
    basis_list_init(&(marg->qnlist_full));

    si = 0.5;
    sj = 0.5;
    sk = 0.5;
    t1 = 0;
    t2 = 0;
    t3 = 0;
    tij = 0;
    T = 0;

    for (lrho = 0; lrho <= Lmax; lrho++)
    {
        for (llam = 0; llam <= Lmax; llam++)
        {
            for (L = abs(lrho - llam); L <= lrho + llam; L++)
            {
                for (sij = 0; sij <= 1; sij++)
                {
                    for (S = fabs(sij - sk); S <= fabs(sij + sk); S++)
                    {
                        if (lrho + llam <= Lmax && P == pow(-1, lrho + llam) && fabs(S - L) <= J && J <= S + L)
                        {
                            if (1 == f12 * pow(-1, 1 + sij + lrho))
                            {
                                c = 1;
                                basis_list_push(&(marg->qnlist_spfy), 1, -1, -1, 1.0, m1, m2, m3, si, sj, sk, t1, t2, t3, tij, T, c, lrho, llam, L, sij, S, J, 0, 0, 0, 0);
                            }
                            c = 2;
                            basis_list_push(&(marg->qnlist_spfy), 1, -1, -1, 1.0, m1, m2, m3, si, sj, sk, t1, t2, t3, tij, T, c, lrho, llam, L, sij, S, J, 0, 0, 0, 0);
                            c = 3;
                            basis_list_push(&(marg->qnlist_spfy), 0, -1, -1, f12 * pow(-1, 1 + sij + lrho), m1, m2, m3, si, sj, sk, t1, t2, t3, tij, T, c, lrho, llam, L, sij, S, J, 0, 0, 0, 0);
                        }
                    }
                }
            }
        }
    }

    basis_list_push_full(&(marg->qnlist_spfy), &(marg->qnlist_full), rmin, rmax, nmax);

    marg->mlsj = (matrix_t **)malloc(sizeof(matrix_t *) * (marg->qnlist_spfy.len_list));
    for (i = 0; i < marg->qnlist_spfy.len_list; i++)
    {
        marg->mlsj[i] = (matrix_t *)malloc(sizeof(matrix_t) * (marg->qnlist_spfy.len_part[i]));
    }
    basis_mlsj_sl(marg);
}

void *getmfi(void *args)
{
    mt_args *mtarg = (mt_args *)args;
    margs *arg = (margs *)mtarg->p;
    int nf = mtarg->ith;
    int nfp, ni, nip;
    int mapf1, mapf2, mapi1, mapi2;

    for (ni = 0; ni < arg->qnlist_full.len_list; ni++)
    {
        arg->Nfi.p[nf][ni] = 0;
        arg->Hfi.p[nf][ni] = 0;
        for (nfp = 0; nfp < arg->qnlist_full.len_part[nf]; nfp++)
        {
            for (nip = 0; nip < arg->qnlist_full.len_part[ni]; nip++)
            {
                mapf1 = arg->qnlist_full.qnum[nf][nfp].map1;
                mapf2 = arg->qnlist_full.qnum[nf][nfp].map2;
                mapi1 = arg->qnlist_full.qnum[ni][nip].map1;
                mapi2 = arg->qnlist_full.qnum[ni][nip].map2;

                arg->Nfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vNfi, 1);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vConLine, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vConLine, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vConLine, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCoul, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCoul, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCoul, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cont_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCont, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cont_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCont, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_cent, arg->scdk_cont_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vOgeCont, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_tens, arg->scdk_tens_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTens, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_tens, arg->scdk_tens_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTens, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_tens, arg->scdk_tens_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTens, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorp, arg->scdk_sorp_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorp, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorp, arg->scdk_sorp_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorp, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorp, arg->scdk_sorp_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorp, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorn, arg->scdk_sorn_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorn, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorn, arg->scdk_sorn_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorn, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(t1r_sorn, arg->scdk_sorn_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vSorn, 3);

                arg->Hfi.p[nf][ni] += inteVcenPart(tpi_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTi, 1);
                arg->Hfi.p[nf][ni] += inteVcenPart(tpi_cent, arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTi, 2);
                arg->Hfi.p[nf][ni] += inteVcenPart(tpi_cent, arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vTi, 3);

                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t1p_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT1, 1);
                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t2p_cent, arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT2, 1);

                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t1p_cent, arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT1, 2);
                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t2p_cent, arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT2, 2);

                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t1p_cent, arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT1, 3);
                arg->Hfi.p[nf][ni] += 0 * inteVcenPart(t2p_cent, arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2], arg->qnlist_full.qnum[nf][nfp], arg->qnlist_full.qnum[ni][nip], arg->varg, vT2, 3);
            }
        }
        arg->Hfi.p[nf][ni] += (arg->qnlist_full.qnum[nf][0].m1 + arg->qnlist_full.qnum[nf][0].m2 + arg->qnlist_full.qnum[nf][0].m3) * arg->Nfi.p[nf][ni];
    }
    return NULL;
}