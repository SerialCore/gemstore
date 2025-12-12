/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * 
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef EIGEN_H
#define EIGEN_H

#include <unistd.h>
#include <pthread.h>

typedef struct
{
    int l;

    int *cal_pth;

    pthread_mutex_t *mutex_subpro;
    pthread_mutex_t *mutex_parent;

    pthread_cond_t *cond_subpro;
    pthread_cond_t *cond_parent;

    int *shared_data_subpro;
    int *shared_data_parent;

    int *lth_use;
    int *lth_cal;

    int *type;
    int *stList;
    int *edList;
    double **a;
    int *i;
    double *uu;
    double *beta;
    int *jList;
    double *cList;
    double *sList;
    int *lx;
    double *vc;
    double *vs;
} trid_args;

typedef struct
{
    int l;

    int *cal_pth;

    pthread_mutex_t *mutex_subpro;
    pthread_mutex_t *mutex_parent;

    pthread_cond_t *cond_subpro;
    pthread_cond_t *cond_parent;

    int *shared_data_subpro;
    int *shared_data_parent;

    int *lth_use;
    int *lth_cal;

    int *type;
    int *stList;
    int *edList;
    double **a;
    double **b;
    double **vt;
    double **G;
    double **IG;
    double **IGA;
    double **S;
    int *i;
    int *ii;
    int *j;
    int n;
    int lt;
} eig2_args;

double drand();

void printArrayD1(double *a, int n);

void printArrayE1(double *a, int n);

void printArrayW1(double *a, int n);

void printArrayD2(double **a, int n, int m);

void sort(double *a, int n);

void trideigv(double **a, int n, double *d, double *e, double *dt, double *et, int lt);

void eigv1(double **a, int n, double *d, double *dt, double **vt, int lt);

void eigv2(double **a, double **b, int n, double *d, double *dt, double **vt, int lt, int *info);

void allocateThreads(int st, int ed, int lth, int min_lnum, int *stList, int *edList, int *lth_use);

void *trid_fun(void *args);

void trideigvMul(double **a, int n, double *d, double *e, double *dt, double *et, int lt);

void eigv1Mul(double **a, int n, double *d, double *dt, double **vt, int lt);

void *eig2_fun(void *args);

void eigv2Mul(double **a, double **b, int n, double *d, double *dt, double **vt, int lt, int *info);


#endif