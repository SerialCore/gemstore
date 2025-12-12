/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <eigen.h>
#include <parallel.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

double drand()
{
    static unsigned char c0 = 100;
    static unsigned char c1 = 110;
    static unsigned char c2 = 180;
    static unsigned char c3 = 177;
    static unsigned char c4 = 202;
    static unsigned char c5 = 94;
    static unsigned char c6 = 75;
    double x;
    unsigned char *c;

    x = 1;
    c = (unsigned char *)&x;

    c0 = (c0 * c6 + c5) % 256;
    c1 = (c1 * c0 + c6) % 256;
    c2 = (c2 * c1 + c0) % 256;
    c3 = (c3 * c2 + c1) % 256;
    c4 = (c4 * c3 + c2) % 256;
    c5 = (c5 * c4 + c3) % 256;
    c6 = (c6 * 117 + 57) % 256;
    c6 = c6 % 16;

    if (0 == c[0]) {
        c[0] = c0;
        c[1] = c1;
        c[2] = c2;
        c[3] = c3;
        c[4] = c4;
        c[5] = c5;
        c[6] += c6;
    }
    else {
        c[6] = c0;
        c[5] = c1;
        c[4] = c2;
        c[3] = c3;
        c[2] = c4;
        c[1] = c5;
        c[0] += c6;
    }

    x = x - 1;
    return x;
}

void printArrayD1(double *a, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%15.10f  ", a[i]);
    }
    printf("\n");
}

void printArrayE1(double *a, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%15.6E  ", a[i]);
    }
    printf("\n");
}

void printArrayW1(double *a, int n)
{
    int i;
    printf("{");
    for (i = 0; i < n; i++) {
        printf("%15.10f", a[i]);
        if (i < n - 1) {
            printf(",");
        }
    }
    printf("}\n");
}

void printArrayD2(double **a, int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%10.6f  ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void sort(double *a, int n)
{
    int i, j;
    double temp;
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (a[j] > a[j + 1]) {
                temp = a[j];
                a[j] = a[j + 1];
                a[j + 1] = temp;
            }
        }
    }
}

void trideigv(double **a, int n, double *d, double *e, double *dt, double *et, int lt)
{
    int i, j, k;
    double sigma, beta, uu, eps, vmax;

    double r, c, s, t, ks;
    double d0, d1, e0, e1;
    double vc, vs;
    double *dd, *ee, *ds, *es, *as;
    int *od;
    double temp, g;

    eps = 1E-17;
    e[0] = 0;
    for (i = n - 1; i >= 1; i--) {
        sigma = 0;
        for (j = 0; j <= i - 1; j++) {
            sigma += a[i][j] * a[i][j];
        }
        sigma = sqrt(sigma);
        if (a[i][i - 1] < 0) {
            sigma *= -1;
        }
        d[i] = a[i][i];
        e[i] = -sigma;
        beta = sigma * (sigma + a[i][i - 1]);
        a[i][i - 1] += sigma;

        if (fabs(beta) > eps * fabs(a[i][i])) {
            for (j = 0; j < i; j++) {
                a[j][i] = 0;
                for (k = 0; k < i; k++) {
                    a[j][i] += a[j][k] * a[i][k];
                }
            }

            uu = 0;
            for (j = 0; j < i; j++) {
                uu += a[i][j] * a[j][i];
            }
            uu = uu / beta / beta;

            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    a[j][k] = a[j][k] - a[i][j] * a[k][i] / beta - a[i][k] * a[j][i] / beta + a[i][j] * a[i][k] * uu;
                }
            }
            a[i][i] = beta;
        }
        else
        {
            a[i][i] = 0;
        }
    }
    d[0] = a[0][0];

    if (lt > 0) {
        a[0][0] = 1;
        for (i = 1; i < n; i++) {
            if (0 != a[i][i]) {
                for (j = 0; j < i; j++) {
                    a[j][i] = 0;
                    for (k = 0; k < i; k++) {
                        a[j][i] += a[j][k] * a[i][k];
                    }
                    a[j][i] /= a[i][i];
                }
                for (j = 0; j < i; j++) {
                    for (k = 0; k < i; k++) {
                        a[j][k] -= a[j][i] * a[i][k];
                    }
                }
                for (j = 0; j < i; j++) {
                    a[i][j] = 0;
                    a[j][i] = 0;
                }
            }
            a[i][i] = 1;
        }
    }

    dd = (double *)malloc(sizeof(double) * n);
    ee = (double *)malloc(sizeof(double) * n);
    ds = (double *)malloc(sizeof(double) * n);
    es = (double *)malloc(sizeof(double) * n);
    as = (double *)malloc(sizeof(double) * n);
    od = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
        ds[i] = d[i];
        es[i] = e[i];
        od[i] = i;
    }

    t = 0;
    e[0] = 0;
    for (i = 0; i < n - 1; i++) {
        while (fabs(e[i + 1]) > eps * (fabs(d[i]) + fabs(d[i + 1]))) {
            g = (d[i + 1] - d[i]) / (2 * e[i + 1]);
            if (g >= 0) {
                ks = d[i] - e[i + 1] / (g + sqrt(1 + g * g));
            }
            else {
                ks = d[i] - e[i + 1] / (g - sqrt(1 + g * g));
            }
            for (j = n - 1; j > i; j--) {
                r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                if (0 == r) {
                    continue;
                }
                c = (d[j] - ks) / r;
                s = -e[j] / r;

                d0 = d[j - 1];
                d1 = d[j];

                e0 = e[j - 1];
                e1 = e[j];

                d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                e[j - 1] = c * e0;
                e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                if (j < n - 1) {
                    e[j + 1] = -s * t + c * e[j + 1];
                }

                t = -s * e0;

                if (lt >= n - 1) {
                    for (k = 0; k < n; k++) {
                        vc = a[j - 1][k];
                        vs = a[j][k];
                        a[j - 1][k] = c * vc + s * vs;
                        a[j][k] = -s * vc + c * vs;
                    }
                }
            }
        }
    }

    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (d[j] > d[j + 1]) {
                temp = d[j];
                d[j] = d[j + 1];
                d[j + 1] = temp;

                k = od[j];
                od[j] = od[j + 1];
                od[j + 1] = k;
            }
        }
    }

    for (i = 0; i < n; i++) {
        dd[i] = d[i];
        ee[i] = e[i];
    }

    if (lt >= n - 1) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                a[i][j] = as[od[i]];
            }
        }
    }
    else {
        for (i = 0; i < n; i++) {
            od[i] = -1;
            d[i] = ds[i];
            e[i] = es[i];
        }

        t = 0;
        e[0] = 0;
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                d[j] -= dd[i];
            }
            while (1) {
                int bk = 0;
                for (j = 0; j < n - 1; j++) {
                    if (fabs(e[j + 1]) <= eps * (fabs(d[j]) + fabs(d[j + 1])) && fabs(d[j]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[j] < 0) {
                        od[j] = i;
                        bk = 1;
                        break;
                    }
                }
                if (fabs(d[n - 1]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[n - 1] < 0 && 0 == bk) {
                    od[n - 1] = i;
                    bk = 1;
                }
                if (1 == bk) {
                    break;
                }

                ks = 0;
                for (j = n - 1; j > 0; j--) {
                    if (od[j - 1] < 0) {
                        r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                        if (0 == r) {
                            continue;
                        }
                        c = (d[j] - ks) / r;
                        s = -e[j] / r;

                        d0 = d[j - 1];
                        d1 = d[j];

                        e0 = e[j - 1];
                        e1 = e[j];

                        d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                        d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                        e[j - 1] = c * e0;
                        e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                        if (j < n - 1) {
                            e[j + 1] = -s * t + c * e[j + 1];
                        }

                        t = -s * e0;

                        for (k = 0; k < n; k++) {
                            vc = a[j - 1][k];
                            vs = a[j][k];
                            a[j - 1][k] = c * vc + s * vs;
                            a[j][k] = -s * vc + c * vs;
                        }
                    }
                }
            }
            for (j = 0; j < n; j++) {
                d[j] += dd[i];
            }
        }

        for (i = 0; i < n; i++) {
        }

        for (i = 0; i < n; i++) {
            as[i] = d[i];
        }

        for (i = 0; i < n; i++) {
            if (od[i] >= 0) {
                d[od[i]] = as[i];
            }
        }
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                if (od[i] >= 0) {
                    a[od[i]][j] = as[i];
                }
            }
        }
    }

    for (i = lt; i < n; i++) {
        for (j = 0; j < n; j++) {
            a[i][j] = 0;
        }
    }
    for (i = 0; i < lt; i++) {
        dt[i] = d[i];
        et[i] = e[i];
    }
    for (i = 0; i < n; i++) {
        d[i] = dd[i];
        e[i] = ee[i];
    }
    for (i = 0; i < lt; i++) {
        vmax = 0;
        for (j = 0; j < n; j++) {
            if (fabs(a[i][j]) > fabs(vmax)) {
                vmax = a[i][j];
            }
        }
        if (vmax < 0) {
            for (j = 0; j < n; j++) {
                a[i][j] *= -1;
            }
        }
    }

    free(dd);
    free(ee);
    free(ds);
    free(es);
    free(as);
    free(od);
}

void eigv1(double **a, int n, double *d, double *dt, double **vt, int lt)
{
    double **aa, *e, *et;
    int i, j;
    aa = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    for (i = 0; i < n; i++) {
        aa[i] = (double *)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
    }
    trideigv(aa, n, d, e, dt, et, lt);
    for (i = 0; i < lt; i++) {
        for (j = 0; j < n; j++) {
            vt[i][j] = aa[i][j];
        }
        free(aa[i]);
    }
    free(aa);
    free(e);
    free(et);
}

void eigv2(double **a, double **b, int n, double *d, double *dt, double **vt, int lt, int *info)
{

    double **G, **IG, **IGA, **S, *e, *et;
    int i, j, k, ii;
    double s, ds;

    G = (double **)malloc(sizeof(double *) * n);
    IG = (double **)malloc(sizeof(double *) * n);
    IGA = (double **)malloc(sizeof(double *) * n);
    S = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);

    for (i = 0; i < n; i++) {
        G[i] = (double *)malloc(sizeof(double) * n);
        IG[i] = (double *)malloc(sizeof(double) * n);
        IGA[i] = (double *)malloc(sizeof(double) * n);
        S[i] = (double *)malloc(sizeof(double) * n);
    }

    for (j = 0; j < n; j++) {
        s = 0;
        for (k = 0; k <= j - 1; k++) {
            s = s + G[j][k] * G[j][k];
        }
        ds = b[j][j] - s;
        if (ds <= 0) {
            printf("error_cholesky\n");
            *info = -1;
            return;
        }
        ds = fabs(ds);
        G[j][j] = sqrt(ds);

        for (i = j + 1; i < n; i++) {
            s = 0;
            for (k = 0; k <= j - 1; k++) {
                s = s + G[i][k] * G[j][k];
            }
            G[i][j] = (b[i][j] - s) / G[j][j];
        }
    }

    for (ii = 0; ii < n; ii++) {
        for (i = 0; i < n; i++) {
            IG[i][ii] = 0;
        }
        IG[ii][ii] = 1;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                IG[j][ii] -= G[j][i] / G[i][i] * IG[i][ii];
            }
        }
        for (i = ii; i < n; i++) {
            IG[i][ii] /= G[i][i];
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            IGA[i][j] = 0;
            for (k = 0; k <= i; k++) {
                IGA[i][j] += IG[i][k] * a[k][j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            S[i][j] = 0;
            for (k = 0; k <= j; k++) {
                S[i][j] += IGA[i][k] * IG[j][k];
            }
            S[j][i] = S[i][j];
        }
    }

    trideigv(S, n, d, e, dt, et, lt);

    for (j = 0; j < n; j++) {
        for (i = 0; i < lt; i++) {
            vt[i][j] = 0;
            for (k = j; k < n; k++) {
                vt[i][j] += S[i][k] * IG[k][j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        free(G[i]);
        free(IG[i]);
        free(IGA[i]);
        free(S[i]);
    }
    free(G);
    free(IG);
    free(IGA);
    free(S);
    free(e);
    free(et);
}

void allocateThreads(int st, int ed, int lth, int min_lnum, int *stList, int *edList, int *lth_use)
{
    int numTot = ed - st + 1;
    int i, nuse;
    if (lth * min_lnum <= numTot)
    {
        nuse = numTot / lth;
        if (numTot % lth != 0)
        {
            nuse++;
        }
    }
    else
    {
        nuse = min_lnum;
    }
    stList[0] = st;
    for (i = 0; i < lth; i++)
    {
        edList[i] = stList[i] + nuse - 1;
        if (edList[i] >= ed)
        {
            edList[i] = ed;
            *lth_use = i + 1;
            break;
        }
        if (i < lth - 1)
        {
            stList[i + 1] = edList[i] + 1;
        }
    }
    for (i = *lth_use; i < lth; i++)
    {
        stList[i] = 0;
        edList[i] = -1;
    }
}

void *trid_fun(void *args)
{
    trid_args *arg = (trid_args *)args;
    int l = arg->l;

    int *cal_pth = arg->cal_pth;
    pthread_mutex_t *mutex_subpro = arg->mutex_subpro;
    pthread_mutex_t *mutex_parent = arg->mutex_parent;

    pthread_cond_t *cond_subpro = arg->cond_subpro;
    pthread_cond_t *cond_parent = arg->cond_parent;

    int *shared_data_subpro = arg->shared_data_subpro;
    int *shared_data_parent = arg->shared_data_parent;

    int *lth_use = arg->lth_use;
    int *lth_cal = arg->lth_cal;

    int *type = arg->type;
    int *stList = arg->stList;
    int *edList = arg->edList;
    double **a = arg->a;
    int *ptr_i = arg->i;
    double *ptr_uu = arg->uu;
    double *ptr_beta = arg->beta;
    int *jList = arg->jList;
    double *cList = arg->cList;
    double *sList = arg->sList;
    int *ptr_lx = arg->lx;

    int i, j, k, lx;
    double uu, beta;
    double c, s, vc, vs;

    while (*cal_pth)
    {
        pthread_mutex_lock(&mutex_subpro[l]);

        while (0 == shared_data_subpro[l])
        {
            pthread_mutex_lock(mutex_parent);
            lth_cal[0]++;
            if (lth_cal[0] == lth_use[0])
            {
                shared_data_parent[0] = 1;
                pthread_cond_signal(cond_parent);
            }
            pthread_mutex_unlock(mutex_parent);
            pthread_cond_wait(&cond_subpro[l], &mutex_subpro[l]);
        }
        shared_data_subpro[l] = 0;

        i = *ptr_i;
        uu = *ptr_uu;
        beta = *ptr_beta;
        lx = *ptr_lx;

        switch (*type)
        {
        case 1:
            for (j = stList[l]; j <= edList[l]; j++)
            {
                a[j][i] = 0;
                for (k = 0; k < i; k++)
                {
                    a[j][i] += a[j][k] * a[i][k];
                }
            }
            break;

        case 2:
            for (j = stList[l]; j <= edList[l]; j++)
            {
                for (k = 0; k < i; k++)
                {
                    a[j][k] = a[j][k] - a[i][j] * a[k][i] / beta - a[i][k] * a[j][i] / beta + a[i][j] * a[i][k] * uu;
                }
            }
            break;

        case 3:
            for (j = stList[l]; j <= edList[l]; j++)
            {
                a[j][i] = 0;
                for (k = 0; k < i; k++)
                {
                    a[j][i] += a[j][k] * a[i][k];
                }
                a[j][i] /= a[i][i];
            }
            break;

        case 4:
            for (j = stList[l]; j <= edList[l]; j++)
            {
                for (k = 0; k < i; k++)
                {
                    a[j][k] -= a[j][i] * a[i][k];
                }
            }
            break;

        case 5:
            for (i = 0; i <= lx; i++)
            {
                j = jList[i];
                c = cList[i];
                s = sList[i];
                for (k = stList[l]; k <= edList[l]; k++)
                {
                    vc = a[j - 1][k];
                    vs = a[j][k];
                    a[j - 1][k] = c * vc + s * vs;
                    a[j][k] = -s * vc + c * vs;
                }
            }
            break;

        case -1:
            break;

        case 0:
            break;

        default:
            break;
        }
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    return NULL;
}

void trideigvMul(double **a, int n, double *d, double *e, double *dt, double *et, int lt)
{
    int l;
    int cal_pth;
    pthread_mutex_t *mutex_subpro;
    pthread_mutex_t mutex_parent;

    pthread_cond_t *cond_subpro;
    pthread_cond_t cond_parent;

    int *shared_data_subpro;
    int shared_data_parent;

    int lth, min_lnum, lth_use, lth_cal;

    int type, *stList, *edList, *jList;
    double *cList, *sList;
    int lx_max, lx;

    trid_args *arg;
    pthread_t *tid;

    int i, j, k;
    double sigma, beta, uu, eps, vmax;

    double r, c, s, t, ks;
    double d0, d1, e0, e1;
    double vc, vs;
    double *dd, *ee, *ds, *es, *as;
    int *od;
    double temp, g;

    min_lnum = 1;

    lth = getNumProcessors();

    mutex_subpro = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t) * lth);
    cond_subpro = (pthread_cond_t *)malloc(sizeof(pthread_cond_t) * lth);
    shared_data_subpro = (int *)malloc(sizeof(int) * lth);
    stList = (int *)malloc(sizeof(int) * lth);
    edList = (int *)malloc(sizeof(int) * lth);
    tid = (pthread_t *)malloc(sizeof(pthread_t) * lth);
    arg = (trid_args *)malloc(sizeof(trid_args) * lth);

    lx_max = 16 * 16 * 16 * 16;
    jList = (int *)malloc(sizeof(int) * lx_max);
    cList = (double *)malloc(sizeof(double) * lx_max);
    sList = (double *)malloc(sizeof(double) * lx_max);

    for (l = 0; l < lth; l++)
    {
        arg[l].l = l;
        arg[l].cal_pth = &cal_pth;
        arg[l].mutex_subpro = mutex_subpro;
        arg[l].mutex_parent = &mutex_parent;
        arg[l].cond_subpro = cond_subpro;
        arg[l].cond_parent = &cond_parent;
        arg[l].shared_data_subpro = shared_data_subpro;
        arg[l].shared_data_parent = &shared_data_parent;
        arg[l].lth_use = &lth_use;
        arg[l].lth_cal = &lth_cal;

        arg[l].type = &type;
        arg[l].stList = stList;
        arg[l].edList = edList;
        arg[l].a = a;
        arg[l].i = &i;
        arg[l].uu = &uu;
        arg[l].beta = &beta;
        arg[l].jList = jList;
        arg[l].cList = cList;
        arg[l].sList = sList;
        arg[l].lx = &lx;
        arg[l].vc = &vc;
        arg[l].vs = &vs;

        pthread_mutex_init(&mutex_subpro[l], NULL);
        pthread_cond_init(&cond_subpro[l], NULL);
        shared_data_subpro[l] = 0;
    }

    cal_pth = 1;
    shared_data_parent = 0;
    pthread_mutex_init(&mutex_parent, NULL);
    pthread_cond_init(&cond_parent, NULL);

    type = 0;
    lth_cal = 0;
    lth_use = lth;
    for (l = 0; l < lth; l++)
    {
        pthread_create(&tid[l], NULL, trid_fun, &arg[l]);
    }

    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    type = -1;
    lth_cal = 0;
    lth_use = lth;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }

    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    eps = 1E-17;
    e[0] = 0;
    for (i = n - 1; i >= 1; i--)
    {
        sigma = 0;
        for (j = 0; j <= i - 1; j++)
        {
            sigma += a[i][j] * a[i][j];
        }
        sigma = sqrt(sigma);
        if (a[i][i - 1] < 0)
        {
            sigma *= -1;
        }
        d[i] = a[i][i];
        e[i] = -sigma;
        beta = sigma * (sigma + a[i][i - 1]);
        a[i][i - 1] += sigma;

        if (fabs(beta) > eps * fabs(a[i][i]))
        {
            allocateThreads(0, i - 1, lth, min_lnum, stList, edList, &lth_use);

            type = 1;
            lth_cal = 0;
            for (l = 0; l < lth_use; l++)
            {
                pthread_mutex_lock(&mutex_subpro[l]);
                shared_data_subpro[l] = 1;
                pthread_cond_signal(&cond_subpro[l]);
                pthread_mutex_unlock(&mutex_subpro[l]);
            }
            pthread_mutex_lock(&mutex_parent);
            while (0 == shared_data_parent)
            {
                pthread_cond_wait(&cond_parent, &mutex_parent);
            }
            shared_data_parent = 0;
            pthread_mutex_unlock(&mutex_parent);

            /*
                        for(j=0;j<i;j++)
                        {
                            a[j][i]=0;
                            for(k=0;k<i;k++)
                            {
                                a[j][i]+=a[j][k]*a[i][k];
                            }
                        }
            */
            uu = 0;
            for (j = 0; j < i; j++)
            {
                uu += a[i][j] * a[j][i];
            }
            uu = uu / beta / beta;

            type = 2;
            lth_cal = 0;
            for (l = 0; l < lth_use; l++)
            {
                pthread_mutex_lock(&mutex_subpro[l]);
                shared_data_subpro[l] = 1;
                pthread_cond_signal(&cond_subpro[l]);
                pthread_mutex_unlock(&mutex_subpro[l]);
            }
            pthread_mutex_lock(&mutex_parent);
            while (0 == shared_data_parent)
            {
                pthread_cond_wait(&cond_parent, &mutex_parent);
            }
            shared_data_parent = 0;
            pthread_mutex_unlock(&mutex_parent);

            /*
                        for(j=0;j<i;j++)
                        {
                            for(k=0;k<i;k++)
                            {
                                a[j][k]=a[j][k]-a[i][j]*a[k][i]/beta-a[i][k]*a[j][i]/beta+a[i][j]*a[i][k]*uu;
                            }
                        }
            */
            a[i][i] = beta;
        }
        else
        {
            a[i][i] = 0;
        }
    }
    d[0] = a[0][0];

    if (lt > 0)
    {
        a[0][0] = 1;
        for (i = 1; i < n; i++)
        {
            if (0 != a[i][i])
            {

                allocateThreads(0, i - 1, lth, min_lnum, stList, edList, &lth_use);

                type = 3;
                lth_cal = 0;
                for (l = 0; l < lth_use; l++)
                {
                    pthread_mutex_lock(&mutex_subpro[l]);
                    shared_data_subpro[l] = 1;
                    pthread_cond_signal(&cond_subpro[l]);
                    pthread_mutex_unlock(&mutex_subpro[l]);
                }
                pthread_mutex_lock(&mutex_parent);
                while (0 == shared_data_parent)
                {
                    pthread_cond_wait(&cond_parent, &mutex_parent);
                }
                shared_data_parent = 0;
                pthread_mutex_unlock(&mutex_parent);

                type = 4;
                lth_cal = 0;
                for (l = 0; l < lth_use; l++)
                {
                    pthread_mutex_lock(&mutex_subpro[l]);
                    shared_data_subpro[l] = 1;
                    pthread_cond_signal(&cond_subpro[l]);
                    pthread_mutex_unlock(&mutex_subpro[l]);
                }
                pthread_mutex_lock(&mutex_parent);
                while (0 == shared_data_parent)
                {
                    pthread_cond_wait(&cond_parent, &mutex_parent);
                }
                shared_data_parent = 0;
                pthread_mutex_unlock(&mutex_parent);

                /*				for(j=0;j<i;j++)
                                {
                                    a[j][i]=0;
                                    for(k=0;k<i;k++)
                                    {
                                        a[j][i]+=a[j][k]*a[i][k];
                                    }
                                    a[j][i]/=a[i][i];
                                }
                                for(j=0;j<i;j++)
                                {
                                    for(k=0;k<i;k++)
                                    {
                                        a[j][k]-=a[j][i]*a[i][k];
                                    }
                                }
                */
                for (j = 0; j < i; j++)
                {
                    a[i][j] = 0;
                    a[j][i] = 0;
                }
            }
            a[i][i] = 1;
        }
    }

    dd = (double *)malloc(sizeof(double) * n);
    ee = (double *)malloc(sizeof(double) * n);
    ds = (double *)malloc(sizeof(double) * n);
    es = (double *)malloc(sizeof(double) * n);
    as = (double *)malloc(sizeof(double) * n);
    od = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
    {
        ds[i] = d[i];
        es[i] = e[i];
        od[i] = i;
    }

    allocateThreads(0, n - 1, lth, min_lnum, stList, edList, &lth_use);
    lx = -1;

    t = 0;
    e[0] = 0;
    for (i = 0; i < n - 1; i++)
    {
        while (fabs(e[i + 1]) > eps * (fabs(d[i]) + fabs(d[i + 1])))
        {
            g = (d[i + 1] - d[i]) / (2 * e[i + 1]);
            if (g >= 0)
            {
                ks = d[i] - e[i + 1] / (g + sqrt(1 + g * g));
            }
            else
            {
                ks = d[i] - e[i + 1] / (g - sqrt(1 + g * g));
            }
            for (j = n - 1; j > i; j--)
            {
                r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                if (0 == r)
                {
                    continue;
                }
                c = (d[j] - ks) / r;
                s = -e[j] / r;

                d0 = d[j - 1];
                d1 = d[j];

                e0 = e[j - 1];
                e1 = e[j];

                d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                e[j - 1] = c * e0;
                e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                if (j < n - 1)
                {
                    e[j + 1] = -s * t + c * e[j + 1];
                }

                t = -s * e0;

                if (lt >= n - 1)
                {
                    lx++;
                    jList[lx] = j;
                    cList[lx] = c;
                    sList[lx] = s;
                    if (lx == lx_max - 1)
                    {
                        type = 5;
                        lth_cal = 0;
                        for (l = 0; l < lth_use; l++)
                        {
                            pthread_mutex_lock(&mutex_subpro[l]);
                            shared_data_subpro[l] = 1;
                            pthread_cond_signal(&cond_subpro[l]);
                            pthread_mutex_unlock(&mutex_subpro[l]);
                        }
                        pthread_mutex_lock(&mutex_parent);
                        while (0 == shared_data_parent)
                        {
                            pthread_cond_wait(&cond_parent, &mutex_parent);
                        }
                        shared_data_parent = 0;
                        pthread_mutex_unlock(&mutex_parent);

                        lx = -1;
                    }
                    /*
                                        for(k=0;k<n;k++)
                                        {
                                            vc=a[j-1][k];
                                            vs=a[j  ][k];
                                            a[j-1][k]= c*vc+s*vs;
                                            a[j  ][k]=-s*vc+c*vs;
                                        }
                    */
                }
            }
        }
    }

    if (lt >= n - 1)
    {
        type = 5;
        lth_cal = 0;
        for (l = 0; l < lth_use; l++)
        {
            pthread_mutex_lock(&mutex_subpro[l]);
            shared_data_subpro[l] = 1;
            pthread_cond_signal(&cond_subpro[l]);
            pthread_mutex_unlock(&mutex_subpro[l]);
        }
        pthread_mutex_lock(&mutex_parent);
        while (0 == shared_data_parent)
        {
            pthread_cond_wait(&cond_parent, &mutex_parent);
        }
        shared_data_parent = 0;
        pthread_mutex_unlock(&mutex_parent);
    }

    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (d[j] > d[j + 1])
            {
                temp = d[j];
                d[j] = d[j + 1];
                d[j + 1] = temp;

                k = od[j];
                od[j] = od[j + 1];
                od[j + 1] = k;
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        dd[i] = d[i];
        ee[i] = e[i];
    }

    if (lt >= n - 1)
    {

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++)
            {
                a[i][j] = as[od[i]];
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            od[i] = -1;
            d[i] = ds[i];
            e[i] = es[i];
        }

        lx = -1;
        t = 0;
        e[0] = 0;
        for (i = 0; i < lt; i++)
        {
            for (j = 0; j < n; j++)
            {
                d[j] -= dd[i];
            }
            while (1)
            {
                int bk = 0;
                for (j = 0; j < n - 1; j++)
                {
                    if (fabs(e[j + 1]) <= eps * (fabs(d[j]) + fabs(d[j + 1])) && fabs(d[j]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[j] < 0)
                    {
                        od[j] = i;
                        bk = 1;
                        break;
                    }
                }
                if (fabs(d[n - 1]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[n - 1] < 0 && 0 == bk)
                {
                    od[n - 1] = i;
                    bk = 1;
                }
                if (1 == bk)
                {
                    break;
                }

                ks = 0;
                for (j = n - 1; j > 0; j--)
                {
                    if (od[j - 1] < 0)
                    {
                        r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                        if (0 == r)
                        {
                            continue;
                        }
                        c = (d[j] - ks) / r;
                        s = -e[j] / r;

                        d0 = d[j - 1];
                        d1 = d[j];

                        e0 = e[j - 1];
                        e1 = e[j];

                        d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                        d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;

                        e[j - 1] = c * e0;
                        e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                        if (j < n - 1)
                        {
                            e[j + 1] = -s * t + c * e[j + 1];
                        }

                        t = -s * e0;

                        lx++;
                        jList[lx] = j;
                        cList[lx] = c;
                        sList[lx] = s;
                        if (lx == lx_max - 1)
                        {
                            type = 5;
                            lth_cal = 0;
                            for (l = 0; l < lth_use; l++)
                            {
                                pthread_mutex_lock(&mutex_subpro[l]);
                                shared_data_subpro[l] = 1;
                                pthread_cond_signal(&cond_subpro[l]);
                                pthread_mutex_unlock(&mutex_subpro[l]);
                            }
                            pthread_mutex_lock(&mutex_parent);
                            while (0 == shared_data_parent)
                            {
                                pthread_cond_wait(&cond_parent, &mutex_parent);
                            }
                            shared_data_parent = 0;
                            pthread_mutex_unlock(&mutex_parent);

                            lx = -1;
                        }
                        /*
                                            for(k=0;k<n;k++)
                                            {
                                                vc=a[j-1][k];
                                                vs=a[j  ][k];
                                                a[j-1][k]= c*vc+s*vs;
                                                a[j  ][k]=-s*vc+c*vs;
                                            }
                        */
                    }
                }
            }
            for (j = 0; j < n; j++)
            {
                d[j] += dd[i];
            }
        }

        type = 5;
        lth_cal = 0;
        for (l = 0; l < lth_use; l++)
        {
            pthread_mutex_lock(&mutex_subpro[l]);
            shared_data_subpro[l] = 1;
            pthread_cond_signal(&cond_subpro[l]);
            pthread_mutex_unlock(&mutex_subpro[l]);
        }
        pthread_mutex_lock(&mutex_parent);
        while (0 == shared_data_parent)
        {
            pthread_cond_wait(&cond_parent, &mutex_parent);
        }
        shared_data_parent = 0;
        pthread_mutex_unlock(&mutex_parent);

        for (i = 0; i < n; i++)
        {
        }

        for (i = 0; i < n; i++)
        {
            as[i] = d[i];
        }

        for (i = 0; i < n; i++)
        {
            if (od[i] >= 0)
            {
                d[od[i]] = as[i];
            }
        }
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++)
            {
                if (od[i] >= 0)
                {
                    a[od[i]][j] = as[i];
                }
            }
        }
    }

    for (i = lt; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            a[i][j] = 0;
        }
    }
    for (i = 0; i < lt; i++)
    {
        dt[i] = d[i];
        et[i] = e[i];
    }
    for (i = 0; i < n; i++)
    {
        d[i] = dd[i];
        e[i] = ee[i];
    }
    for (i = 0; i < lt; i++)
    {
        vmax = 0;
        for (j = 0; j < n; j++)
        {
            if (fabs(a[i][j]) > fabs(vmax))
            {
                vmax = a[i][j];
            }
        }
        if (vmax < 0)
        {
            for (j = 0; j < n; j++)
            {
                a[i][j] *= -1;
            }
        }
    }

    cal_pth = 0;
    type = 0;
    lth_cal = 0;
    lth_use = lth;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }

    for (l = 0; l < lth; l++)
    {
        pthread_join(tid[l], NULL);
    }

    free(dd);
    free(ee);
    free(ds);
    free(es);
    free(as);
    free(od);

    for (l = 0; l < lth; l++)
    {
        pthread_mutex_destroy(&mutex_subpro[l]);
        pthread_cond_destroy(&cond_subpro[l]);
    }
    free(mutex_subpro);
    free(cond_subpro);
    free(shared_data_subpro);
    free(stList);
    free(edList);
    free(jList);
    free(cList);
    free(sList);
    free(arg);
    free(tid);
}

void eigv1Mul(double **a, int n, double *d, double *dt, double **vt, int lt)
{
    double **aa, *e, *et;
    int i, j;
    aa = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    for (i = 0; i < n; i++)
    {
        aa[i] = (double *)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++)
        {
            aa[i][j] = a[i][j];
        }
    }
    trideigvMul(aa, n, d, e, dt, et, lt);
    for (i = 0; i < lt; i++)
    {
        for (j = 0; j < n; j++)
        {
            vt[i][j] = aa[i][j];
        }
        free(aa[i]);
    }
    free(aa);
    free(e);
    free(et);
}

void *eig2_fun(void *args)
{
    eig2_args *arg = (eig2_args *)args;
    int l = arg->l;

    int *cal_pth = arg->cal_pth;
    pthread_mutex_t *mutex_subpro = arg->mutex_subpro;
    pthread_mutex_t *mutex_parent = arg->mutex_parent;

    pthread_cond_t *cond_subpro = arg->cond_subpro;
    pthread_cond_t *cond_parent = arg->cond_parent;

    int *shared_data_subpro = arg->shared_data_subpro;
    int *shared_data_parent = arg->shared_data_parent;

    int *lth_use = arg->lth_use;
    int *lth_cal = arg->lth_cal;

    int *type = arg->type;
    int *stList = arg->stList;
    int *edList = arg->edList;
    double **a = arg->a;
    double **b = arg->b;
    double **vt = arg->vt;
    double **G = arg->G;
    double **IG = arg->IG;
    double **IGA = arg->IGA;
    double **S = arg->S;
    int *ptr_j = arg->j;
    int n = arg->n;
    int lt = arg->lt;

    int i, j, k, ii;
    double s;

    while (*cal_pth)
    {
        pthread_mutex_lock(&mutex_subpro[l]);

        while (0 == shared_data_subpro[l])
        {
            pthread_mutex_lock(mutex_parent);
            lth_cal[0]++;
            if (lth_cal[0] == lth_use[0])
            {
                shared_data_parent[0] = 1;
                pthread_cond_signal(cond_parent);
            }
            pthread_mutex_unlock(mutex_parent);
            pthread_cond_wait(&cond_subpro[l], &mutex_subpro[l]);
        }
        shared_data_subpro[l] = 0;

        switch (*type)
        {
        case 1:
            j = *ptr_j;
            for (i = stList[l]; i <= edList[l]; i++)
            {
                s = 0;
                for (k = 0; k <= j - 1; k++)
                {
                    s = s + G[i][k] * G[j][k];
                }
                G[i][j] = (b[i][j] - s) / G[j][j];
            }
            break;

        case 2:
            for (ii = stList[l]; ii <= edList[l]; ii++)
            {
                for (i = 0; i < n; i++)
                {
                    IG[i][ii] = 0;
                }
                IG[ii][ii] = 1;
                for (i = 0; i < n; i++)
                {
                    for (j = i + 1; j < n; j++)
                    {
                        IG[j][ii] -= G[j][i] / G[i][i] * IG[i][ii];
                    }
                }
                for (i = ii; i < n; i++)
                {
                    IG[i][ii] /= G[i][i];
                }
            }
            break;

        case 3:
            for (i = stList[l]; i <= edList[l]; i++)
            {
                for (j = 0; j < n; j++)
                {
                    IGA[i][j] = 0;
                    for (k = 0; k <= i; k++)
                    {
                        IGA[i][j] += IG[i][k] * a[k][j];
                    }
                }
            }
            break;

        case 4:
            for (i = stList[l]; i <= edList[l]; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    S[i][j] = 0;
                    for (k = 0; k <= j; k++)
                    {
                        S[i][j] += IGA[i][k] * IG[j][k];
                    }
                    S[j][i] = S[i][j];
                }
            }
            break;

        case 5:
            for (j = stList[l]; j <= edList[l]; j++)
            {
                for (i = 0; i < lt; i++)
                {
                    vt[i][j] = 0;
                    for (k = j; k < n; k++)
                    {
                        vt[i][j] += S[i][k] * IG[k][j];
                    }
                }
            }

        case -1:
            break;

        case 0:
            break;

        default:
            break;
        }
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    return NULL;
}

void eigv2Mul(double **a, double **b, int n, double *d, double *dt, double **vt, int lt, int *info)
{

    int l;
    int cal_pth;
    pthread_mutex_t *mutex_subpro;
    pthread_mutex_t mutex_parent;

    pthread_cond_t *cond_subpro;
    pthread_cond_t cond_parent;

    int *shared_data_subpro;
    int shared_data_parent;

    int lth, min_lnum, lth_use, lth_cal;
    int type, *stList, *edList;

    eig2_args *arg;
    pthread_t *tid;

    double **G, **IG, **IGA, **S, *e, *et;
    int i, j, k, ii;
    double s, ds;

    G = (double **)malloc(sizeof(double *) * n);
    IG = (double **)malloc(sizeof(double *) * n);
    IGA = (double **)malloc(sizeof(double *) * n);
    S = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);

    for (i = 0; i < n; i++)
    {
        G[i] = (double *)malloc(sizeof(double) * n);
        IG[i] = (double *)malloc(sizeof(double) * n);
        IGA[i] = (double *)malloc(sizeof(double) * n);
        S[i] = (double *)malloc(sizeof(double) * n);
    }

    min_lnum = 1;

    lth = getNumProcessors();

    mutex_subpro = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t) * lth);
    cond_subpro = (pthread_cond_t *)malloc(sizeof(pthread_cond_t) * lth);
    shared_data_subpro = (int *)malloc(sizeof(int) * lth);
    stList = (int *)malloc(sizeof(int) * lth);
    edList = (int *)malloc(sizeof(int) * lth);
    tid = (pthread_t *)malloc(sizeof(pthread_t) * lth);
    arg = (eig2_args *)malloc(sizeof(eig2_args) * lth);

    for (l = 0; l < lth; l++)
    {
        arg[l].l = l;
        arg[l].cal_pth = &cal_pth;
        arg[l].mutex_subpro = mutex_subpro;
        arg[l].mutex_parent = &mutex_parent;
        arg[l].cond_subpro = cond_subpro;
        arg[l].cond_parent = &cond_parent;
        arg[l].shared_data_subpro = shared_data_subpro;
        arg[l].shared_data_parent = &shared_data_parent;
        arg[l].lth_use = &lth_use;
        arg[l].lth_cal = &lth_cal;

        arg[l].type = &type;
        arg[l].stList = stList;
        arg[l].edList = edList;
        arg[l].a = a;
        arg[l].b = b;
        arg[l].vt = vt;
        arg[l].G = G;
        arg[l].IG = IG;
        arg[l].IGA = IGA;
        arg[l].S = S;
        arg[l].i = &i;
        arg[l].ii = &ii;
        arg[l].j = &j;
        arg[l].n = n;
        arg[l].lt = lt;

        pthread_mutex_init(&mutex_subpro[l], NULL);
        pthread_cond_init(&cond_subpro[l], NULL);
        shared_data_subpro[l] = 0;
    }

    cal_pth = 1;
    shared_data_parent = 0;
    pthread_mutex_init(&mutex_parent, NULL);
    pthread_cond_init(&cond_parent, NULL);

    type = -1;
    lth_cal = 0;
    lth_use = lth;
    for (l = 0; l < lth; l++)
    {
        pthread_create(&tid[l], NULL, eig2_fun, &arg[l]);
    }

    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    *info = 1;
    for (j = 0; j < n; j++)
    {
        s = 0;
        for (k = 0; k <= j - 1; k++)
        {
            s = s + G[j][k] * G[j][k];
        }
        ds = b[j][j] - s;
        if (ds <= 0)
        {
            printf("error_cholesky\n");
            *info = -1;
            return;
        }
        ds = fabs(ds);
        G[j][j] = sqrt(ds);

        allocateThreads(j + 1, n - 1, lth, min_lnum, stList, edList, &lth_use);
        type = 1;
        lth_cal = 0;
        for (l = 0; l < lth_use; l++)
        {
            pthread_mutex_lock(&mutex_subpro[l]);
            shared_data_subpro[l] = 1;
            pthread_cond_signal(&cond_subpro[l]);
            pthread_mutex_unlock(&mutex_subpro[l]);
        }
        pthread_mutex_lock(&mutex_parent);
        while (0 == shared_data_parent)
        {
            pthread_cond_wait(&cond_parent, &mutex_parent);
        }
        shared_data_parent = 0;
        pthread_mutex_unlock(&mutex_parent);

        /*
                for(i=j+1;i<n;i++)
                {
                    s=0;
                    for(k=0;k<=j-1;k++)
                    {
                        s=s+G[i][k]*G[j][k];
                    }
                    G[i][j]=(b[i][j]-s)/G[j][j];
                }
        */
    }

    allocateThreads(0, n - 1, lth, min_lnum, stList, edList, &lth_use);
    type = 2;
    lth_cal = 0;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    allocateThreads(0, n - 1, lth, min_lnum, stList, edList, &lth_use);
    type = 3;
    lth_cal = 0;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    allocateThreads(0, n - 1, lth, min_lnum, stList, edList, &lth_use);
    type = 4;
    lth_cal = 0;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    /*
        for(ii=0;ii<n;ii++)
        {
            for(i=0;i<n;i++)
            {
                IG[i][ii]=0;
            }
            IG[ii][ii]=1;
            for(i=0;i<n;i++)
            {
                for(j=i+1;j<n;j++)
                {
                    IG[j][ii]-=G[j][i]/G[i][i]*IG[i][ii];
                }
            }
            for(i=ii;i<n;i++)
            {
                IG[i][ii]/=G[i][i];
            }
        }

        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                IGA[i][j]=0;
                for(k=0;k<=i;k++)
                {
                    IGA[i][j]+=IG[i][k]*a[k][j];
                }
            }
        }

        for(i=0;i<n;i++)
        {
            for(j=0;j<=i;j++)
            {
                S[i][j]=0;
                for(k=0;k<=j;k++)
                {
                    S[i][j]+=IGA[i][k]*IG[j][k];
                }
                S[j][i]=S[i][j];
            }
        }
    */
    trideigvMul(S, n, d, e, dt, et, lt);

    allocateThreads(0, n - 1, lth, min_lnum, stList, edList, &lth_use);
    type = 5;
    lth_cal = 0;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }
    pthread_mutex_lock(&mutex_parent);
    while (0 == shared_data_parent)
    {
        pthread_cond_wait(&cond_parent, &mutex_parent);
    }
    shared_data_parent = 0;
    pthread_mutex_unlock(&mutex_parent);

    /*
        for(j=0;j<n;j++)
        {
            for(i=0;i<lt;i++)
            {
                vt[i][j]=0;
                for(k=j;k<n;k++)
                {
                    vt[i][j]+=S[i][k]*IG[k][j];
                }
            }
        }
    */

    cal_pth = 0;
    type = 0;
    lth_cal = 0;
    lth_use = lth;
    for (l = 0; l < lth_use; l++)
    {
        pthread_mutex_lock(&mutex_subpro[l]);
        shared_data_subpro[l] = 1;
        pthread_cond_signal(&cond_subpro[l]);
        pthread_mutex_unlock(&mutex_subpro[l]);
    }

    for (l = 0; l < lth; l++)
    {
        pthread_join(tid[l], NULL);
    }

    for (i = 0; i < n; i++)
    {
        free(G[i]);
        free(IG[i]);
        free(IGA[i]);
        free(S[i]);
    }
    free(G);
    free(IG);
    free(IGA);
    free(S);
    free(e);
    free(et);

    for (l = 0; l < lth; l++)
    {
        pthread_mutex_destroy(&mutex_subpro[l]);
        pthread_cond_destroy(&cond_subpro[l]);
    }
    free(mutex_subpro);
    free(cond_subpro);
    free(shared_data_subpro);
    free(stList);
    free(edList);
    free(arg);
    free(tid);
}
