/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/threebody.h>
#include <gemstore/math/eigen.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int threebody_pair_identical(int id1, int id2, int id3, int c)
{
    switch (c) {
        case 1: return id1 == id2;
        case 2: return id3 == id1;
        case 3: return id2 == id3;
        default: return 0;
    }
}

double threebody_exchange_eta(int f12, double sij, int lrho)
{
    int phase = ((1 + (int)sij + lrho) % 2 == 0) ? 1 : -1;
    return (double)(f12 * phase);
}

void threebody_scdk_table_alloc(int *len_part, int len_list, sumckdk_scdk *****tab)
{
    *tab = (sumckdk_scdk ****)malloc(sizeof(sumckdk_scdk ***) * (size_t)len_list);
    for (int i = 0; i < len_list; i++) {
        (*tab)[i] = (sumckdk_scdk ***)malloc(sizeof(sumckdk_scdk **) * (size_t)len_part[i]);
        for (int j = 0; j < len_part[i]; j++) {
            (*tab)[i][j] = (sumckdk_scdk **)malloc(sizeof(sumckdk_scdk *) * (size_t)len_list);
            for (int k = 0; k < len_list; k++) {
                (*tab)[i][j][k] = (sumckdk_scdk *)malloc(sizeof(sumckdk_scdk) * (size_t)len_part[k]);
            }
        }
    }
}

void threebody_scdk_table_free(int *len_part, int len_list, sumckdk_scdk *****tab)
{
    for (int i = 0; i < len_list; i++) {
        for (int j = 0; j < len_part[i]; j++) {
            for (int k = 0; k < len_list; k++) {
                for (int l = 0; l < len_part[k]; l++) {
                    sumckdk_scdk_free(&((*tab)[i][j][k][l]));
                }
                free((*tab)[i][j][k]);
            }
            free((*tab)[i][j]);
        }
        free((*tab)[i]);
    }
    free(*tab);
    *tab = NULL;
}

int threebody_overlap_basis(const matrix_t *N, double rel_thresh, matrix_t *vt)
{
    int n = N->row;
    matrix_t ncopy = matrix_init(n, n);
    matrix_copy(&ncopy, N);

    double ntrace = 0.0;
    for (int i = 0; i < n; i++) {
        ntrace += N->value[i][i];
    }
    double ridge = 1e-12 * (fabs(ntrace) / (n > 0 ? n : 1) + 1.0);
    for (int i = 0; i < n; i++) {
        ncopy.value[i][i] += ridge;
    }

    array_t n_e = array_init(n);
    matrix_t n_vec = matrix_init(n, n);
    eigen_standard(ncopy.value, n, n_e.value, n_vec.value, n);

    double n_emax = 0.0;
    for (int i = 0; i < n; i++) {
        if (n_e.value[i] > n_emax) {
            n_emax = n_e.value[i];
        }
    }
    double n_thresh = rel_thresh * (n_emax > 0.0 ? n_emax : 1.0);
    if (n_thresh < 1e-12) {
        n_thresh = 1e-12;
    }

    int nkeep = 0;
    for (int i = 0; i < n; i++) {
        if (n_e.value[i] > n_thresh) {
            nkeep++;
        }
    }
    if (nkeep <= 0) {
        matrix_free(&ncopy);
        matrix_free(&n_vec);
        array_free(&n_e);
        vt->value = NULL;
        vt->row = 0;
        vt->col = 0;
        return 0;
    }

    printf("  overlap N: rank %d / %d  (drop λ ≤ %.2e, λ_max=%.3e)\n",
        nkeep, n, n_thresh, n_emax);
    fflush(stdout);

    *vt = matrix_init(nkeep, n);
    int row = 0;
    for (int i = 0; i < n; i++) {
        if (n_e.value[i] <= n_thresh) {
            continue;
        }
        double scale = 1.0 / sqrt(n_e.value[i]);
        for (int j = 0; j < n; j++) {
            vt->value[row][j] = n_vec.value[i][j] * scale;
        }
        row++;
    }

    matrix_free(&ncopy);
    matrix_free(&n_vec);
    array_free(&n_e);
    return nkeep;
}
