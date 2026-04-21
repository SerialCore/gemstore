/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/gaussnode.h>

#include <stdlib.h>
#include <math.h>
#include <string.h>

/** @addtogroup GaussNode
 * @{
 */

/**
 * Evaluate Hermite polynomial H_n(x) using recurrence.
 * Recurrence: H_{n+1}(x) = 2x*H_n(x) - 2n*H_{n-1}(x)
 */
long double gaussnode_hermite_poly(int n, long double x)
{
    if (n < 0) return 0.0L;
    if (n == 0) return 1.0L;
    if (n == 1) return 2.0L * x;
    
    long double h_prev2 = 1.0L;      /* H_0(x) */
    long double h_prev1 = 2.0L * x;  /* H_1(x) */
    
    for (int k = 1; k < n; k++) {
        long double h_curr = 2.0L * x * h_prev1 - 2.0L * k * h_prev2;
        h_prev2 = h_prev1;
        h_prev1 = h_curr;
    }
    
    return h_prev1;
}

/**
 * Evaluate derivative of Hermite polynomial: H'_n(x) = 2n * H_{n-1}(x)
 */
long double gaussnode_hermite_poly_derivative(int n, long double x)
{
    if (n == 0) return 0.0L;
    return 2.0L * n * gaussnode_hermite_poly(n - 1, x);
}

/**
 * Compute moments: m_k = ∫₀^∞ x^k exp(-x²)dx
 * m_0 = √π/2, m_1 = 1/2, m_k = (k-1)/2 * m_{k-2}
 */
int gaussnode_compute_moments(int max_order, long double *moments)
{
    if (moments == NULL || max_order < 0) {
        return -1;
    }
    
    moments[0] = 1.7724538509055160272981674833411L;  /* √π/2 */
    if (max_order >= 1) {
        moments[1] = 0.5L;
    }
    
    for (int k = 2; k <= max_order; k++) {
        moments[k] = (k - 1) * 0.5L * moments[k - 2];
    }
    
    return 0;
}

/**
 * Find root of H_n(x) using Newton-Raphson method
 */
int gaussnode_find_root(
    long double *x,
    int degree,
    int max_iter,
    long double tol)
{
    if (x == NULL || degree <= 0) {
        return -1;
    }
    
    for (int iter = 0; iter < max_iter; iter++) {
        long double f = gaussnode_hermite_poly(degree, *x);
        long double df = gaussnode_hermite_poly_derivative(degree, *x);
        
        if (fabsl(df) < 1.0e-30L) {
            return -1;
        }
        
        long double dx = f / df;
        *x = *x - dx;
        
        if (fabsl(dx) < tol) {
            return 0;
        }
    }
    
    return -1;
}

/**
 * Solve Vandermonde system for quadrature weights
 */
int gaussnode_solve_vandermonde(
    int n,
    const long double *nodes,
    const long double *moments,
    long double *weights)
{
    if (n <= 0 || nodes == NULL || moments == NULL || weights == NULL) {
        return -1;
    }
    
    /* Build and solve Vandermonde system Ax=b via LU factorization */
    long double *A = (long double *)malloc(n * n * sizeof(long double));
    long double *b = (long double *)malloc(n * sizeof(long double));
    int *piv = (int *)malloc(n * sizeof(int));
    
    if (A == NULL || b == NULL || piv == NULL) {
        free(A);
        free(b);
        free(piv);
        return -1;
    }
    
    /* Build Vandermonde matrix A[i][j] = x_i^j */
    for (int i = 0; i < n; i++) {
        long double power = 1.0L;
        for (int j = 0; j < n; j++) {
            A[i * n + j] = power;
            power *= nodes[i];
        }
        b[i] = moments[i];
    }
    
    /* LU factorization with partial pivoting */
    for (int k = 0; k < n; k++) {
        int pivot_row = k;
        long double max_val = fabsl(A[k * n + k]);
        
        for (int i = k + 1; i < n; i++) {
            if (fabsl(A[i * n + k]) > max_val) {
                max_val = fabsl(A[i * n + k]);
                pivot_row = i;
            }
        }
        
        if (max_val < 1.0e-30L) {
            free(A);
            free(b);
            free(piv);
            return -1;
        }
        
        piv[k] = pivot_row;
        if (pivot_row != k) {
            for (int j = k; j < n; j++) {
                long double temp = A[k * n + j];
                A[k * n + j] = A[pivot_row * n + j];
                A[pivot_row * n + j] = temp;
            }
            long double temp = b[k];
            b[k] = b[pivot_row];
            b[pivot_row] = temp;
        }
        
        /* Eliminate column */
        for (int i = k + 1; i < n; i++) {
            long double factor = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        long double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * weights[j];
        }
        weights[i] = sum / A[i * n + i];
    }
    
    free(A);
    free(b);
    free(piv);
    
    return 0;
}

/**
 * Robust root finding using bracketing + refinement
 */
static int gaussnode_find_roots_robust(
    int n,
    long double **nodes,
    int max_iter,
    long double tol)
{
    *nodes = (long double *)malloc(n * sizeof(long double));
    if (*nodes == NULL) {
        return -1;
    }
    
    int roots_found = 0;
    long double x_max = sqrtl((long double)n) + 5.0L;
    int n_samples = 500 * n;
    
    long double prev_x = -x_max;
    long double prev_f = gaussnode_hermite_poly(n, prev_x);
    
    for (int sample = 1; sample < n_samples && roots_found < n; sample++) {
        long double x = -x_max + 2.0L * x_max * (long double)sample / (long double)n_samples;
        long double f = gaussnode_hermite_poly(n, x);
        
        if (prev_f * f < 0.0L) {
            long double x_root = 0.5L * (prev_x + x);
            
            if (gaussnode_find_root(&x_root, n, max_iter, tol) == 0) {
                /* Take absolute value for half-range integral */
                x_root = fabsl(x_root);
                
                /* Check for duplicates */
                int is_dup = 0;
                for (int j = 0; j < roots_found; j++) {
                    if (fabsl((*nodes)[j] - x_root) < 10.0L * tol) {
                        is_dup = 1;
                        break;
                    }
                }
                
                if (!is_dup) {
                    (*nodes)[roots_found++] = x_root;
                }
            }
        }
        
        prev_x = x;
        prev_f = f;
    }
    
    if (roots_found < n) {
        free(*nodes);
        *nodes = NULL;
        return -1;
    }
    
    /* Sort */
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if ((*nodes)[i] > (*nodes)[j]) {
                long double temp = (*nodes)[i];
                (*nodes)[i] = (*nodes)[j];
                (*nodes)[j] = temp;
            }
        }
    }
    
    return 0;
}

/* Public API */

int gaussnode_compute(int n, long double **nodes, long double **weights)
{
    return gaussnode_compute_advanced(n, 1000, 1.0e-18L, nodes, weights);
}

int gaussnode_compute_advanced(
    int n,
    int max_iter,
    long double tol,
    long double **nodes,
    long double **weights)
{
    if (n <= 0 || nodes == NULL || weights == NULL) {
        return -1;
    }
    
    *weights = (long double *)malloc(n * sizeof(long double));
    long double *moments = (long double *)malloc((n + 1) * sizeof(long double));
    
    if (*weights == NULL || moments == NULL) {
        free(*weights);
        free(moments);
        *nodes = NULL;
        *weights = NULL;
        return -1;
    }
    
    if (gaussnode_compute_moments(n, moments) != 0) {
        free(*weights);
        free(moments);
        *weights = NULL;
        return -1;
    }
    
    if (gaussnode_find_roots_robust(n, nodes, max_iter, tol) != 0) {
        free(*weights);
        free(moments);
        *weights = NULL;
        return -2;
    }
    
    if (gaussnode_solve_vandermonde(n, *nodes, moments, *weights) != 0) {
        free(*nodes);
        free(*weights);
        free(moments);
        *nodes = NULL;
        *weights = NULL;
        return -1;
    }
    
    free(moments);
    return 0;
}

void gaussnode_free(long double *nodes, long double *weights)
{
    free(nodes);
    free(weights);
}

/** @} */
