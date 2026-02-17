/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_NUMERICAL_EIGEN
#define GEMSTORE_NUMERICAL_EIGEN

/* Single-threaded version: Householder tridiagonalization + implicit QR
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * e: Output / working array: subdiagonal elements of the tridiagonal matrix (length ≥ n)
 * et: Output subdiagonal elements corresponding to the selected eigenvalues (length ≥ lt)
 * lt: Input: number of eigenvectors requested (how many columns/rows to keep) */
void eigen_tridiagonal(double **a, int n, double *d, double *e, double *et, int lt);

/* Multi-threaded version of eigen_tridiagonal()
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * e: Output / working array: subdiagonal elements of the tridiagonal matrix (length ≥ n)
 * et: Output subdiagonal elements corresponding to the selected eigenvalues (length ≥ lt)
 * lt: Input: number of eigenvectors requested (how many columns/rows to keep) */
void eigen_tridiagonal_thread(double **a, int n, double *d, double *e, double *et, int lt);

/* Wrapper: copies matrix, calls eigen_tridiagonal(), extracts eigenvectors
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return */
void eigen_standard(double **a, int n, double *d, double **vt, int lt);

/* Multi-threaded version of eigen_standard()
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return */
void eigen_standard_thread(double **a, int n, double *d, double **vt, int lt);

/* Generalized symmetric eigenproblem A x = λ B x using Cholesky + reduction
 * a: Input symmetric matrix A (n × n)
 * b: Input symmetric positive definite matrix B (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return
 * (info): Output error code: error_cholesky */
void eigen_general(double **a, double **b, int n, double *d, double **vt, int lt);

/* Multi-threaded version of eigen_general()
 * a: Input symmetric matrix A (n × n)
 * b: Input symmetric positive definite matrix B (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return
 * (info): Output error code: error_cholesky */
void eigen_general_thread(double **a, double **b, int n, double *d, double **vt, int lt);

#endif
