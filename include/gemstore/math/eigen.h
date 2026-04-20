/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * 
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_EIGEN
#define GEMSTORE_MATH_EIGEN

#include <complex.h>

/* Householder tridiagonalization + implicit QR
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * e: Output / working array: subdiagonal elements of the tridiagonal matrix (length ≥ n)
 * et: Output subdiagonal elements corresponding to the selected eigenvalues (length ≥ lt)
 * lt: Input: number of eigenvectors requested (how many columns/rows to keep) */
void eigen_tridiagonal(double **a, int n, double *d, double *e, double *et, int lt);

/* Wrapper: copies matrix, calls eigen_tridiagonal(), extracts eigenvectors
 * a: Input symmetric matrix A (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return */
void eigen_standard(double **a, int n, double *d, double **vt, int lt);

/* Generalized symmetric eigenproblem A x = λ B x using Cholesky + reduction
 * a: Input symmetric matrix A (n × n)
 * b: Input symmetric positive definite matrix B (n × n)
 * n: Dimension of the matrices A and B (order of the problem)
 * d: Output array of eigenvalues (length at least n)
 * vt: Output matrix of selected eigenvectors (lt rows × n columns)
 * lt: Number of largest (or requested) eigenvectors to compute and return
 * (info): Output error code: error_cholesky */
void eigen_general(double **a, double **b, int n, double *d, double **vt, int lt);

/* ============ COMPLEX EIGENPROBLEM SOLVERS (for CGEM) ============ */

/* Complex Hermitian tridiagonalization + implicit QR
 * a: Input complex Hermitian matrix A (n × n)
 * n: Dimension of the matrices
 * d: Output array of REAL eigenvalues (length at least n)
 * e: Output / working array: REAL subdiagonal elements (length ≥ n)
 * et: Output REAL subdiagonal elements for selected eigenvalues
 * lt: Number of eigenvectors requested */
void eigen_tridiagonal_complex(double complex **a, int n, double *d, double *e, double *et, int lt);

/* Wrapper for complex Hermitian eigenproblem: copies matrix, calls tridiagonal, extracts eigenvectors
 * a: Input complex Hermitian matrix A (n × n)
 * n: Dimension
 * d: Output array of REAL eigenvalues
 * vt: Output matrix of selected COMPLEX eigenvectors (lt rows × n columns)
 * lt: Number of eigenvectors */
void eigen_standard_complex(double complex **a, int n, double *d, double complex **vt, int lt);

/* Generalized complex Hermitian eigenproblem: A x = λ B x using complex Cholesky + reduction
 * a: Input complex Hermitian matrix A (n × n)
 * b: Input complex Hermitian positive-definite matrix B (n × n)
 * n: Dimension
 * d: Output array of REAL eigenvalues
 * vt: Output matrix of selected COMPLEX eigenvectors (lt rows × n columns)
 * lt: Number of eigenvectors */
void eigen_general_complex(double complex **a, double complex **b, int n, double *d, double complex **vt, int lt);

#ifdef LAPACKE

void lapack_general(double **a, double **b, int n, double *e, double **vt, int lt);

void lapack_general_complex(double complex **a, double complex **b, int n, double *e, double complex **vt, int lt);

#endif

#endif