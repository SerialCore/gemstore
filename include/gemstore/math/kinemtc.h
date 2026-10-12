/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Kinematic weights of one pair matrix element between two GEM Gaussians.
 * The angular polynomial is sumckdk. The pair potential is the callback iv.
 *
 * A basis Gaussian is written in its own Jacobi chart. The operator acts
 * on the separation r of channel c, with spectator coordinate R. getT*
 * maps the bra chart a = qf.c and the ket chart b = qi.c onto c, completes
 * the square, and shifts R. The spectator integral is then analytic and the
 * matrix element factors as
 *
 *   Σ coe * (kinematic powers) * ∫ r^n exp(-b11 r^2) V(r) dr.
 *
 * t[] carries those weights:
 *   t[0 .. lent)   bases raised to nfgh[0 .. lent)
 *   t[18]          b11, the coefficient of r^2 after the shift
 *   t[19]          constant produced by the R integral
 *   t[20]          product of the four GEM normalizations
 * lent is 10 for a central operator, 14 for a tensor, 18 for spin-orbit.
 * nfgh[18] is the radial power n, written by spinang, not by getT*.
 *
 * Which vector is expanded:
 *   t1r   jacobi_3b_r_1    pair separation ρ
 *   tir   r_1, then quark i of the pair
 *   tjr   r_1, then quark j of the pair
 *   t2r   jacobi_3b_r_2    spectator λ
 *   t1p   jacobi_3b_p_1    momentum conjugate to ρ
 *   t2p   jacobi_3b_p_2    momentum conjugate to λ
 *   tpi   jacobi_3b_p_i    momentum of one quark
 *
 * The vt argument of getT* chooses the operator vectors stored in t[10..17]:
 *   1 central, 2 tensor, 3 spin-orbit on quark i, 4 on quark j,
 *   5 along ρ, 6 along λ.
 * A wrapper only passes the values that its getT* implements.
 */

#ifndef GEMSTORE_MATH_KINEMTC
#define GEMSTORE_MATH_KINEMTC

#include <gemstore/basis/basis.h>
#include <gemstore/math/sumckdk.h>

/* Fills t[] for one bra/ket pair. lent is the number of polynomial bases. */
typedef void (*txrp_vtype)(basis_3b, basis_3b, double *, int *, int);

/* ∫ r^n exp(-b11 r^2) V(r) dr. ctx is the model parameter block. */
typedef double (*kin_radial_fn)(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);

void t1r_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t2r_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1p_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t2p_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tpi_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_tens(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_soii(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_soij(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_soji(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_sojj(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_sorr(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void t1r_sorl(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tir_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tir_tens(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tir_soii(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tir_soij(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tjr_cent(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tjr_tens(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tjr_sojj(basis_3b qf, basis_3b qi, double *t, int *lent, int c);
void tjr_soji(basis_3b qf, basis_3b qi, double *t, int *lent, int c);

/* Contract the angular polynomial with t[] and the radial callback. */
double inteVcenPartA(txrp_vtype get_trxp_vtype, sumckdk_scdk scdk, basis_3b qf, basis_3b qi,
    kin_radial_fn iv, void *ctx, int c);

#endif
