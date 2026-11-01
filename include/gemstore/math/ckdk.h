/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Angular generating function for the three-body GEM.
 * Jacobi maps live in jacobi.c. The kinematic weights and the radial
 * contraction live in kinemtc.c.
 *
 * Y_lm is stored as monomials c (d · n̂)^l. The bra harmonics (l_ρ, l_λ),
 * the ket harmonics, and up to two operator harmonics are multiplied and
 * reduced to one polynomial. Term i is
 *
 *   coe[i] * coef^{nf} * f1^{nf1} * f2^{nf2} * f3^{nf3} * f4^{nf4}
 *          * coeg^{ng} * g1^{ng1} * g2^{ng2} * g3^{ng3} * g4^{ng4}
 *          * h11^{nh11} * ... * h24^{nh24}.
 *
 * nfgh[i] holds those 18 exponents, then one extra slot:
 *   0        nf              power of the pair scale t[0] (coef)
 *   1..4     nf1..nf4        bra ρ, bra λ, ket ρ, ket λ along the pair
 *   5        ng              power of the spectator scale t[5] (coeg)
 *   6..9     ng1..ng4        the same four harmonics along the spectator
 *   10..13   nh11..nh14      first operator harmonic, same four factors
 *   14..17   nh21..nh24      second operator harmonic
 *   18                       radial power of the pair. This file leaves it
 *                            at -1. ckdk_scdk_vtype in spinang.c writes it
 *                            after calc returns, and inteVcenPartA reads it.
 *
 * inteVcenPartA evaluates coe * Π_j t[j]^{nfgh[j]}.
 */

#ifndef GEMSTORE_MATH_CKDK
#define GEMSTORE_MATH_CKDK

#include <pthread.h>

/* Polynomial whose exponents are the 18 kinematic powers above. */
typedef struct {
    int len;
    double *coe;
    int **nfgh;
} ckdk_scdk;

void ckdk_scdk_init(ckdk_scdk *scdk);
void ckdk_scdk_free(ckdk_scdk *scdk);

/* Reduce Y_l1m1(bra ρ) Y_l2m2(bra λ) Y_lh1mh1 Y_lh2mh2 Y_l3m3(ket ρ) Y_l4m4(ket λ)
 * and add cit times that polynomial into *scdk.
 * l1+l2+l3+l4-lh1-lh2 must be even. Central operators use lh = 0,
 * the tensor uses (lh1, lh2) = (2, 0), and spin-orbit uses (1, 1).
 * When lock is 1, mutex covers the Y_lm builds and each push.
 * The final compact always locks mutex, so the pointer must be valid either way. */
void ckdk_scdk_calc(ckdk_scdk *scdk, double cit,
    int l1, int l2, int l3, int l4, int lh1, int lh2,
    int m1, int m2, int m3, int m4, int mh1, int mh2,
    pthread_mutex_t *mutex, int lock);

#endif
