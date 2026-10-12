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
 *                            at -1. sumckdk_scdk_vtype in spinang.c writes it
 *                            after calc returns, and inteVcenPartA reads it.
 *
 * inteVcenPartA evaluates coe * Π_j t[j]^{nfgh[j]}. sumckdk_scdk_vals is the
 * recycle evaluator and also multiplies by sin(nf) cos(ng); the baryon
 * contraction does not call it.
 */

#ifndef GEMSTORE_MATH_SUMCKDK
#define GEMSTORE_MATH_SUMCKDK

#include <pthread.h>

/* One j term in the power series of Y_lm.
 * factorial is n! for n = 0..170, stored as double. No range check. */
double almj(int l, int m, int j);

/* C(n, m) = n! / ((n-m)! m!) from that same table. Needs 0 <= m <= n <= 170. */
double binomial(int n, int m);

/* One expansion of Y_lm. Term i is c[i] (d[i] · n̂)^l, with d[i] = (x, y, z). */
typedef struct {
    int len;
    double *c;
    double **d;
} sumckdk_cdlm;

void sumckdk_cdlm_init(sumckdk_cdlm *cdlm);
void sumckdk_cdlm_free(sumckdk_cdlm *cdlm);
void sumckdk_cdlm_push(sumckdk_cdlm *cdlm, double c, double dx, double dy, double dz);

/* Ket expansion of Y_lm. l = m = 0 is the constant 1/sqrt(4π). */
void sumckdk_cdlm_cali(sumckdk_cdlm *cdlm, int l, int m);

/* Bra expansion. Same coefficients as cali; d_y has the opposite sign. */
void sumckdk_cdlm_calf(sumckdk_cdlm *cdlm, int l, int m);

/* Polynomial whose exponents are the 18 kinematic powers above. */
typedef struct {
    int len;
    double *coe;
    int **nfgh;
} sumckdk_scdk;

void sumckdk_scdk_init(sumckdk_scdk *scdk);
void sumckdk_scdk_free(sumckdk_scdk *scdk);

/* Add coe onto an existing exponent vector, or append a new term.
 * The new term's nfgh[18] is -1. */
void sumckdk_scdk_push(sumckdk_scdk *scdk, double coe,
    int nf, int nf1, int nf2, int nf3, int nf4,
    int ng, int ng1, int ng2, int ng3, int ng4,
    int nh11, int nh12, int nh13, int nh14,
    int nh21, int nh22, int nh23, int nh24);

/* Drop terms with |coe| <= 2^-53 and compact the list.
 * Exponents 0..17 are copied. Slot 18 is not. */
void sumckdk_scdk_spfy(sumckdk_scdk *scdk);

/* Reduce Y_l1m1(bra ρ) Y_l2m2(bra λ) Y_lh1mh1 Y_lh2mh2 Y_l3m3(ket ρ) Y_l4m4(ket λ)
 * and add cit times that polynomial into *scdk.
 * l1+l2+l3+l4-lh1-lh2 must be even. Central operators use lh = 0,
 * the tensor uses (lh1, lh2) = (2, 0), and spin-orbit uses (1, 1).
 * When lock is 1, mutex covers the Y_lm builds and each push.
 * The final compact always locks mutex, so the pointer must be valid either way. */
void sumckdk_scdk_calc(sumckdk_scdk *scdk, double cit,
    int l1, int l2, int l3, int l4, int lh1, int lh2,
    int m1, int m2, int m3, int m4, int mh1, int mh2,
    pthread_mutex_t *mutex, int lock);

/* Recycle evaluator: the 18 powers, times sin(nf) cos(ng).
 * Not used by inteVcenPartA. */
double sumckdk_scdk_vals(sumckdk_scdk scdk,
    double coef, double f1, double f2, double f3, double f4,
    double coeg, double g1, double g2, double g3, double g4,
    double h11, double h12, double h13, double h14,
    double h21, double h22, double h23, double h24);

#endif
