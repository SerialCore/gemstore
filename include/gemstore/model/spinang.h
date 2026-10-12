/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Spin-angular coupling of one quark pair (recycle/vtype.h).
 * V(r) and the Gaussian widths are not used here. The radial integral
 * is giradial.c, and the polynomial itself is built by sumckdk.
 *
 * An mlsj row is (coefficient, m_s1, m_s2, m_s3, m_ρ, m_λ), with quark
 * labels 1, 2, 3. sumckdk_scdk_vtype permutes the spin projections onto
 * channel c before calling the operator.
 *
 * l1, l2, ml1, ml2 are the bra l_ρ, l_λ and their projections.
 * l3, l4, ml3, ml4 are the ket. The operator harmonics passed to sumckdk
 * are (0,0) for a central or contact term, (2,0) for the tensor, and
 * (1,1) for spin-orbit.
 *
 * fg then writes the radial power nfgh[18] and a factorial weight.
 * nf = nfgh[0] is the pair-contraction count, mg = nfgh[5] the spectator count:
 *   1  central and contact    power 2*nf+2,  weight 1/(2*nf+1)! / mg!
 *   2  tensor                 power 2*nf+4,  weight (nf+2)! / ((2*nf+5)! nf! mg!)
 *   3  Vso^{ii}
 *   4  Vso^{jj}               power 2*nf+4,  weight (nf+1)! / ((2*nf+3)! nf! mg!)
 *   5  Vso^{ji}
 *   6  Vso^{ij}
 */

#ifndef GEMSTORE_MODEL_SPINANG
#define GEMSTORE_MODEL_SPINANG

#include <gemstore/basis/basis.h>
#include <gemstore/math/sumckdk.h>
#include <gemstore/math/matrix.h>

#include <pthread.h>

/* One spin-angular operator. The nine doubles are the bra and ket CG weights
 * and the six spin projections on channel c. */
typedef void (*vtype)(sumckdk_scdk *, double, double, double, double, double, double, double, double,
    int, int, int, int, int, int, int, int, basis_3b, basis_3b, pthread_mutex_t *, int, int);

/* Sum operator over the mlsj projections of bra (nf,nfp) and ket (ni,nip),
 * then apply the fg weight. mutex is the rule documented on sumckdk_scdk_calc. */
void sumckdk_scdk_vtype(sumckdk_scdk *scdk, vtype vsodt, int fg, bases_3b qnlist, matrix_t **mlsj,
    int nf, int nfp, int ni, int nip, pthread_mutex_t *mutex, int lock, int c);

/* Spin scalar. All three m_s must match. */
void vcent(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Contact spin-spin. The spectator projection m_sk is unchanged. */
void vcont(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Tensor. Operator harmonic of rank 2 on the pair direction. */
void vtens(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Vector spin-orbit on quark i. Two rank-1 harmonics. */
void vsoii(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Vso^{ji}. The CG is the same as vsoii. cbaryon distinguishes it by
 * calling t1r_soji rather than tir_soii. */
void vsoji(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Vector spin-orbit on quark j. */
void vsojj(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);
/* Vso^{ij}. The CG is the same as vsojj. cbaryon pairs it with tir_soij. */
void vsoij(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c);

#endif
