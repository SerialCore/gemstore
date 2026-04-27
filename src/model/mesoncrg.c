/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/mesoncrg.h>
#include <gemstore/model/gimodel.h>

#include <gemstore/basis/orbit.h>

#include <gemstore/math/soc.h>
#include <gemstore/math/ceigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/cmatrix.h>
#include <gemstore/math/integral.h>

#include <gemstore/param/argset.h>

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

static inline double getmq(int index, const argsGIModel_t *args_model)
{
    switch (index) {
        case 1: return args_model->mn;
        case 2: return args_model->ms;
        case 3: return args_model->mc;
        case 4: return args_model->mb;
        default: return args_model->mn;
    }
}

void spectra_meson_CRG(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    double omega = args_input->omega;
    int f1 = args_input->f1, f2 = args_input->f2;
    double S = args_input->S, L = args_input->L, J = args_input->J;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
        basis[i].param = omega;
    }

    /* construct matrices */
    cmatrix_t mT = cmatrix_init(nmax, nmax);
    cmatrix_t mbetaijCoul = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltaijCont = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltaiiSov = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltajjSov = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltaijSov = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltaiiSos = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltajjSos = cmatrix_init(nmax, nmax);
    cmatrix_t mdeltaijTens = cmatrix_init(nmax, nmax);
    cmatrix_t mVcoul = cmatrix_init(nmax, nmax);
    cmatrix_t mVconf = cmatrix_init(nmax, nmax);
    cmatrix_t mVcont = cmatrix_init(nmax, nmax);
    cmatrix_t mVsovi = cmatrix_init(nmax, nmax);
    cmatrix_t mVsovj = cmatrix_init(nmax, nmax);
    cmatrix_t mVsovij = cmatrix_init(nmax, nmax);
    cmatrix_t mVsosi = cmatrix_init(nmax, nmax);
    cmatrix_t mVsosj = cmatrix_init(nmax, nmax);
    cmatrix_t mVtens = cmatrix_init(nmax, nmax);
    cmatrix_t tmT = cmatrix_init(nmax, nmax);
    cmatrix_t tmbetaijCoul = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltaijCont = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltaiiSov = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltajjSov = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltaijSov = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltaiiSos = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltajjSos = cmatrix_init(nmax, nmax);
    cmatrix_t tmdeltaijTens = cmatrix_init(nmax, nmax);
    cmatrix_t tmVcoul = cmatrix_init(nmax, nmax);
    cmatrix_t tmVconf = cmatrix_init(nmax, nmax);
    cmatrix_t tmVcont = cmatrix_init(nmax, nmax);
    cmatrix_t tmVsovi = cmatrix_init(nmax, nmax);
    cmatrix_t tmVsovj = cmatrix_init(nmax, nmax);
    cmatrix_t tmVsovij = cmatrix_init(nmax, nmax);
    cmatrix_t tmVsosi = cmatrix_init(nmax, nmax);
    cmatrix_t tmVsosj = cmatrix_init(nmax, nmax);
    cmatrix_t tmVtens = cmatrix_init(nmax, nmax);
    cmatrix_t Hfi = cmatrix_init(nmax, nmax);
    cmatrix_t Nfi = cmatrix_init(nmax, nmax);

    /* prepare variables */
    double factor_r;
    double factor_p;
    double s1 = 0.5, s2 = 0.5;
    double m1 = getmq(f1, args_model);
    double m2 = getmq(f2, args_model);
    double C12 = -4.0 / 3.0;
    double sigmaij = sigma_ij(m1, m2, args_model->sigma_0, args_model->s);
    args_dynmc->mi = m1;
    args_dynmc->mj = m2;
    args_dynmc->Cij = C12;
    args_dynmc->Sigij = sigmaij;
    sigma_k_ij(sigmaij, args_dynmc->Sigkij);

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor_r = 1 / sqrt(basis[i].scale + basis[j].scale);
            factor_p = sqrt(4 * basis[i].scale * basis[j].scale * (1 + omega * omega) / (basis[i].scale + basis[j].scale));

            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_crg_hamilton(CGRnlp, GIVt, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mbetaijCoul.value[i][j] = integral_crg_hamilton(CGRnlp, GIVbetaijcoul, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijCont.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltaijcont, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSov.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltaiisov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSov.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltajjsov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijSov.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltaijsov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSos.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltaiisos, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSos.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltajjsos, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijTens.value[i][j] = integral_crg_hamilton(CGRnlp, GIVdeltaijtens, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mVcoul.value[i][j] = integral_crg_hamilton(CGRnlr, GIVcoul, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVconf.value[i][j] = integral_crg_hamilton(CGRnlr, GIVconf, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVcont.value[i][j] = integral_crg_hamilton(CGRnlr, GIVcont, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovi.value[i][j] = integral_crg_hamilton(CGRnlr, GIVsovi, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovj.value[i][j] = integral_crg_hamilton(CGRnlr, GIVsovj, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovij.value[i][j] = integral_crg_hamilton(CGRnlr, GIVsovij, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosi.value[i][j] = integral_crg_hamilton(CGRnlr, GIVsosi, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosj.value[i][j] = integral_crg_hamilton(CGRnlr, GIVsosj, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVtens.value[i][j] = integral_crg_hamilton(CGRnlr, GIVtens, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            Nfi.value[i][j] = integral_crg_overlap(CGRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    /* Cholesky decomposition */
    cmatrix_t mL = cmatrix_init(nmax, nmax);
    cmatrix_t mLinv = cmatrix_init(nmax, nmax);
    cmatrix_cholesky_decomp(&Nfi, &mL);          /* S = L * L^T */
    cmatrix_inverse_lowertri(&mL, &mLinv);       /* Linv = L^{-1} */

    /* transform Hamiltonian matrices in orthogonal basis */
    cmatrix_productT(&mLinv, &mT, &tmT);
    cmatrix_productT(&mLinv, &mbetaijCoul, &tmbetaijCoul);
    cmatrix_productT(&mLinv, &mdeltaijCont, &tmdeltaijCont);
    cmatrix_productT(&mLinv, &mdeltaiiSov, &tmdeltaiiSov);
    cmatrix_productT(&mLinv, &mdeltajjSov, &tmdeltajjSov);
    cmatrix_productT(&mLinv, &mdeltaijSov, &tmdeltaijSov);
    cmatrix_productT(&mLinv, &mdeltaiiSos, &tmdeltaiiSos);
    cmatrix_productT(&mLinv, &mdeltajjSos, &tmdeltajjSos);
    cmatrix_productT(&mLinv, &mdeltaijTens, &tmdeltaijTens);
    cmatrix_productT(&mLinv, &mVcoul, &tmVcoul);
    cmatrix_productT(&mLinv, &mVconf, &tmVconf);
    cmatrix_productT(&mLinv, &mVcont, &tmVcont);
    cmatrix_productT(&mLinv, &mVsovi, &tmVsovi);
    cmatrix_productT(&mLinv, &mVsovj, &tmVsovj);
    cmatrix_productT(&mLinv, &mVsovij, &tmVsovij);
    cmatrix_productT(&mLinv, &mVsosi, &tmVsosi);
    cmatrix_productT(&mLinv, &mVsosj, &tmVsosj);
    cmatrix_productT(&mLinv, &mVtens, &tmVtens);

    /* construct Hamiltonian matrix */
    cmatrix_t temp = cmatrix_init(nmax, nmax);
    cmatrix_sum(&tmT, &tmVconf, &Hfi);
    cmatrix_productT(&tmbetaijCoul, &tmVcoul, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltaijCont, &tmVcont, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltaiiSov, &tmVsovi, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltajjSov, &tmVsovj, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltaijSov, &tmVsovij, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltaiiSos, &tmVsosi, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltajjSos, &tmVsosj, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&tmdeltaijTens, &tmVtens, &temp);
    cmatrix_sum(&Hfi, &temp, &Hfi);
    cmatrix_productT(&mLinv, &Nfi, &temp);

/* In GI model, the final basis should be orthogonal.
 * Therefore this could be general eigen system with orthogonal basis and diagonal Nfi,
 * or just standard eigen system directly. */
#ifdef LAPACKE
    lapack_general_complex(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, v_len);
#else
    eigen_general_complex(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, v_len);
#endif

    free(basis);
    cmatrix_free(&temp);
    cmatrix_free(&mT);
    cmatrix_free(&mbetaijCoul);
    cmatrix_free(&mdeltaijCont);
    cmatrix_free(&mdeltaiiSov);
    cmatrix_free(&mdeltajjSov);
    cmatrix_free(&mdeltaijSov);
    cmatrix_free(&mdeltaiiSos);
    cmatrix_free(&mdeltajjSos);
    cmatrix_free(&mdeltaijTens);
    cmatrix_free(&mVcoul);
    cmatrix_free(&mVconf);
    cmatrix_free(&mVcont);
    cmatrix_free(&mVsovi);
    cmatrix_free(&mVsovj);
    cmatrix_free(&mVsovij);
    cmatrix_free(&mVsosi);
    cmatrix_free(&mVsosj);
    cmatrix_free(&mVtens);
    cmatrix_free(&tmT);
    cmatrix_free(&tmbetaijCoul);
    cmatrix_free(&tmdeltaijCont);
    cmatrix_free(&tmdeltaiiSov);
    cmatrix_free(&tmdeltajjSov);
    cmatrix_free(&tmdeltaijSov);
    cmatrix_free(&tmdeltaiiSos);
    cmatrix_free(&tmdeltajjSos);
    cmatrix_free(&tmdeltaijTens);
    cmatrix_free(&tmVcoul);
    cmatrix_free(&tmVconf);
    cmatrix_free(&tmVcont);
    cmatrix_free(&tmVsovi);
    cmatrix_free(&tmVsovj);
    cmatrix_free(&tmVsovij);
    cmatrix_free(&tmVsosi);
    cmatrix_free(&tmVsosj);
    cmatrix_free(&tmVtens);
    cmatrix_free(&Hfi);
    cmatrix_free(&Nfi);
}

void radius_meson_CRG(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len)
{
    int nmax = input->nmax;
    double rmax = input->rmax;
    double rmin = input->rmin;
    int L = (int)input->L;
    double omega = input->omega;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
        basis[i].param = omega;
    }

    /* construct matrices */
    matrix_t mR2 = matrix_init(nmax, nmax);
    matrix_t mOver = matrix_init(nmax, nmax);

    /* prepare variables */
    double factor;
    double coef;
    complex r2sum;
    complex oversum;
    double norm;
    complex rms2;
    double fm = 5.06773093854369882649;

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);
            mR2.value[i][j] = creal(integral_crg_radius(CGRnlr, factor, &basis[i], &basis[j]));
            mOver.value[i][j] = creal(integral_crg_overlap(CGRnlr, factor, &basis[i], &basis[j]));
        }
    }

    /* calculate rms radius with orthogonalized coefficients */
    for (int n = 0; n < len; n++) {
        r2sum = 0.0 + 0.0 * I;
        oversum = 0.0 + 0.0 * I;
        norm = 0.0;

        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                coef = vector->value[n][i] * vector->value[n][j];
                r2sum += coef * mR2.value[i][j];
                oversum += coef * mOver.value[i][j];
            }
            norm += vector->value[n][i] * vector->value[n][i];
        }

        rms2 = (norm > 1e-12) ? r2sum / oversum : 0.0 + 0.0 * I;
        radius->value[n] = sqrt(cabs(rms2)) / fm;
    }

    free(basis);
    matrix_free(&mR2);
    matrix_free(&mOver);
}