/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/mesongem.h>
#include <gemstore/model/gimodel.h>

#include <gemstore/basis/orbit.h>

#include <gemstore/math/soc.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>

#include <gemstore/param/argset.h>

#include <stdio.h>
#include <stdlib.h>

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

void spectra_meson_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    int f1 = args_input->f1, f2 = args_input->f2;
    double S = args_input->S, L = args_input->L, J = args_input->J;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
    }

    /* construct matrices */
    matrix_t mT = matrix_init(nmax, nmax);
    matrix_t mbetaijCoul = matrix_init(nmax, nmax);
    matrix_t mdeltaijCont = matrix_init(nmax, nmax);
    matrix_t mdeltaiiSov = matrix_init(nmax, nmax);
    matrix_t mdeltajjSov = matrix_init(nmax, nmax);
    matrix_t mdeltaijSov = matrix_init(nmax, nmax);
    matrix_t mdeltaiiSos = matrix_init(nmax, nmax);
    matrix_t mdeltajjSos = matrix_init(nmax, nmax);
    matrix_t mdeltaijTens = matrix_init(nmax, nmax);
    matrix_t mVcoul = matrix_init(nmax, nmax);
    matrix_t mVconf = matrix_init(nmax, nmax);
    matrix_t mVcont = matrix_init(nmax, nmax);
    matrix_t mVsovi = matrix_init(nmax, nmax);
    matrix_t mVsovj = matrix_init(nmax, nmax);
    matrix_t mVsovij = matrix_init(nmax, nmax);
    matrix_t mVsosi = matrix_init(nmax, nmax);
    matrix_t mVsosj = matrix_init(nmax, nmax);
    matrix_t mVtens = matrix_init(nmax, nmax);
    matrix_t tmT = matrix_init(nmax, nmax);
    matrix_t tmbetaijCoul = matrix_init(nmax, nmax);
    matrix_t tmdeltaijCont = matrix_init(nmax, nmax);
    matrix_t tmdeltaiiSov = matrix_init(nmax, nmax);
    matrix_t tmdeltajjSov = matrix_init(nmax, nmax);
    matrix_t tmdeltaijSov = matrix_init(nmax, nmax);
    matrix_t tmdeltaiiSos = matrix_init(nmax, nmax);
    matrix_t tmdeltajjSos = matrix_init(nmax, nmax);
    matrix_t tmdeltaijTens = matrix_init(nmax, nmax);
    matrix_t tmVcoul = matrix_init(nmax, nmax);
    matrix_t tmVconf = matrix_init(nmax, nmax);
    matrix_t tmVcont = matrix_init(nmax, nmax);
    matrix_t tmVsovi = matrix_init(nmax, nmax);
    matrix_t tmVsovj = matrix_init(nmax, nmax);
    matrix_t tmVsovij = matrix_init(nmax, nmax);
    matrix_t tmVsosi = matrix_init(nmax, nmax);
    matrix_t tmVsosj = matrix_init(nmax, nmax);
    matrix_t tmVtens = matrix_init(nmax, nmax);
    matrix_t Hfi = matrix_init(nmax, nmax);
    matrix_t Nfi = matrix_init(nmax, nmax);

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
            factor_p =  sqrt(4 * basis[i].scale * basis[j].scale / (basis[i].scale + basis[j].scale));

            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_nlp_hamilton(GRnlp, GIVt, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mbetaijCoul.value[i][j] = integral_nlp_hamilton(GRnlp, GIVbetaijcoul, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijCont.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijcont, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaiisov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltajjsov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijsov, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSos.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaiisos, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSos.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltajjsos, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijTens.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijtens, factor_p, &basis[i], &basis[j], args_model, args_dynmc);
            mVcoul.value[i][j] = integral_nlr_hamilton(GRnlr, GIVcoul, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVconf.value[i][j] = integral_nlr_hamilton(GRnlr, GIVconf, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVcont.value[i][j] = integral_nlr_hamilton(GRnlr, GIVcont, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovi.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovi, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovj.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovj, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovij.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovij, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosi.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsosi, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosj.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsosj, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            mVtens.value[i][j] = integral_nlr_hamilton(GRnlr, GIVtens, factor_r, &basis[i], &basis[j], args_model, args_dynmc);
            Nfi.value[i][j] = integral_nlr_overlap(GRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    /* Cholesky decomposition */
    matrix_t mL = matrix_init(nmax, nmax);
    matrix_t mLinv = matrix_init(nmax, nmax);
    matrix_cholesky_decomp(&Nfi, &mL);          /* S = L * L^T */
    matrix_inverse_lowertri(&mL, &mLinv);       /* Linv = L^{-1} */

    /* transform Hamiltonian matrices in orthogonal basis */
    matrix_productT(&mLinv, &mT, &tmT);
    matrix_productT(&mLinv, &mbetaijCoul, &tmbetaijCoul);
    matrix_productT(&mLinv, &mdeltaijCont, &tmdeltaijCont);
    matrix_productT(&mLinv, &mdeltaiiSov, &tmdeltaiiSov);
    matrix_productT(&mLinv, &mdeltajjSov, &tmdeltajjSov);
    matrix_productT(&mLinv, &mdeltaijSov, &tmdeltaijSov);
    matrix_productT(&mLinv, &mdeltaiiSos, &tmdeltaiiSos);
    matrix_productT(&mLinv, &mdeltajjSos, &tmdeltajjSos);
    matrix_productT(&mLinv, &mdeltaijTens, &tmdeltaijTens);
    matrix_productT(&mLinv, &mVcoul, &tmVcoul);
    matrix_productT(&mLinv, &mVconf, &tmVconf);
    matrix_productT(&mLinv, &mVcont, &tmVcont);
    matrix_productT(&mLinv, &mVsovi, &tmVsovi);
    matrix_productT(&mLinv, &mVsovj, &tmVsovj);
    matrix_productT(&mLinv, &mVsovij, &tmVsovij);
    matrix_productT(&mLinv, &mVsosi, &tmVsosi);
    matrix_productT(&mLinv, &mVsosj, &tmVsosj);
    matrix_productT(&mLinv, &mVtens, &tmVtens);

    /* construct Hamiltonian matrix */
    matrix_t temp = matrix_init(nmax, nmax);
    matrix_sum(&tmT, &tmVconf, &Hfi);
    matrix_productT(&tmbetaijCoul, &tmVcoul, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijCont, &tmVcont, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaiiSov, &tmVsovi, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltajjSov, &tmVsovj, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijSov, &tmVsovij, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaiiSos, &tmVsosi, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltajjSos, &tmVsosj, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijTens, &tmVtens, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&mLinv, &Nfi, &temp);

/* In GI model, the final basis should be orthogonal.
 * Therefore this could be general eigen system with orthogonal basis and diagonal Nfi,
 * or just standard eigen system directly. */
#ifdef LAPACKE
    lapack_general(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, v_len);
#else
    eigen_general(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, v_len);
#endif

    free(basis);
    matrix_free(&temp);
    matrix_free(&mT);
    matrix_free(&mbetaijCoul);
    matrix_free(&mdeltaijCont);
    matrix_free(&mdeltaiiSov);
    matrix_free(&mdeltajjSov);
    matrix_free(&mdeltaijSov);
    matrix_free(&mdeltaiiSos);
    matrix_free(&mdeltajjSos);
    matrix_free(&mdeltaijTens);
    matrix_free(&mVcoul);
    matrix_free(&mVconf);
    matrix_free(&mVcont);
    matrix_free(&mVsovi);
    matrix_free(&mVsovj);
    matrix_free(&mVsovij);
    matrix_free(&mVsosi);
    matrix_free(&mVsosj);
    matrix_free(&mVtens);
    matrix_free(&tmT);
    matrix_free(&tmbetaijCoul);
    matrix_free(&tmdeltaijCont);
    matrix_free(&tmdeltaiiSov);
    matrix_free(&tmdeltajjSov);
    matrix_free(&tmdeltaijSov);
    matrix_free(&tmdeltaiiSos);
    matrix_free(&tmdeltajjSos);
    matrix_free(&tmdeltaijTens);
    matrix_free(&tmVcoul);
    matrix_free(&tmVconf);
    matrix_free(&tmVcont);
    matrix_free(&tmVsovi);
    matrix_free(&tmVsovj);
    matrix_free(&tmVsovij);
    matrix_free(&tmVsosi);
    matrix_free(&tmVsosj);
    matrix_free(&tmVtens);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void radius_meson_GEM(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len)
{
    int nmax = input->nmax;
    double rmax = input->rmax;
    double rmin = input->rmin;
    int L = (int)input->L;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
    }

    /* construct matrices */
    matrix_t mR2 = matrix_init(nmax, nmax);
    matrix_t mOver = matrix_init(nmax, nmax);

    /* prepare variables */
    double factor;
    double coef;
    double r2sum;
    double oversum;
    double norm;
    double rms2;
    double fm = 5.06773093854369882649;

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);
            mR2.value[i][j] = integral_nlr_radius(GRnlr, factor, &basis[i], &basis[j]);
            mOver.value[i][j] = integral_nlr_overlap(GRnlr, factor, &basis[i], &basis[j]);
        }
    }

    /* calculate rms radius with orthogonalized coefficients */
    for (int n = 0; n < len; n++) {
        r2sum = 0.0;
        oversum = 0.0;
        norm = 0.0;

        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                coef = vector->value[n][i] * vector->value[n][j];
                r2sum += coef * mR2.value[i][j];
                oversum += coef * mOver.value[i][j];
            }
            norm += vector->value[n][i] * vector->value[n][i];
        }

        rms2 = (norm > 1e-12) ? r2sum / norm : 0.0;
        radius->value[n] = sqrt(rms2) / fm;
    }

    free(basis);
    matrix_free(&mR2);
    matrix_free(&mOver);
}