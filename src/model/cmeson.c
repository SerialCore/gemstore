/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/cmeson.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/model/nrmodel.h>

#include <gemstore/basis/orbit.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>
#include <gemstore/param/argset.h>

#include <gemstore/input.h>

#include <stdio.h>
#include <stdlib.h>

static void add_reduced(matrix_t *H, const matrix_t *vt, const matrix_t *M, matrix_t *tmp)
{
    matrix_productT(vt, M, tmp);
    matrix_sum(H, tmp, H);
}

static void add_sandwiched(matrix_t *H, const matrix_t *vt, const matrix_t *V, const matrix_t *P,
    matrix_t *tV, matrix_t *tP, matrix_t *tmp)
{
    matrix_productT(vt, V, tV);
    matrix_productT(vt, P, tP);
	matrix_productT(tP, tV, tmp);
    matrix_sum(H, tmp, H);
}

/* RMS from (c^T r^2 c) / (c^T N c) in the original basis. Allocates *rms_out when both pointers are set. */
static void meson_radius_RMS(const argsInput_t *input, const argsOrbit_t *basis, int nmax,
    const matrix_t *Nfi, matrix_t *v_out, array_t *rms_out)
{
    if (rms_out == NULL || v_out == NULL || Nfi == NULL) {
        return;
    }

    matrix_t mR2 = matrix_init(nmax, nmax);
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            if (input->basis == BASIS_GEM) {
                double factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);
                mR2.value[i][j] = integral_nlr_radius(GRnlr, factor, &basis[i], &basis[j]);
            }
            else {
                double factor = sqrt(2.0 / (basis[i].scale * basis[i].scale + basis[j].scale * basis[j].scale));
                mR2.value[i][j] = integral_nlr_radius(SRnlr, factor, &basis[i], &basis[j]);
            }
        }
    }

    *rms_out = array_init(v_out->row);
    for (int n = 0; n < v_out->row; n++) {
        double nrm = matrix_expect(v_out, n, Nfi);
        double r2 = matrix_expect(v_out, n, &mR2);
        if (nrm > 1e-16 && r2 > 0.0) {
            rms_out->value[n] = sqrt(r2 / nrm) / GEMSTORE_FM;
        }
        else {
            rms_out->value[n] = 0.0;
        }
    }
    matrix_free(&mR2);
}

void meson_gimodel_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
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
    matrix_t Nfi = matrix_init(nmax, nmax);
    matrix_t Hfi = matrix_init(nmax, nmax);

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

    gi_pot_ctx_t ctx = { args_model, args_dynmc };

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
                
            mT.value[i][j] = integral_nlp_hamilton(GRnlp, GIVt, factor_p, &basis[i], &basis[j], &ctx);
            mbetaijCoul.value[i][j] = integral_nlp_hamilton(GRnlp, GIVbetaijcoul, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijCont.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijcont, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaiiSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaiisov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltajjSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltajjsov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijSov.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijsov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaiiSos.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaiisos, factor_p, &basis[i], &basis[j], &ctx);
            mdeltajjSos.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltajjsos, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijTens.value[i][j] = integral_nlp_hamilton(GRnlp, GIVdeltaijtens, factor_p, &basis[i], &basis[j], &ctx);
            mVcoul.value[i][j] = integral_nlr_hamilton(GRnlr, GIVcoul, factor_r, &basis[i], &basis[j], &ctx);
            mVconf.value[i][j] = integral_nlr_hamilton(GRnlr, GIVconf, factor_r, &basis[i], &basis[j], &ctx);
            mVcont.value[i][j] = integral_nlr_hamilton(GRnlr, GIVcont, factor_r, &basis[i], &basis[j], &ctx);
            mVsovi.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovi, factor_r, &basis[i], &basis[j], &ctx);
            mVsovj.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovj, factor_r, &basis[i], &basis[j], &ctx);
            mVsovij.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsovij, factor_r, &basis[i], &basis[j], &ctx);
            mVsosi.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsosi, factor_r, &basis[i], &basis[j], &ctx);
            mVsosj.value[i][j] = integral_nlr_hamilton(GRnlr, GIVsosj, factor_r, &basis[i], &basis[j], &ctx);
            mVtens.value[i][j] = integral_nlr_hamilton(GRnlr, GIVtens, factor_r, &basis[i], &basis[j], &ctx);
            Nfi.value[i][j] = integral_nlr_overlap(GRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    matrix_symmetrize(&mT);
    matrix_symmetrize(&mbetaijCoul);
    matrix_symmetrize(&mdeltaijCont);
    matrix_symmetrize(&mdeltaiiSov);
    matrix_symmetrize(&mdeltajjSov);
    matrix_symmetrize(&mdeltaijSov);
    matrix_symmetrize(&mdeltaiiSos);
    matrix_symmetrize(&mdeltajjSos);
    matrix_symmetrize(&mdeltaijTens);
    matrix_symmetrize(&mVcoul);
    matrix_symmetrize(&mVconf);
    matrix_symmetrize(&mVcont);
    matrix_symmetrize(&mVsovi);
    matrix_symmetrize(&mVsovj);
    matrix_symmetrize(&mVsovij);
    matrix_symmetrize(&mVsosi);
    matrix_symmetrize(&mVsosj);
    matrix_symmetrize(&mVtens);
    matrix_symmetrize(&Nfi);

    /* construct orthogonal basis */
    matrix_t vt;
    int nkeep = matrix_orthogonal(&Nfi, 1e-8, &vt);
    if (nkeep <= 0) {
        fprintf(stderr, "Error: meson overlap N has no positive eigenvalues.\n");
        exit(1);
    }

    /* construct Hamiltonian matrix */
    matrix_t temp = matrix_init(nkeep, nkeep);
    matrix_t tV = matrix_init(nkeep, nkeep);
    matrix_t tP = matrix_init(nkeep, nkeep);
    
    add_reduced(&Hfi, &vt, &mT, &temp);
    add_reduced(&Hfi, &vt, &mVconf, &temp);
    add_sandwiched(&Hfi, &vt, &mVcoul, &mbetaijCoul, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVcont, &mdeltaijCont, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVsovi, &mdeltaiiSov, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVsovj, &mdeltajjSov, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVsovij, &mdeltaijSov, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVsosi, &mdeltaiiSos, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVsosj, &mdeltajjSos, &tV, &tP, &temp);
    add_sandwiched(&Hfi, &vt, &mVtens, &mdeltaijTens, &tV, &tP, &temp);
    
    if (e_out == NULL || e_out->value == NULL || e_out->len < nkeep) {
        fprintf(stderr, "Error: meson eigenvalue buffer cannot hold rank %d.\n", nkeep);
        exit(1);
    }

    matrix_t ut = {0};
    if (v_out != NULL) {
        ut = matrix_init(nkeep, nkeep);
    }

#ifdef LAPACKE
    lapack_standard(Hfi.value, nkeep, e_out->value, (v_out != NULL) ? ut.value : NULL, (v_out != NULL) ? nkeep : 0);
#else
    eigen_standard(Hfi.value, nkeep, e_out->value, (v_out != NULL) ? ut.value : NULL, (v_out != NULL) ? nkeep : 0);
#endif

    e_out->len = nkeep;

    if (v_out != NULL) {
        for (int i = nkeep; i < v_out->row; i++) {
            free(v_out->value[i]);
        }
        v_out->row = nkeep;
        matrix_product(&ut, &vt, v_out);
        matrix_free(&ut);
    }
    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    meson_radius_RMS(args_input, basis, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&temp);
    matrix_free(&tV);
    matrix_free(&tP);
    matrix_free(&vt);
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
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void meson_gimodel_SHO(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double beta = args_input->beta;
    int f1 = args_input->f1, f2 = args_input->f2;
    double S = args_input->S, L = args_input->L, J = args_input->J;

    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i;
        basis[i].l = L;
        basis[i].scale = beta;
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

    gi_pot_ctx_t ctx = { args_model, args_dynmc };

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor_r = sqrt(2 / (basis[i].scale * basis[i].scale + basis[j].scale * basis[j].scale));
            factor_p = sqrt(2 * basis[i].scale * basis[i].scale * basis[j].scale * basis[j].scale 
                / (basis[i].scale * basis[i].scale + basis[j].scale * basis[j].scale));

            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_nlp_hamilton(SRnlp, GIVt, factor_p, &basis[i], &basis[j], &ctx);
            mbetaijCoul.value[i][j] = integral_nlp_hamilton(SRnlp, GIVbetaijcoul, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijCont.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltaijcont, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaiiSov.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltaiisov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltajjSov.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltajjsov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijSov.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltaijsov, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaiiSos.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltaiisos, factor_p, &basis[i], &basis[j], &ctx);
            mdeltajjSos.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltajjsos, factor_p, &basis[i], &basis[j], &ctx);
            mdeltaijTens.value[i][j] = integral_nlp_hamilton(SRnlp, GIVdeltaijtens, factor_p, &basis[i], &basis[j], &ctx);
            mVcoul.value[i][j] = integral_nlr_hamilton(SRnlr, GIVcoul, factor_r, &basis[i], &basis[j], &ctx);
            mVconf.value[i][j] = integral_nlr_hamilton(SRnlr, GIVconf, factor_r, &basis[i], &basis[j], &ctx);
            mVcont.value[i][j] = integral_nlr_hamilton(SRnlr, GIVcont, factor_r, &basis[i], &basis[j], &ctx);
            mVsovi.value[i][j] = integral_nlr_hamilton(SRnlr, GIVsovi, factor_r, &basis[i], &basis[j], &ctx);
            mVsovj.value[i][j] = integral_nlr_hamilton(SRnlr, GIVsovj, factor_r, &basis[i], &basis[j], &ctx);
            mVsovij.value[i][j] = integral_nlr_hamilton(SRnlr, GIVsovij, factor_r, &basis[i], &basis[j], &ctx);
            mVsosi.value[i][j] = integral_nlr_hamilton(SRnlr, GIVsosi, factor_r, &basis[i], &basis[j], &ctx);
            mVsosj.value[i][j] = integral_nlr_hamilton(SRnlr, GIVsosj, factor_r, &basis[i], &basis[j], &ctx);
            mVtens.value[i][j] = integral_nlr_hamilton(SRnlr, GIVtens, factor_r, &basis[i], &basis[j], &ctx);
            Nfi.value[i][j] = integral_nlr_overlap(SRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    matrix_symmetrize(&mT);
    matrix_symmetrize(&mbetaijCoul);
    matrix_symmetrize(&mdeltaijCont);
    matrix_symmetrize(&mdeltaiiSov);
    matrix_symmetrize(&mdeltajjSov);
    matrix_symmetrize(&mdeltaijSov);
    matrix_symmetrize(&mdeltaiiSos);
    matrix_symmetrize(&mdeltajjSos);
    matrix_symmetrize(&mdeltaijTens);
    matrix_symmetrize(&mVcoul);
    matrix_symmetrize(&mVconf);
    matrix_symmetrize(&mVcont);
    matrix_symmetrize(&mVsovi);
    matrix_symmetrize(&mVsovj);
    matrix_symmetrize(&mVsovij);
    matrix_symmetrize(&mVsosi);
    matrix_symmetrize(&mVsosj);
    matrix_symmetrize(&mVtens);
    matrix_symmetrize(&Nfi);

    /* construct Hamiltonian matrix */
    matrix_t temp = matrix_init(nmax, nmax);
    matrix_sum(&mT, &mVconf, &Hfi);
    add_reduced(&Hfi, &mbetaijCoul, &mVcoul, &temp);
    add_reduced(&Hfi, &mdeltaijCont, &mVcont, &temp);
    add_reduced(&Hfi, &mdeltaiiSov, &mVsovi, &temp);
    add_reduced(&Hfi, &mdeltajjSov, &mVsovj, &temp);
    add_reduced(&Hfi, &mdeltaijSov, &mVsovij, &temp);
    add_reduced(&Hfi, &mdeltaiiSos, &mVsosi, &temp);
    add_reduced(&Hfi, &mdeltajjSos, &mVsosj, &temp);
    add_reduced(&Hfi, &mdeltaijTens, &mVtens, &temp);

#ifdef LAPACKE
    lapack_standard(Hfi.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, nmax);
#else
    eigen_standard(Hfi.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, nmax);
#endif

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    meson_radius_RMS(args_input, basis, nmax, &Nfi, v_out, rms_out);

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
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void meson_nrmodel_GEM(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    int f1 = args_input->f1, f2 = args_input->f2;
    double L = args_input->L;

    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
    }

    matrix_t mT = matrix_init(nmax, nmax);
    matrix_t mVcoul = matrix_init(nmax, nmax);
    matrix_t mVconf = matrix_init(nmax, nmax);
    matrix_t mVcont = matrix_init(nmax, nmax);
    matrix_t mVsov = matrix_init(nmax, nmax);
    matrix_t mVsos = matrix_init(nmax, nmax);
    matrix_t mVtens = matrix_init(nmax, nmax);
    matrix_t Hfi = matrix_init(nmax, nmax);
    matrix_t Nfi = matrix_init(nmax, nmax);

    double s1 = 0.5, s2 = 0.5;
    double S = args_input->S, J = args_input->J;
    args_dynmc->mi = getmq(f1, args_model);
    args_dynmc->mj = getmq(f2, args_model);
    args_dynmc->Cij = -4.0 / 3.0;
    args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);

    nr_pot_ctx_t ctx = { args_model, args_dynmc };

    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double factor_r = 1 / sqrt(basis[i].scale + basis[j].scale);
            double factor_p = sqrt(4 * basis[i].scale * basis[j].scale / (basis[i].scale + basis[j].scale));

            mT.value[i][j] = integral_nlp_hamilton(GRnlp, NRVt, factor_p, &basis[i], &basis[j], &ctx);
            mVcoul.value[i][j] = integral_nlr_hamilton(GRnlr, NRVcoul, factor_r, &basis[i], &basis[j], &ctx);
            mVconf.value[i][j] = integral_nlr_hamilton(GRnlr, NRVconf, factor_r, &basis[i], &basis[j], &ctx);
            mVcont.value[i][j] = integral_nlr_hamilton(GRnlr, NRVcont, factor_r, &basis[i], &basis[j], &ctx);
            mVsov.value[i][j] = integral_nlr_hamilton(GRnlr, NRVsov, factor_r, &basis[i], &basis[j], &ctx);
            mVsos.value[i][j] = integral_nlr_hamilton(GRnlr, NRVsos, factor_r, &basis[i], &basis[j], &ctx);
            mVtens.value[i][j] = integral_nlr_hamilton(GRnlr, NRVtens, factor_r, &basis[i], &basis[j], &ctx);
            Nfi.value[i][j] = integral_nlr_overlap(GRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    matrix_symmetrize(&mT);
    matrix_symmetrize(&mVcoul);
    matrix_symmetrize(&mVconf);
    matrix_symmetrize(&mVcont);
    matrix_symmetrize(&mVsov);
    matrix_symmetrize(&mVsos);
    matrix_symmetrize(&mVtens);
    matrix_symmetrize(&Nfi);

    matrix_sum(&mT, &mVcoul, &Hfi);
    matrix_sum(&Hfi, &mVconf, &Hfi);
    matrix_sum(&Hfi, &mVcont, &Hfi);
    matrix_sum(&Hfi, &mVsov, &Hfi);
    matrix_sum(&Hfi, &mVsos, &Hfi);
    matrix_sum(&Hfi, &mVtens, &Hfi);

#ifdef LAPACKE
    lapack_general(Hfi.value, Nfi.value, nmax, e_out->value, (v_out == NULL) ? NULL : v_out->value, nmax);
#else
    eigen_general(Hfi.value, Nfi.value, nmax, e_out->value, (v_out == NULL) ? NULL : v_out->value, nmax);
#endif

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    meson_radius_RMS(args_input, basis, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&mT);
    matrix_free(&mVcoul);
    matrix_free(&mVconf);
    matrix_free(&mVcont);
    matrix_free(&mVsov);
    matrix_free(&mVsos);
    matrix_free(&mVtens);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void meson_nrmodel_SHO(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double beta = args_input->beta;
    int f1 = args_input->f1, f2 = args_input->f2;
    double L = args_input->L;

    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i;
        basis[i].l = L;
        basis[i].scale = beta;
    }

    matrix_t mT = matrix_init(nmax, nmax);
    matrix_t mVcoul = matrix_init(nmax, nmax);
    matrix_t mVconf = matrix_init(nmax, nmax);
    matrix_t mVcont = matrix_init(nmax, nmax);
    matrix_t mVsov = matrix_init(nmax, nmax);
    matrix_t mVsos = matrix_init(nmax, nmax);
    matrix_t mVtens = matrix_init(nmax, nmax);
    matrix_t Hfi = matrix_init(nmax, nmax);
    matrix_t Nfi = matrix_init(nmax, nmax);

    double s1 = 0.5, s2 = 0.5;
    double S = args_input->S, J = args_input->J;
    args_dynmc->mi = getmq(f1, args_model);
    args_dynmc->mj = getmq(f2, args_model);
    args_dynmc->Cij = -4.0 / 3.0;
    args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);

    nr_pot_ctx_t ctx = { args_model, args_dynmc };

    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double factor_r = sqrt(2 / (basis[i].scale * basis[i].scale + basis[j].scale * basis[j].scale));
            double factor_p = sqrt(2 * basis[i].scale * basis[i].scale * basis[j].scale * basis[j].scale
                / (basis[i].scale * basis[i].scale + basis[j].scale * basis[j].scale));

            mT.value[i][j] = integral_nlp_hamilton(SRnlp, NRVt, factor_p, &basis[i], &basis[j], &ctx);
            mVcoul.value[i][j] = integral_nlr_hamilton(SRnlr, NRVcoul, factor_r, &basis[i], &basis[j], &ctx);
            mVconf.value[i][j] = integral_nlr_hamilton(SRnlr, NRVconf, factor_r, &basis[i], &basis[j], &ctx);
            mVcont.value[i][j] = integral_nlr_hamilton(SRnlr, NRVcont, factor_r, &basis[i], &basis[j], &ctx);
            mVsov.value[i][j] = integral_nlr_hamilton(SRnlr, NRVsov, factor_r, &basis[i], &basis[j], &ctx);
            mVsos.value[i][j] = integral_nlr_hamilton(SRnlr, NRVsos, factor_r, &basis[i], &basis[j], &ctx);
            mVtens.value[i][j] = integral_nlr_hamilton(SRnlr, NRVtens, factor_r, &basis[i], &basis[j], &ctx);
            Nfi.value[i][j] = integral_nlr_overlap(SRnlr, factor_r, &basis[i], &basis[j]);
        }
    }

    matrix_symmetrize(&mT);
    matrix_symmetrize(&mVcoul);
    matrix_symmetrize(&mVconf);
    matrix_symmetrize(&mVcont);
    matrix_symmetrize(&mVsov);
    matrix_symmetrize(&mVsos);
    matrix_symmetrize(&mVtens);
    matrix_symmetrize(&Nfi);

    matrix_sum(&mT, &mVcoul, &Hfi);
    matrix_sum(&Hfi, &mVconf, &Hfi);
    matrix_sum(&Hfi, &mVcont, &Hfi);
    matrix_sum(&Hfi, &mVsov, &Hfi);
    matrix_sum(&Hfi, &mVsos, &Hfi);
    matrix_sum(&Hfi, &mVtens, &Hfi);

#ifdef LAPACKE
    lapack_standard(Hfi.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, nmax);
#else
    eigen_standard(Hfi.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, nmax);
#endif

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    meson_radius_RMS(args_input, basis, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&mT);
    matrix_free(&mVcoul);
    matrix_free(&mVconf);
    matrix_free(&mVcont);
    matrix_free(&mVsov);
    matrix_free(&mVsos);
    matrix_free(&mVtens);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}