/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/spectra.h>
#include <gemstore/model/gimodel.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/orbit.h>

#include <gemstore/math/matrix.h>
#include <gemstore/math/cmatrix.h>
#include <gemstore/math/integral.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/ceigen.h>
#include <gemstore/math/soc.h>

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
    double factor;
    double factor_complex;
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
            factor = 1 / sqrt(basis[i].scale + basis[j].scale);
            factor_complex =  sqrt(4 * basis[i].scale * basis[j].scale / (basis[i].scale + basis[j].scale));

            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVt, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mbetaijCoul.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVbetaijcoul, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijCont.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltaijcont, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSov.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltaiisov, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSov.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltajjsov, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijSov.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltaijsov, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaiiSos.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltaiisos, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltajjSos.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltajjsos, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mdeltaijTens.value[i][j] = integral_wfn_hamilton_complex(GRnlp, GIVdeltaijtens, factor_complex, &basis[i], &basis[j], args_model, args_dynmc);
            mVcoul.value[i][j] = integral_wfn_hamilton(GRnlr, GIVcoul, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVconf.value[i][j] = integral_wfn_hamilton(GRnlr, GIVconf, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVcont.value[i][j] = integral_wfn_hamilton(GRnlr, GIVcont, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovi.value[i][j] = integral_wfn_hamilton(GRnlr, GIVsovi, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovj.value[i][j] = integral_wfn_hamilton(GRnlr, GIVsovj, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVsovij.value[i][j] = integral_wfn_hamilton(GRnlr, GIVsovij, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosi.value[i][j] = integral_wfn_hamilton(GRnlr, GIVsosi, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVsosj.value[i][j] = integral_wfn_hamilton(GRnlr, GIVsosj, factor, &basis[i], &basis[j], args_model, args_dynmc);
            mVtens.value[i][j] = integral_wfn_hamilton(GRnlr, GIVtens, factor, &basis[i], &basis[j], args_model, args_dynmc);
            Nfi.value[i][j] = integral_wfn_overlap(GRnlr, factor, &basis[i], &basis[j]);
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

void spectra_meson_CRG(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    double theta = args_input->theta;
    int f1 = args_input->f1, f2 = args_input->f2;
    double S = args_input->S, L = args_input->L, J = args_input->J;
    
    fprintf(stderr, "[CSM] Constructing %dx%d complex Hamiltonian (nmax=%d, theta=%.6f rad)\n", 
            nmax, nmax, nmax, theta);
    
    /* Construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, rmax, rmin);
    }
    
    /* ========================================================================
     * Compute all 17 matrices (9 operator + 8 potential) - COMPLEX
     * ======================================================================== */
    
    cmatrix_t cmT = cmatrix_init(nmax, nmax);
    cmatrix_t cmbetaijCoul = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltaijCont = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltaiiSov = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltajjSov = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltaijSov = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltaiiSos = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltajjSos = cmatrix_init(nmax, nmax);
    cmatrix_t cmdeltaijTens = cmatrix_init(nmax, nmax);
    
    cmatrix_t cmVcoul = cmatrix_init(nmax, nmax);
    cmatrix_t cmVconf = cmatrix_init(nmax, nmax);
    cmatrix_t cmVcont = cmatrix_init(nmax, nmax);
    cmatrix_t cmVsovi = cmatrix_init(nmax, nmax);
    cmatrix_t cmVsovj = cmatrix_init(nmax, nmax);
    cmatrix_t cmVsovij = cmatrix_init(nmax, nmax);
    cmatrix_t cmVsosi = cmatrix_init(nmax, nmax);
    cmatrix_t cmVsosj = cmatrix_init(nmax, nmax);
    cmatrix_t cmVtens = cmatrix_init(nmax, nmax);
    
    cmatrix_t cN_overlap = cmatrix_init(nmax, nmax);
    
    /* Prepare variables */
    double factor, factor_complex;
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
    
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            factor = 1 / sqrt(basis[i].scale + basis[j].scale);
            factor_complex = sqrt(4 * basis[i].scale * basis[j].scale / (basis[i].scale + basis[j].scale));
            
            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
            
             /* DEBUG: Compare first element with GEM (if theta=0) */
             if (i == 0 && j == 0 && theta == 0.0) {
                 complex double csm_t = integral_csm_hamiltonian_complex(GIVt, theta, factor_complex,
                                                                         &basis[i], &basis[j],
                                                                         args_model, args_dynmc);
                 double gem_t = integral_matrix_element_complex(GRnlp, GIVt, factor_complex,
                                                               &basis[i], &basis[j], args_model, args_dynmc);
                 fprintf(stderr, "[CSM-DEBUG] T[0,0]: CSM=%.6e + %.6ei, GEM=%.6e\n", 
                         creal(csm_t), cimag(csm_t), gem_t);
             }
             if (i == 1 && j == 1 && theta == 0.0) {
                 complex double csm_t = integral_csm_hamiltonian_complex(GIVt, theta, factor_complex,
                                                                         &basis[i], &basis[j],
                                                                         args_model, args_dynmc);
                 double gem_t = integral_matrix_element_complex(GRnlp, GIVt, factor_complex,
                                                               &basis[i], &basis[j], args_model, args_dynmc);
                 fprintf(stderr, "[CSM-DEBUG] T[1,1]: CSM=%.6e + %.6ei, GEM=%.6e\n", 
                         creal(csm_t), cimag(csm_t), gem_t);
             }
             if (i == 0 && j == 1 && theta == 0.0) {
                 complex double csm_t = integral_csm_hamiltonian_complex(GIVt, theta, factor_complex,
                                                                         &basis[i], &basis[j],
                                                                         args_model, args_dynmc);
                 double gem_t = integral_matrix_element_complex(GRnlp, GIVt, factor_complex,
                                                               &basis[i], &basis[j], args_model, args_dynmc);
                 fprintf(stderr, "[CSM-DEBUG] T[0,1]: CSM=%.6e + %.6ei, GEM=%.6e\n", 
                         creal(csm_t), cimag(csm_t), gem_t);
             }
            
            /* 9 Operator matrices with complex scaling */
            cmT.value[i][j] = integral_csm_hamiltonian_complex(GIVt, theta, factor_complex,
                                                               &basis[i], &basis[j],
                                                               args_model, args_dynmc);
            cmbetaijCoul.value[i][j] = integral_csm_hamiltonian_complex(GIVbetaijcoul, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltaijCont.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltaijcont, theta, factor_complex,
                                                                         &basis[i], &basis[j],
                                                                         args_model, args_dynmc);
            cmdeltaiiSov.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltaiisov, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltajjSov.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltajjsov, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltaijSov.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltaijsov, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltaiiSos.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltaiisos, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltajjSos.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltajjsos, theta, factor_complex,
                                                                        &basis[i], &basis[j],
                                                                        args_model, args_dynmc);
            cmdeltaijTens.value[i][j] = integral_csm_hamiltonian_complex(GIVdeltaijtens, theta, factor_complex,
                                                                         &basis[i], &basis[j],
                                                                         args_model, args_dynmc);
             
            /* 8 Potential matrices with complex scaling (use factor_complex like GEM's complex integrals) */
             cmVcoul.value[i][j] = integral_csm_potential_complex(GIVcoul, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVconf.value[i][j] = integral_csm_potential_complex(GIVconf, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVcont.value[i][j] = integral_csm_potential_complex(GIVcont, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVsovi.value[i][j] = integral_csm_potential_complex(GIVsovi, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVsovj.value[i][j] = integral_csm_potential_complex(GIVsovj, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVsovij.value[i][j] = integral_csm_potential_complex(GIVsovij, theta, factor,
                                                                    &basis[i], &basis[j],
                                                                    args_model, args_dynmc);
             cmVsosi.value[i][j] = integral_csm_potential_complex(GIVsosi, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVsosj.value[i][j] = integral_csm_potential_complex(GIVsosj, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
             cmVtens.value[i][j] = integral_csm_potential_complex(GIVtens, theta, factor,
                                                                   &basis[i], &basis[j],
                                                                   args_model, args_dynmc);
            
            /* Overlap matrix - COMPLEX (use factor_complex to match Hamiltonian integrals) */
            cN_overlap.value[i][j] = integral_csm_overlap_complex(theta, factor_complex, &basis[i], &basis[j]);
        }
    }
    
    /* DEBUG: Print first few overlap matrix elements */
    fprintf(stderr, "[CSM-OVERLAP] Matrix N[i][j] (first 3x3):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            fprintf(stderr, "[%.6e, %.6e]  ", creal(cN_overlap.value[i][j]), cimag(cN_overlap.value[i][j]));
        }
        fprintf(stderr, "\n");
    }
    
    /* DEBUG: Check if matrix is Hermitian */
    int hermitian = 1;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double complex diff = cN_overlap.value[i][j] - conj(cN_overlap.value[j][i]);
            if (cabs(diff) > 1e-10) {
                hermitian = 0;
                fprintf(stderr, "[CSM-ERROR] N is not Hermitian at [%d,%d]: diff=%.6e\n", i, j, cabs(diff));
            }
        }
    }
    if (hermitian) fprintf(stderr, "[CSM-OK] Overlap matrix is Hermitian\n");
    
    /* ========================================================================
     * Complex Cholesky decomposition and transformation to orthogonal basis
     * ======================================================================== */
    
    cmatrix_t cmL = cmatrix_init(nmax, nmax);
    cmatrix_t cmLinv = cmatrix_init(nmax, nmax);
    cmatrix_cholesky_decomp(&cN_overlap, &cmL);
    cmatrix_inverse_lowertri(&cmL, &cmLinv);
    
    /* Transform all matrices to orthogonal basis using complex operations */
    cmatrix_t ctmT = cmatrix_init(nmax, nmax);
    cmatrix_t ctmbetaijCoul = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltaijCont = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltaiiSov = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltajjSov = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltaijSov = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltaiiSos = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltajjSos = cmatrix_init(nmax, nmax);
    cmatrix_t ctmdeltaijTens = cmatrix_init(nmax, nmax);
    
    cmatrix_t ctmVcoul = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVconf = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVcont = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVsovi = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVsovj = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVsovij = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVsosi = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVsosj = cmatrix_init(nmax, nmax);
    cmatrix_t ctmVtens = cmatrix_init(nmax, nmax);
    
    cmatrix_productT(&cmLinv, &cmT, &ctmT);
    cmatrix_productT(&cmLinv, &cmbetaijCoul, &ctmbetaijCoul);
    cmatrix_productT(&cmLinv, &cmdeltaijCont, &ctmdeltaijCont);
    cmatrix_productT(&cmLinv, &cmdeltaiiSov, &ctmdeltaiiSov);
    cmatrix_productT(&cmLinv, &cmdeltajjSov, &ctmdeltajjSov);
    cmatrix_productT(&cmLinv, &cmdeltaijSov, &ctmdeltaijSov);
    cmatrix_productT(&cmLinv, &cmdeltaiiSos, &ctmdeltaiiSos);
    cmatrix_productT(&cmLinv, &cmdeltajjSos, &ctmdeltajjSos);
    cmatrix_productT(&cmLinv, &cmdeltaijTens, &ctmdeltaijTens);
    
    cmatrix_productT(&cmLinv, &cmVcoul, &ctmVcoul);
    cmatrix_productT(&cmLinv, &cmVconf, &ctmVconf);
    cmatrix_productT(&cmLinv, &cmVcont, &ctmVcont);
    cmatrix_productT(&cmLinv, &cmVsovi, &ctmVsovi);
    cmatrix_productT(&cmLinv, &cmVsovj, &ctmVsovj);
    cmatrix_productT(&cmLinv, &cmVsovij, &ctmVsovij);
    cmatrix_productT(&cmLinv, &cmVsosi, &ctmVsosi);
    cmatrix_productT(&cmLinv, &cmVsosj, &ctmVsosj);
    cmatrix_productT(&cmLinv, &cmVtens, &ctmVtens);
    
    /* ========================================================================
     * Construct Hamiltonian from GI operators (same as GEM pattern)
     * H = T + Vconf + β·Vcoul + δ_cont·Vcont + δ_i_sov·Vsovi + δ_j_sov·Vsovj + 
     *     δ_ij_sov·Vsovij + δ_i_sos·Vsosi + δ_j_sos·Vsosj + tensor·Vtens
     * ======================================================================== */
    
    /* DEBUG: Print first element of transformed matrix */
    fprintf(stderr, "[CSM-TRANSFORMED] ctmT[0,0] = %.6e + %.6ei\n", creal(ctmT.value[0][0]), cimag(ctmT.value[0][0]));
    fprintf(stderr, "[CSM-TRANSFORMED] ctmT[0,1] = %.6e + %.6ei\n", creal(ctmT.value[0][1]), cimag(ctmT.value[0][1]));
    fprintf(stderr, "[CSM-TRANSFORMED] ctmT[1,1] = %.6e + %.6ei\n", creal(ctmT.value[1][1]), cimag(ctmT.value[1][1]));
    
    cmatrix_t cH_final = cmatrix_init(nmax, nmax);
    cmatrix_t ctemp = cmatrix_init(nmax, nmax);
    
    cmatrix_sum(&ctmT, &ctmVconf, &cH_final);
    cmatrix_productT(&ctmbetaijCoul, &ctmVcoul, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltaijCont, &ctmVcont, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltaiiSov, &ctmVsovi, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltajjSov, &ctmVsovj, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltaijSov, &ctmVsovij, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltaiiSos, &ctmVsosi, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltajjSos, &ctmVsosj, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    cmatrix_productT(&ctmdeltaijTens, &ctmVtens, &ctemp);
    cmatrix_sum(&cH_final, &ctemp, &cH_final);
    
    /* Compute transformed overlap for eigenvalue problem: L^{-1} * N * (L^{-1})^† */
    cmatrix_t cN_final = cmatrix_init(nmax, nmax);
    cmatrix_productT(&cmLinv, &cN_overlap, &cN_final);
    
    /* ========================================================================
     * Solve complex generalized eigenvalue problem: H x = λ S x
     * Returns complex eigenvalues (poles in complex plane for resonances)
     * ======================================================================== */
    
    carray_t c_eigenvalues = carray_init(nmax);
    cmatrix_t c_eigenvectors = cmatrix_init((v_out == NULL) ? 0 : v_len, nmax);
    
    /* DEBUG: Check if transformed overlap matrix is identity */
    fprintf(stderr, "[CSM-NFINAL] cN_final[i][j] (first 3x3):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            fprintf(stderr, "[%.6e, %.6e]  ", creal(cN_final.value[i][j]), cimag(cN_final.value[i][j]));
        }
        fprintf(stderr, "\n");
    }
    
    /* DEBUG: Show final Hamiltonian diagonal elements */
    fprintf(stderr, "[CSM-HFINAL] cH_final[i][i] (first 5):\n");
    for (int i = 0; i < 5; i++) {
        fprintf(stderr, "  [%d,%-d] = %.6e + %.6ei\n", i, i, creal(cH_final.value[i][i]), cimag(cH_final.value[i][i]));
    }
    
    #ifdef LAPACKE
        lapack_general_complex(cH_final.value, cN_final.value, nmax, c_eigenvalues.value,
                              (v_out == NULL) ? NULL : c_eigenvectors.value, 
                              (v_out == NULL) ? 0 : v_len);
    #else
        eigen_general_complex(cH_final.value, cN_final.value, nmax, c_eigenvalues.value,
                             (v_out == NULL) ? NULL : c_eigenvectors.value, 
                             (v_out == NULL) ? 0 : v_len);
    #endif
    
    /* Extract real parts of complex eigenvalues for output
     * Full complex eigenvalues contain resonance pole information:
     * Re(λ) = energy position, Im(λ) = -width/2 (width of resonance)
     */
    for (int i = 0; i < nmax; i++) {
        e_out->value[i] = creal(c_eigenvalues.value[i]);
    }
    
    /* If eigenvectors requested, convert complex to real (take real part) */
    if (v_out != NULL) {
        for (int i = 0; i < v_len && i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                v_out->value[i][j] = creal(c_eigenvectors.value[i][j]);
            }
        }
    }
    
    fprintf(stderr, "[CSM] Complex Hamiltonian solved successfully (theta=%.6f rad)\n", theta);
    fprintf(stderr, "[CSM] Eigenvalues are complex (resonance poles):\n");
    for (int i = 0; i < ((nmax < 5) ? nmax : 5); i++) {
        fprintf(stderr, "      λ_%d = %.6f + %.6fi (width = %.6f)\n", 
                i+1, creal(c_eigenvalues.value[i]), cimag(c_eigenvalues.value[i]),
                -2.0 * cimag(c_eigenvalues.value[i]));
    }
    if (nmax > 5) fprintf(stderr, "      ... and %d more\n", nmax - 5);
    
    /* DEBUG: Print final Hamiltonian diagonal */
    fprintf(stderr, "[CSM] Final Hamiltonian diagonal elements:\n");
    for (int i = 0; i < ((nmax < 5) ? nmax : 5); i++) {
        fprintf(stderr, "      H[%d,%d] = %.6f + %.6fi\n", i, i, creal(cH_final.value[i][i]), cimag(cH_final.value[i][i]));
    }
    
    /* Cleanup */
    free(basis);
    cmatrix_free(&cmT);
    cmatrix_free(&cmbetaijCoul);
    cmatrix_free(&cmdeltaijCont);
    cmatrix_free(&cmdeltaiiSov);
    cmatrix_free(&cmdeltajjSov);
    cmatrix_free(&cmdeltaijSov);
    cmatrix_free(&cmdeltaiiSos);
    cmatrix_free(&cmdeltajjSos);
    cmatrix_free(&cmdeltaijTens);
    cmatrix_free(&cmVcoul);
    cmatrix_free(&cmVconf);
    cmatrix_free(&cmVcont);
    cmatrix_free(&cmVsovi);
    cmatrix_free(&cmVsovj);
    cmatrix_free(&cmVsovij);
    cmatrix_free(&cmVsosi);
    cmatrix_free(&cmVsosj);
    cmatrix_free(&cmVtens);
    cmatrix_free(&cN_overlap);
    cmatrix_free(&cmL);
    cmatrix_free(&cmLinv);
    cmatrix_free(&ctmT);
    cmatrix_free(&ctmbetaijCoul);
    cmatrix_free(&ctmdeltaijCont);
    cmatrix_free(&ctmdeltaiiSov);
    cmatrix_free(&ctmdeltajjSov);
    cmatrix_free(&ctmdeltaijSov);
    cmatrix_free(&ctmdeltaiiSos);
    cmatrix_free(&ctmdeltajjSos);
    cmatrix_free(&ctmdeltaijTens);
    cmatrix_free(&ctmVcoul);
    cmatrix_free(&ctmVconf);
    cmatrix_free(&ctmVcont);
    cmatrix_free(&ctmVsovi);
    cmatrix_free(&ctmVsovj);
    cmatrix_free(&ctmVsovij);
    cmatrix_free(&ctmVsosi);
    cmatrix_free(&ctmVsosj);
    cmatrix_free(&ctmVtens);
    cmatrix_free(&cH_final);
    cmatrix_free(&ctemp);
    cmatrix_free(&cN_final);
    carray_free(&c_eigenvalues);
    cmatrix_free(&c_eigenvectors);
 }