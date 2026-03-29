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
#include <gemstore/math/integral.h>
#include <gemstore/math/eigen.h>
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

void spectra_meson_GI(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
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
    argsOrbit_t args_bra;
    argsOrbit_t args_ket;
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
            args_bra = basis[i];
            args_ket = basis[j];

            factor = 1 / sqrt(args_bra.scale + args_ket.scale);
            factor_complex =  sqrt(4 * args_bra.scale * args_ket.scale / (args_bra.scale + args_ket.scale));

            args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_matrix_element_complex(GRnlp, GIVt, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mbetaijCoul.value[i][j] = integral_matrix_element_complex(GRnlp, GIVbetaijcoul, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltaijCont.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltaijcont, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltaiiSov.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltaiisov, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltajjSov.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltajjsov, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltaijSov.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltaijsov, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltaiiSos.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltaiisos, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltajjSos.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltajjsos, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mdeltaijTens.value[i][j] = integral_matrix_element_complex(GRnlp, GIVdeltaijtens, factor_complex, &args_bra, &args_ket, args_model, args_dynmc);
            mVcoul.value[i][j] = integral_matrix_element(GRnlr, GIVcoul, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVconf.value[i][j] = integral_matrix_element(GRnlr, GIVconf, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVcont.value[i][j] = integral_matrix_element(GRnlr, GIVcont, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVsovi.value[i][j] = integral_matrix_element(GRnlr, GIVsovi, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVsovj.value[i][j] = integral_matrix_element(GRnlr, GIVsovj, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVsovij.value[i][j] = integral_matrix_element(GRnlr, GIVsovij, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVsosi.value[i][j] = integral_matrix_element(GRnlr, GIVsosi, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVsosj.value[i][j] = integral_matrix_element(GRnlr, GIVsosj, factor, &args_bra, &args_ket, args_model, args_dynmc);
            mVtens.value[i][j] = integral_matrix_element(GRnlr, GIVtens, factor, &args_bra, &args_ket, args_model, args_dynmc);
            Nfi.value[i][j] = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_ket);
        }
    }

    /* prepare a random symmetric matrix */
    matrix_t rand = matrix_random(nmax, nmax);
    matrix_t temp = matrix_init(nmax, nmax);
    matrix_transpose(&rand, &temp);
    matrix_sum(&rand, &temp, &rand);
    
    matrix_t vt = matrix_init(nmax, nmax);

#ifdef LAPACKE
    lapack_general(rand.value, Nfi.value, nmax, e_out->value, vt.value, nmax);
#else
    eigen_general(rand.value, Nfi.value, nmax, e_out->value, vt.value, nmax);
#endif

    /* construct new orthogonal basis */
    matrix_productT(&vt, &Nfi, &temp);
    for (int k = 0; k < nmax; k++) {
        double norm = sqrt(fabs(temp.value[k][k]));  /* sqrt(Diagonal[k]) */
        if (norm > 1e-10) {
            for (int i = 0; i < nmax; i++) {
                vt.value[k][i] /= norm;             /* vt[i][k] /= norm for col major, vt[k][i] /= norm for row major */
            }
        }
        else {
            printf("Warning: singular vector %d, norm=%.2e\n", k, norm);
        }
    }
    
    /* transform Hamiltonian matrices in new basis */
    matrix_productT(&vt, &mT, &tmT);
    matrix_productT(&vt, &mbetaijCoul, &tmbetaijCoul);
    matrix_productT(&vt, &mdeltaijCont, &tmdeltaijCont);
    matrix_productT(&vt, &mdeltaiiSov, &tmdeltaiiSov);
    matrix_productT(&vt, &mdeltajjSov, &tmdeltajjSov);
    matrix_productT(&vt, &mdeltaijSov, &tmdeltaijSov);
    matrix_productT(&vt, &mdeltaiiSos, &tmdeltaiiSos);
    matrix_productT(&vt, &mdeltajjSos, &tmdeltajjSos);
    matrix_productT(&vt, &mdeltaijTens, &tmdeltaijTens);
    matrix_productT(&vt, &mVcoul, &tmVcoul);
    matrix_productT(&vt, &mVconf, &tmVconf);
    matrix_productT(&vt, &mVcont, &tmVcont);
    matrix_productT(&vt, &mVsovi, &tmVsovi);
    matrix_productT(&vt, &mVsovj, &tmVsovj);
    matrix_productT(&vt, &mVsovij, &tmVsovij);
    matrix_productT(&vt, &mVsosi, &tmVsosi);
    matrix_productT(&vt, &mVsosj, &tmVsosj);
    matrix_productT(&vt, &mVtens, &tmVtens);

    /* construct Hamiltonian matrix */
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
    matrix_productT(&vt, &Nfi, &temp);

    /* final eigen system */
    matrix_t ut = matrix_init(nmax, nmax);

#ifdef LAPACKE
    lapack_general(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : ut.value, v_len);
#else
    eigen_general(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : ut.value, v_len);
#endif

    if (v_out != NULL) {
        matrix_product(&ut, &vt, v_out);
    }

    matrix_free(&rand);
    matrix_free(&temp);
    matrix_free(&vt);
    matrix_free(&ut);

    free(basis);
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