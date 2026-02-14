/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/spectra.h>
#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/integral.h>
#include <gemstore/numerical/eigen.h>

#include <gemstore/model/model.h>
#include <gemstore/model/NRScreen.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/basis/soc.h>

#include <stdlib.h>

static inline double getmq(int index, argsModel_t *args_model)
{
    double mq;

    switch (index)
    {
    case 1:
        mq = args_model->mn;
        break;
    case 2:
        mq = args_model->ms;
        break;
    case 3:
        mq = args_model->mc;
        break;
    case 4:
        mq = args_model->mb;
        break;
    case 5:
        mq = args_model->mt;
        break;
    default:
        mq = 0;
        break;
    }

    return mq;
}

void spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, array_t *e_out, matrix_t *v_out, int v_len)
{
    /* construct basis */
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    double nu;
    for (int i = 0; i < nmax; i++) {
        nu = getnu(i + 1, nmax, rmax, rmin);
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = nu;
    }

    /* construc matrices */
    matrix_t mT = matrix_init(nmax, nmax);
    matrix_t mVconf = matrix_init(nmax, nmax);
    matrix_t mVcont = matrix_init(nmax, nmax);
    matrix_t mVsocm = matrix_init(nmax, nmax);
    matrix_t mVsotp = matrix_init(nmax, nmax);
    matrix_t mVtens = matrix_init(nmax, nmax);
    matrix_t Hfi = matrix_init(nmax, nmax);
    matrix_t Nfi = matrix_init(nmax, nmax);

    /* prepare variables */
    argsOrbit_t args_bra;
    argsOrbit_t args_ket;
    double factor;
    double factor_complex;
    argsModel_t args_model = argsNRScreen_meson;
    double s1 = 0.5, s2 = 0.5;
    double m1 = getmq(f1, &args_model);
    double m2 = getmq(f2, &args_model);
    double C12 = -4.0 / 3.0;
    argsModelDy_t args_dynmc = {
        .m1 = m1,
        .m2 = m2,
        .C12 = C12,
    };

    /* calculate matrix elements */
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            args_bra = basis[i];
            args_ket = basis[j];

            factor = 1 / sqrt(args_bra.scale + args_ket.scale);
            factor_complex =  sqrt(4 * args_bra.scale * args_ket.scale / (args_bra.scale + args_ket.scale));

            args_dynmc.OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc.OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc.OLS1 = operator_ldots1_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc.OLS2 = operator_ldots2_sl(s1, s2, S, L, s1, s2, S, L, J);
            args_dynmc.OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
                
            mT.value[i][j] = integral_matrix_element_complex(GRnlp_nonexp, NRScreen_T, factor_complex, &args_bra, &args_ket, &args_model, &args_dynmc);
            mVconf.value[i][j] = integral_matrix_element(GRnlr_nonexp, NRScreen_Vconf, factor, &args_bra, &args_ket, &args_model, &args_dynmc);
            mVcont.value[i][j] = integral_matrix_element(GRnlr_nonexp, NRScreen_Vcont, factor, &args_bra, &args_ket, &args_model, &args_dynmc);
            mVsocm.value[i][j] = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsocm, factor, &args_bra, &args_ket, &args_model, &args_dynmc);
            mVsotp.value[i][j] = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsotp, factor, &args_bra, &args_ket, &args_model, &args_dynmc);
            mVtens.value[i][j] = integral_matrix_element(GRnlr_nonexp, NRScreen_Vtens, factor, &args_bra, &args_ket, &args_model, &args_dynmc);
            Nfi.value[i][j] = integral_wfn_overlap(GRnlr_nonexp, factor, &args_bra, &args_ket);
        }
    }

    matrix_sum(&mT, &mVconf, &Hfi);
    matrix_sum(&Hfi, &mVcont, &Hfi);
    matrix_sum(&Hfi, &mVsocm, &Hfi);
    matrix_sum(&Hfi, &mVsotp, &Hfi);
    matrix_sum(&Hfi, &mVtens, &Hfi);

    int info = 0;    
    eigen_general_thread(Hfi.value, Nfi.value, nmax, e_out->value, v_out->value, v_len, &info);

    free(basis);
    matrix_free(&mT);
    matrix_free(&mVconf);
    matrix_free(&mVcont);
    matrix_free(&mVsocm);
    matrix_free(&mVsotp);
    matrix_free(&mVtens);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}