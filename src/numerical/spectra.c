/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/spectra.h>
#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/integral.h>

#include <gemstore/model/model.h>
#include <gemstore/model/NRScreen.h>

#include <gemstore/basis/orbit.h>
#include <gemstore/basis/soc.h>

void spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, array_t *e_out, matrix_t *v_out)
{
    double s1 = 0.5, s2 = 0.5;

    double nu1 = getnu(1, 30, 25, 0.1);
    argsOrbit_t args_bra = {
        .n = 1,
        .l = 0,
        .scale = nu1
    };

    double nu2 = getnu(2, 30, 25, 0.1);
    argsOrbit_t args_ket = {
        .n = 2,
        .l = 0,
        .scale = nu2
    };

    argsModel_t args_model = argsNRScreen_meson;
    double m1 = args_model.mc;
    double m2 = args_model.mc;
    argsFlavor_t args_flavor = {
        .m1 = m1,
        .m2 = m2
    };

    double cen = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
    double sds = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
    double ls1 = operator_ldots1_sl(s1, s2, S, L, s1, s2, S, L, J);
    double ls2 = operator_ldots2_sl(s1, s2, S, L, s1, s2, S, L, J);
    double ten = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
    argsSOC_t args_soc = {
        .OCent = cen,
        .OSdS = sds,
        .OLS1 = ls1,
        .OLS2 = ls2,
        .OTens = ten
    };

    double factor = 1 / sqrt(nu1 + nu2);
    double factor_complex =  sqrt(4 * nu1 * nu2 / (nu1 + nu2));
    double element;

    element = integral_matrix_element_complex(GRnlp_nonexp, NRScreen_T, factor_complex, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vconf, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vcont, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsocm, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsotp, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vtens, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);

    element++;
}