/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/meson.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_MESON = {
    // B (n b-bar) 2
    {1, 4, 1, 0, 0, 0, 5279.5,   1},   // B(1S)
    {1, 4, 1, 1, 0, 1, 5324.8,   1},   // B*(1S)

    // D (n c-bar) 2
    {1, 3, 1, 0, 0, 0, 1867.2,   2},   // D(1S)
    {1, 3, 1, 1, 0, 1, 2008.6,   2},   // D*(1S)

    // Bs (s b-bar) 2
    {2, 4, 1, 0, 0, 0, 5366.9,   1},   // Bs(1S)
    {2, 4, 1, 1, 0, 1, 5415.4,   2},   // Bs*(1S)

    // Ds (s c-bar) 2
    {2, 3, 1, 0, 0, 0, 1968.4,   1},   // Ds(1S)
    {2, 3, 1, 1, 0, 1, 2112.2,   1},   // Ds*(1S)

    // Bc (c b-bar) 2
    {3, 4, 1, 0, 0, 0, 6274.5,   1},   // Bc(1S)
    {3, 4, 2, 0, 0, 0, 6871.2,   1},   // Bc(2S)

    // charmonium (c c-bar) 17
    {3, 3, 1, 0, 0, 0, 2984.1,   1},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3637.8,   1},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3096.9,   1},   // ψ(1S)
    {3, 3, 2, 1, 0, 1, 3686.1,   1},   // ψ(2S)
    {3, 3, 3, 1, 0, 1, 4039.6,   4},   // ψ(3S)
    {3, 3, 4, 1, 0, 1, 4222.2,   3},   // ψ(4S)
    {3, 3, 5, 1, 0, 1, 4415.0,   5},   // ψ(5S)
    {3, 3, 1, 1, 2, 1, 3778.1,   1},   // ψ(1D)
    {3, 3, 1, 1, 2, 2, 3823.5,   1},   // ψ2(1D)
    {3, 3, 1, 1, 2, 3, 3842.7,   1},   // ψ3(1D)
    {3, 3, 2, 1, 2, 1, 4191.0,   5},   // ψ(2D)
    {3, 3, 3, 1, 2, 1, 4374.0,   7},   // ψ(3D)
    {3, 3, 1, 0, 1, 1, 3525.4,   1},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3414.7,   1},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3510.7,   1},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3556.2,   1},   // χc2(1P)
    {3, 3, 2, 1, 1, 2, 3922.5,   1},   // χc2(2P)

    // bottomonium (b b-bar) 19
    {4, 4, 1, 0, 0, 0, 9398.7,   2},   // ηb(1S)
    {4, 4, 2, 0, 0, 0, 9999.0,   4},   // ηb(2S)
    {4, 4, 1, 1, 0, 1, 9460.4,   1},   // Υ(1S)
    {4, 4, 2, 1, 0, 1, 10023.4,  1},   // Υ(2S)
    {4, 4, 3, 1, 0, 1, 10355.1,  1},   // Υ(3S)
    {4, 4, 4, 1, 0, 1, 10579.4,  2},   // Υ(4S)
    {4, 4, 5, 1, 0, 1, 10885.2,  3},   // Υ(5S)
    {4, 4, 6, 1, 0, 1, 11000.0,  4},   // Υ(6S)
    {4, 4, 1, 1, 2, 2, 10163.7,  2},   // Υ2(1D)
    {4, 4, 1, 0, 1, 1, 9899.3,   1},   // hb(1P)
    {4, 4, 2, 0, 1, 1, 10259.8,  2},   // hb(2P)
    {4, 4, 1, 1, 1, 0, 9859.4,   1},   // χb0(1P)
    {4, 4, 1, 1, 1, 1, 9892.8,   1},   // χb1(1P)
    {4, 4, 1, 1, 1, 2, 9912.2,   1},   // χb2(1P)
    {4, 4, 2, 1, 1, 0, 10232.5,  1},   // χb0(2P)
    {4, 4, 2, 1, 1, 1, 10255.5,  1},   // χb1(2P)
    {4, 4, 2, 1, 1, 2, 10268.7,  1},   // χb2(2P)
    {4, 4, 3, 1, 1, 1, 10513.4,  1},   // χb1(3P)
    {4, 4, 3, 1, 1, 2, 10524.0,  1}    // χb2(3P)
};

void minuit2_meson_GIScreen(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.220, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.419, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.628, 0.01, 1.5, 2.0);
    upar.Add("mb", 4.977, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.18, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_MESON, MODEL_GI_SCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_MESON, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}

void minuit2_meson_GIQuadra(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.220, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.419, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.628, 0.01, 1.5, 2.0);
    upar.Add("mb", 4.977, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.18, 0.01, 0.1, 0.3);
    upar.Add("b2", 0.02, 0.01, 0.0, 0.1);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_MESON, MODEL_GI_QUADRA, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_MESON, min_result.UserParameters().Params(), MODEL_GI_QUADRA, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}