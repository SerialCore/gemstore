/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitmeson.h>
#include <gemstore/param/helper.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> data = {
    // K (u/d s-bar)
    {1, 2, 1, 0, 0, 0, 497.6,    5},   // K(1S)
    {1, 2, 1, 1, 0, 1, 895.6,    5},   // K*(1S)

    // phi (s s-bar)
    {2, 2, 1, 1, 0, 1, 1019.5,   5},   // phi(1S)

    // B (u/d b-bar)
    {1, 4, 1, 0, 0, 0, 5279.6,   5},   // B(1S)
    {1, 4, 1, 1, 0, 1, 5324.8,   5},   // B*(1S)

    // D (u/d c-bar)
    {1, 3, 1, 0, 0, 0, 1864.8,   5},   // D(1S)
    {1, 3, 1, 1, 0, 1, 2006.9,   5},   // D*(1S)

    // Bs (s b-bar)
    {2, 4, 1, 0, 0, 0, 5366.9,   5},   // Bs(1S)
    {2, 4, 1, 1, 0, 1, 5415.4,   5},   // Bs*(1S)

    // Ds (s c-bar)
    {2, 3, 1, 0, 0, 0, 1968.4,   5},   // Ds(1S)
    {2, 3, 1, 1, 0, 1, 2112.2,   5},   // Ds*(1S)

    // Bc (c b-bar)
    {3, 4, 1, 0, 0, 0, 6274.5,   5},   // Bc(1S)
    {3, 4, 2, 0, 0, 0, 6871.2,   5},   // Bc(2S)

    // charmonium (c c-bar)
    {3, 3, 1, 0, 0, 0, 2984.1,   5},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3637.8,   5},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3096.9,   5},   // ψ(1S)
    {3, 3, 2, 1, 0, 1, 3686.1,   5},   // ψ(2S)
    {3, 3, 1, 0, 1, 1, 3525.4,   5},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3414.7,   5},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3510.7,   5},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3556.2,   5},   // χc2(1P)

    // bottomonium (b b-bar)
    {4, 4, 1, 0, 0, 0, 9398.7,   5},   // ηb(1S)
    {4, 4, 2, 0, 0, 0, 9999.0,   5},   // ηb(2S)
    {4, 4, 1, 1, 0, 1, 9460.4,   5},   // Υ(1S)
    {4, 4, 2, 1, 0, 1, 10023.4,  5},   // Υ(2S)
    {4, 4, 3, 1, 0, 1, 10355.1,  5},   // Υ(3S)
    {4, 4, 4, 1, 0, 1, 10579.4,  5},   // Υ(4S)
    {4, 4, 1, 1, 2, 2, 10163.7,  5},   // Υ(1D₂)
    {4, 4, 1, 0, 1, 1, 9899.3,   5},   // hb(1P)
    {4, 4, 2, 0, 1, 1, 10259.8,  5},   // hb(2P)
    {4, 4, 1, 1, 1, 0, 9859.4,   5},   // χb0(1P)
    {4, 4, 1, 1, 1, 1, 9892.8,   5},   // χb1(1P)
    {4, 4, 1, 1, 1, 2, 9912.2,   5},   // χb2(1P)
    {4, 4, 2, 1, 1, 0, 10232.5,  5},   // χb0(2P)
    {4, 4, 2, 1, 1, 1, 10255.5,  5},   // χb1(2P)
    {4, 4, 2, 1, 1, 2, 10268.7,  5},   // χb2(2P)
    {4, 4, 3, 1, 1, 1, 10513.4,  5},   // χb1(3P)
    {4, 4, 3, 1, 1, 2, 10524.0,  5}    // χb2(3P)
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
    int N_DATA = data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(data, MODEL_GI_SCREEN, N_PARAMS, N_DATA, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(data, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
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
    int N_DATA = data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(data, MODEL_GI_QUADRA, N_PARAMS, N_DATA, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(data, min_result.UserParameters().Params(), MODEL_GI_QUADRA, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}