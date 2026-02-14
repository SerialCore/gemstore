/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2026, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/fitting.h>
#include <gemstore/entry.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>

struct State {
    int f1, f2, N, S, L, J;		/* quantum numbers */
	double exp_mass;			/* experimental mass */
	double exp_error;			/* experimental error */
};

const std::vector<State> experimental_data = {
    // Bc (c b-bar)
    {3, 4, 1, 0, 0, 0, 6.2745,   0.001},   // Bc(1S)
    {3, 4, 2, 0, 0, 0, 6.8712,   0.005},   // Bc(2S)  — large error

    // Bs (s b-bar)
    {2, 4, 1, 0, 0, 0, 5.3669,   0.0005},   // Bs(1S)
    {2, 4, 1, 1, 0, 1, 5.4154,   0.001},   // Bs*(1S) vector

    // Ds (s c-bar)
    {2, 3, 1, 0, 0, 0, 1.9684,   0.0005},   // Ds(1S)
    {2, 3, 1, 1, 0, 1, 2.1122,   0.001},   // Ds*(1S)

    // B (u/d b-bar)
    {1, 4, 1, 0, 0, 0, 5.2796,   0.0003},   // B(1S) average
    {1, 4, 1, 1, 0, 1, 5.3248,   0.0008},   // B*(1S)

    // D (u/d c-bar)
    {1, 3, 1, 0, 0, 0, 1.8648,   0.0003},   // D(1S) average
    {1, 3, 1, 1, 0, 1, 2.0069,   0.0008},   // D*(1S)

    // charmonium (c c-bar)
    {3, 3, 1, 0, 0, 0, 2.9841,   0.002},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3.6378,   0.003},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3.0969,   0.0001},   // J/ψ(1S)
    {3, 3, 2, 1, 0, 1, 3.6861,   0.0005},   // ψ(2S)
    {3, 3, 1, 0, 1, 1, 3.5254,   0.001},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3.4147,   0.001},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3.5107,   0.0005},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3.5562,   0.0005},   // χc2(1P)

    // bottomonium (b b-bar)
    {4, 4, 1, 0, 0, 0, 9.3987,   0.002},   // ηb(1S)
    {4, 4, 2, 0, 0, 0, 9.9990,   0.005},   // ηb(2S) — large error
    {4, 4, 1, 1, 0, 1, 9.4604,   0.0005},   // Υ(1S)
    {4, 4, 2, 1, 0, 1, 10.0234,  0.001},   // Υ(2S)
    {4, 4, 3, 1, 0, 1, 10.3551,  0.001},   // Υ(3S)
    {4, 4, 4, 1, 0, 1, 10.5794,  0.002},   // Υ(4S)
    {4, 4, 1, 1, 2, 2, 10.1637,  0.003},   // Υ(1D₂)
    {4, 4, 1, 0, 1, 1, 9.8993,   0.001},   // hb(1P)
    {4, 4, 2, 0, 1, 1, 10.2598,  0.0015},   // hb(2P)
    {4, 4, 1, 1, 1, 0, 9.8594,   0.001},   // χb0(1P)
    {4, 4, 1, 1, 1, 1, 9.8928,   0.0008},   // χb1(1P)
    {4, 4, 1, 1, 1, 2, 9.9122,   0.0008},   // χb2(1P)
    {4, 4, 2, 1, 1, 0, 10.2325,  0.001},   // χb0(2P)
    {4, 4, 2, 1, 1, 1, 10.2555,  0.001},   // χb1(2P)
    {4, 4, 2, 1, 1, 2, 10.2687,  0.001},   // χb2(2P)
    {4, 4, 3, 1, 1, 1, 10.5134,  0.0015},   // χb1(3P)
    {4, 4, 3, 1, 1, 2, 10.5240,  0.002}    // χb2(3P)
};

const int N_PARAMS = 9;		/* number of parameters to be fitted */

double compute_chi2(const std::vector<double>& params)
{
	double chi_square = 0.0;

    for (const auto& state : experimental_data) {
        double e_out = call_fitting_meson_NRScreen(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 10.0, 0.1, params.data());

		//std::cout << "State (" << state.f1 << "," << state.f2 << ") calc = " << e_out << "  exp = " << state.exp_mass << std::endl;
        double diff = e_out - state.exp_mass;
        double weight = (state.exp_error > 0.0) ? 1.0 / (state.exp_error * state.exp_error) : 1.0;
        chi_square += diff * diff * weight;
    }

    std::cout << "Total chi2=" << chi_square << std::endl;
    return chi_square;
}

class Chi2Functor : public ROOT::Minuit2::FCNBase {
public:
	Chi2Functor() : error_def_(1.0), n_params_(N_PARAMS), n_data_(experimental_data.size()) {}
    ~Chi2Functor() {}

	double operator()(const std::vector<double>& params) const override {
		assert(params.size() == n_params_);
		return compute_chi2(params);
	}
	
	double Up() const override { return error_def_; }
	void SetErrorDef(double def) { error_def_ = def; }

private:
	double error_def_;			/* define error */
	size_t n_params_;			/* number of parameters */
	size_t n_data_;				/* number of data, DOF = n_data - n_params */
};

void perform_fit(double *params_out)
{
    srand(time(0));

    Chi2Functor minuit_fit;

    ROOT::Minuit2::MnUserParameters upar;
    double step = 0.01;  	/* initial guess of step */
    /*double err = 0.1;    	 initial guess of error */

    upar.Add("mn", 0.606, step, 0.1, 1.0);
    upar.Add("ms", 0.780, step, 0.3, 1.0);
    upar.Add("mc", 1.984, step, 1.0, 3.0);
    upar.Add("mb", 5.368, step, 4.0, 6.0);
    upar.Add("as", 0.3930, step, 0.1, 0.5);
    upar.Add("b1", 0.2312, step, 0.1, 0.3);
    upar.Add("mu", 0.069, step, 0.01, 0.1);
    upar.Add("c", -1.1711, step, -2.0, 0.0);
    upar.Add("sig", 1.842, step, 1.0, 3.0);

    /* use of Migrad algorithm with strategy 2, high precision */
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    std::cout << "Fit converged: " << min_result.IsValid() << std::endl;
    std::cout << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}