/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2026, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitccbar.h>
#include <gemstore/param/fitting.h>
#include <gemstore/param/typecc.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>

static const std::vector<State> experimental_data = {
    // charmonium (c c-bar)
    {3, 3, 1, 0, 0, 0, 2984.1,   5},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3637.8,   5},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3096.9,   5},   // J/ψ(1S)
    {3, 3, 2, 1, 0, 1, 3686.1,   5},   // ψ(2S)
    {3, 3, 1, 0, 1, 1, 3525.4,   5},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3414.7,   5},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3510.7,   5},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3556.2,   5},   // χc2(1P)
};

static double compute_chi2_GIScreen(const std::vector<double>& params, bool print_details)
{
	double chi_square = 0.0;
    DualStream dual("Fitting.out");

    for (const auto& state : experimental_data) {
        double e_out = call_meson_GIScreen(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 20.0, 0.1, params.data());
        double diff = 1000 * e_out - state.exp_mass;

        if (print_details) {
		    dual << "State (" << state.f1 << "," << state.f2 << "," << state.N << "," << state.S << "," << state.L << "," << state.J << ") calc = " 
            << 1000 * e_out << "  exp = " << state.exp_mass << "  diff = " << diff << std::endl;
        }
        double weight = (state.exp_error > 0.0) ? 1.0 / (state.exp_error * state.exp_error) : 1.0;
        chi_square += diff * diff * weight;
    }

    dual << "Total chi2=" << chi_square << std::endl;
    return chi_square;
}

static double compute_chi2_GIQuadra(const std::vector<double>& params, bool print_details)
{
	double chi_square = 0.0;
    DualStream dual("Fitting.out");

    for (const auto& state : experimental_data) {
        double e_out = call_meson_GIQuadra(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 20.0, 0.1, params.data());
        double diff = 1000 * e_out - state.exp_mass;

        if (print_details) {
		    dual << "State (" << state.f1 << "," << state.f2 << "," << state.N << "," << state.S << "," << state.L << "," << state.J << ") calc = " 
            << 1000 * e_out << "  exp = " << state.exp_mass << "  diff = " << diff << std::endl;
        }
        double weight = (state.exp_error > 0.0) ? 1.0 / (state.exp_error * state.exp_error) : 1.0;
        chi_square += diff * diff * weight;
    }

    dual << "Total chi2=" << chi_square << std::endl;
    return chi_square;
}

class Chi2GIScreen : public ROOT::Minuit2::FCNBase {
public:
	Chi2GIScreen(double error_def, size_t n_params, size_t n_data) {
        this->error_def_ = error_def;
        this->n_params_ = n_params;
        this->n_data_ = n_data;
    }
    ~Chi2GIScreen() {}

	double operator()(const std::vector<double>& params) const override {
		assert(params.size() == n_params_);
		return compute_chi2_GIScreen(params, false);
	}
	
	double Up() const override { return error_def_; }
	void SetErrorDef(double def) { error_def_ = def; }

private:
	double error_def_;			/* define error */
	size_t n_params_;			/* number of parameters */
	size_t n_data_;				/* number of data, DOF = n_data - n_params */
};

class Chi2GIQuadra : public ROOT::Minuit2::FCNBase {
public:
	Chi2GIQuadra(double error_def, size_t n_params, size_t n_data) {
        this->error_def_ = error_def;
        this->n_params_ = n_params;
        this->n_data_ = n_data;
    }
    ~Chi2GIQuadra() {}

	double operator()(const std::vector<double>& params) const override {
		assert(params.size() == n_params_);
		return compute_chi2_GIQuadra(params, false);
	}
	
	double Up() const override { return error_def_; }
	void SetErrorDef(double def) { error_def_ = def; }

private:
	double error_def_;			/* define error */
	size_t n_params_;			/* number of parameters */
	size_t n_data_;				/* number of data, DOF = n_data - n_params */
};

void minuit2_ccbar_GIScreen(double *params_out)
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
    //upar.Add("b2", 0.02, 0.01, 0.0, 0.1);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    int N_PARAMS = upar.Params().size();
    int N_DATA = experimental_data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2GIScreen minuit_fit(1.0, N_PARAMS, N_DATA);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2_GIScreen(min_result.UserParameters().Params(), true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}

void minuit2_ccbar_GIQuadra(double *params_out)
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
    int N_DATA = experimental_data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2GIQuadra minuit_fit(1.0, N_PARAMS, N_DATA);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2_GIQuadra(min_result.UserParameters().Params(), true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}