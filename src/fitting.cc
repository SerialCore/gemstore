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
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>

class DualStream {
public:
    DualStream(const std::string& filename = "Fitting.out")
        : console(std::cout), file(filename, std::ios::out | std::ios::app) {
        if (!file.is_open()) {
            std::cerr << "Cannot open log file: " << filename << std::endl;
        }
    }

    DualStream(const DualStream&) = delete;
    DualStream& operator=(const DualStream&) = delete;
    DualStream(DualStream&&) = delete;
    DualStream& operator=(DualStream&&) = delete;

    ~DualStream() {
        if (file.is_open()) {
            file << std::endl;
            file.close();
        }
    }

    template<typename T>
    DualStream& operator<<(const T& value) {
        console << value;
        if (file.is_open()) {
            file << value;
        }
        return *this;
    }

    DualStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        console << manip;
        if (file.is_open()) {
            file << manip;
        }
        return *this;
    }

private:
    std::ostream& console;
    std::ofstream file;
};

struct State {
    int f1, f2, N, S, L, J;		/* quantum numbers */
	double exp_mass;			/* experimental mass */
	double exp_error;			/* experimental error */
};

const std::vector<State> experimental_data = {
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
    {3, 3, 1, 1, 0, 1, 3096.9,   5},   // J/ψ(1S)
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

double compute_chi2(const std::vector<double>& params, bool print_details)
{
	double chi_square = 0.0;
    DualStream dual("Fitting.out");

    for (const auto& state : experimental_data) {
        double e_out = call_fitting_meson_GIScreen(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 10.0, 0.1, params.data());
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

class Chi2Functor : public ROOT::Minuit2::FCNBase {
public:
	Chi2Functor(double error_def, size_t n_params, size_t n_data) {
        this->error_def_ = error_def;
        this->n_params_ = n_params;
        this->n_data_ = n_data;
    }
    ~Chi2Functor() {}

	double operator()(const std::vector<double>& params) const override {
		assert(params.size() == n_params_);
		return compute_chi2(params, false);
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
    Chi2Functor minuit_fit(1.0, N_PARAMS, N_DATA);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(min_result.UserParameters().Params(), true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}