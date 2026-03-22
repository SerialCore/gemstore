/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_TYPECC
#define GEMSTORE_PARAM_TYPECC

#include <gemstore/param/fitting.h>
#include <gemstore/types.h>

#include <Minuit2/FCNBase.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

struct State {
    int f1, f2, N, S, L, J;		/* quantum numbers */
	double exp_mass;			/* experimental mass */
	double exp_error;			/* experimental error */
};

/* DualStream is to print console and file at once */
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

/* The entrance of chi2 computation with data, parameters and model selection inputed.
 * Should be static to be compiled for each fitting process. */
static double compute_chi2(const std::vector<State> data, const std::vector<double>& params, model_type_t model, bool print_details)
{
	double chi_square = 0.0;
    DualStream dual("Fitting.out");

    for (const auto& state : data) {
        double e_out;
        switch (model)
        {
        case MODEL_GI_SCREEN:
            e_out = call_meson_GIScreen(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 20.0, 0.01, params.data());
            break;
        case MODEL_GI_QUADRA:
            e_out = call_meson_GIQuadra(state.f1, state.f2, state.N, state.S, state.L, state.J, 20, 20.0, 0.01, params.data());
            break;
        default:
            e_out = 0.0;
            break;
        }
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

/* The helper class containing base-data and model selection */
class Chi2Minimizer : public ROOT::Minuit2::FCNBase {
public:
	Chi2Minimizer(std::vector<State> data, model_type_t model, size_t n_params, double error_def) {
        this->data_ = data;
        this->model_ = model;
        this->n_params_ = n_params;
        this->error_def_ = error_def;
    }
    ~Chi2Minimizer() {}

	double operator()(const std::vector<double>& params) const override {
		assert(params.size() == n_params_);
		return compute_chi2(data_, params, model_, false);
	}
	
	double Up() const override { return error_def_; }
	void SetErrorDef(double def) { error_def_ = def; }

private:
    std::vector<State> data_;   /* define input data */
    model_type_t model_;          /* define model type */
	double error_def_;			/* define error */
	size_t n_params_;			/* number of parameters */
};

#endif