/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameters.h>

#include <iostream>
#include <vector>
#include <cstdlib>

namespace {

class Chi2Minimizer : public ROOT::Minuit2::FCNBase {
public:
    Chi2Minimizer(const argsInput_t *job, const fit_target_t *target, double error_def)
        : job_(job), target_(target), error_def_(error_def) {}

    double operator()(const std::vector<double> &params) const override
    {
        if (params.size() != (size_t)target_->nparams) {
            std::cerr << "Minuit parameter vector size does not match the fit table\n";
            std::exit(1);
        }
        return fit_chi2(job_, target_, params.data(), NULL, NULL);
    }

    double Up() const override { return error_def_; }

private:
    const argsInput_t *job_;
    const fit_target_t *target_;
    double error_def_;
};

} /* namespace */

void fit_minuit_run(const argsInput_t *job, const fit_target_t *target, fit_result_t *out)
{
    ROOT::Minuit2::MnUserParameters upar;
    int i;

    out->value = (double *)calloc((size_t)target->nparams, sizeof(double));
    out->error = (double *)calloc((size_t)target->nparams, sizeof(double));
    out->mass_mev = (double *)calloc((size_t)target->nstates, sizeof(double));
    out->residual = (double *)calloc((size_t)target->nstates, sizeof(double));
    if (out->value == NULL || out->error == NULL || out->mass_mev == NULL || out->residual == NULL) {
        std::cerr << "Out of memory\n";
        std::exit(1);
    }

    for (i = 0; i < target->nparams; i++) {
        const fit_param_t *param = &target->params[i];
        upar.Add(param->name, param->value, param->step, param->min, param->max);
        if (param->fixed) {
            upar.Fix(param->name);
        }
    }

    Chi2Minimizer chi2(job, target, 1.0);
    ROOT::Minuit2::MnMigrad migrad(chi2, upar, ROOT::Minuit2::MnStrategy{2});
    ROOT::Minuit2::FunctionMinimum minimum = migrad();
    std::vector<double> values = minimum.UserParameters().Params();

    if (values.size() != (size_t)target->nparams) {
        std::cerr << "Minuit result size does not match the fit table\n";
        std::exit(1);
    }
    for (i = 0; i < target->nparams; i++) {
        out->value[i] = values[(size_t)i];
        out->error[i] = target->params[i].fixed ? 0.0
            : minimum.UserParameters().Error((unsigned int)i);
    }
    out->edm = minimum.Edm();
    out->valid = minimum.IsValid() ? 1 : 0;
    out->chi2 = fit_chi2(job, target, out->value, out->mass_mev, out->residual);
    std::cout << minimum.UserParameters() << std::endl;
}
