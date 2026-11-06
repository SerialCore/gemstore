/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/qmcmodel.h>
#include <gemstore/model/bagmodel.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double QMC_HBARC = 0.1973269718; /* GeV·fm */
static const double QMC_GAMMA = 4.0;          /* spin × isospin of p+n */

/* 1 fm^{-1} = ħc GeV, so 1 fm^{-3} = (ħc)^3 GeV^3. */
static double qmc_rho_gev(double rho_fm)
{
    double h = QMC_HBARC;

    return rho_fm * h * h * h;
}

/* Scalar charge S_q = ∂ω/∂m of one κ = −1 quark. y = R m and Ω = sqrt(x² + y²).
 * This is not bag_momentum, which is the magnetic moment. */
static double qmc_scalar_charge(double m, double x, double r)
{
    double y = r * m;
    double omega = sqrt(x * x + y * y);
    double den = omega * (omega - 1.0) + 0.5 * y;

    if (den == 0.0) return 0.0;
    return (0.5 * omega + y * (omega - 1.0)) / den;
}

/* Color-free bag of three light quarks at this mn and B. The mean field always
 * uses this source. A CMI or strange spectrum is evaluated later, at the solved σ.
 * x is the light-quark boundary root. Pair bindings stay zero. */
static void qmc_source(const argsQMCModel_t *model, double mn, double B, double *mass, double *radius, double *x)
{
    argsInput_t input = {0};
    argsBAGModel_t bag = {0};
    bag_mit_t mit;

    input.nn = 3;
    input.basis = BASIS_NONE;
    bag.mn = mn;
    bag.ms = model->ms;
    bag.mc = model->mc;
    bag.mb = model->mb;
    bag.B = B;
    bag.Z = model->Z;
    bag_variant(&input, &bag, NULL, NULL, 0, NULL, NULL, 0, &mit);
    *mass = mit.mass;
    *radius = mit.R;
    *x = mit.x[0];
}

/* QMCSTAD returns B0 and sets dB = 0. QMCMOD uses
 * B/B0 = (1 − x)^δ with x = 4 g_σ^B σ / (δ M_N), and
 * dB/dσ = −(4 g_σ^B / M_N) B0 (1 − x)^{δ−1}.
 * x ≥ 1 is rejected: the bag constant has vanished. */
static double qmc_bag_constant(const argsInput_t *input, const argsQMCModel_t *model, double sigma, double m_free, double *dB)
{
    double x;
    double factor;

    *dB = 0.0;
    if (input->model != MODEL_QMCMOD) return model->B;
    if (!(m_free > 0.0) || !(model->delta > 0.0)) {
        fprintf(stderr, "QMCMOD needs a positive free nucleon mass and delta > 0.\n");
        exit(1);
    }
    x = 4.0 * model->g_sigma_B * sigma / (model->delta * m_free);
    if (x >= 1.0) {
        fprintf(stderr, "QMCMOD bag constant vanishes at sigma = %.6f GeV.\n", sigma);
        exit(1);
    }
    factor = pow(1.0 - x, model->delta);
    *dB = -(4.0 * model->g_sigma_B / m_free) * model->B * pow(1.0 - x, model->delta - 1.0);
    return model->B * factor;
}

/* Fermi–Dirac occupation at kinetic chemical potential ν. |z| > 40 is saturated.
 * An antiparticle is the same function at −ν. */
static double qmc_fermi(double energy, double nu, double temp)
{
    double z = (energy - nu) / temp;

    if (z > 40.0) return 0.0;
    if (z < -40.0) return 1.0;
    return 1.0 / (exp(z) + 1.0);
}

/* Finite-T sea at fixed ν. Net density uses f − f̄; scalar density, energy and
 * pressure use f + f̄. Simpson weights need an even panel count. The momentum
 * cutoff is max(15 T, k at ν + 15 T). T = 0 is qmc_fermi_zero, so this returns zeros. */
static void qmc_fermi_sea(double temp, double nu, double mass, double *rho, double *rho_s, double *eps, double *press)
{
    double kmax;
    double dk;
    double sum_rho = 0.0;
    double sum_s = 0.0;
    double sum_e = 0.0;
    double sum_p = 0.0;
    int n = 480;
    int i;

    if (temp == 0.0) {
        *rho = 0.0;
        *rho_s = 0.0;
        *eps = 0.0;
        *press = 0.0;
        return;
    }
    kmax = 15.0 * temp;
    if (nu + 15.0 * temp > mass) {
        double k_f = sqrt((nu + 15.0 * temp) * (nu + 15.0 * temp) - mass * mass);

        if (k_f > kmax) kmax = k_f;
    }
    dk = kmax / n;
    for (i = 0; i <= n; i++) {
        double k = i * dk;
        double e = sqrt(k * k + mass * mass);
        double occ = qmc_fermi(e, nu, temp);
        double anti = qmc_fermi(e, -nu, temp);
        double w = (i == 0 || i == n) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        double k2 = k * k;

        sum_rho += w * k2 * (occ - anti);
        sum_s += w * k2 * mass / e * (occ + anti);
        sum_e += w * k2 * e * (occ + anti);
        sum_p += w * k2 * k2 / e * (occ + anti);
    }
    sum_rho *= dk / 3.0;
    sum_s *= dk / 3.0;
    sum_e *= dk / 3.0;
    sum_p *= dk / 3.0;
    *rho = QMC_GAMMA / (2.0 * M_PI * M_PI) * sum_rho;
    *rho_s = QMC_GAMMA / (2.0 * M_PI * M_PI) * sum_s;
    *eps = QMC_GAMMA / (2.0 * M_PI * M_PI) * sum_e;
    *press = QMC_GAMMA / (6.0 * M_PI * M_PI) * sum_p;
}

/* T = 0 Fermi sphere, γ = 4. ef is sqrt(k_F² + M²). A negligible mass uses the
 * ultrarelativistic sea: ρ_s = 0 and P = ε/3, because the massive logarithm diverges. */
static void qmc_fermi_zero(double rho, double mass, double *rho_s, double *eps, double *press, double *ef)
{
    double k_f;
    double e_f;
    double ln;

    if (rho <= 0.0) {
        *rho_s = 0.0;
        *eps = 0.0;
        *press = 0.0;
        *ef = mass;
        return;
    }
    k_f = cbrt(6.0 * M_PI * M_PI * rho / QMC_GAMMA);
    if (mass < 1.0e-4) {
        *rho_s = 0.0;
        *eps = QMC_GAMMA * k_f * k_f * k_f * k_f / (8.0 * M_PI * M_PI);
        *press = (*eps) / 3.0;
        *ef = k_f;
        return;
    }
    e_f = sqrt(k_f * k_f + mass * mass);
    ln = log((e_f + k_f) / mass);
    *rho_s = QMC_GAMMA * mass / (4.0 * M_PI * M_PI) * (k_f * e_f - mass * mass * ln);
    *eps = QMC_GAMMA / (16.0 * M_PI * M_PI) * (k_f * e_f * (2.0 * k_f * k_f + mass * mass) - mass * mass * mass * mass * ln);
    *press = QMC_GAMMA / (48.0 * M_PI * M_PI) * (k_f * e_f * (2.0 * k_f * k_f - 3.0 * mass * mass) + 3.0 * mass * mass * mass * mass * ln);
    *ef = e_f;
}

/* Shared by the σ root and its inner ν root. The inner root reads m_star and
 * must not replace it. One solve uses the context at a time. */
typedef struct qmc_sigma_ctx {
    const argsInput_t *input;
    const argsQMCModel_t *model;
    double m_free;
    double rho;
    double temp;
    double g_n;
    double rho_s;
    double m_star;
    double radius;
} qmc_sigma_ctx_t;

/* ρ(ν) − ρ_B. The scalar density of that sea is kept in the context. */
static double qmc_nu_residual(double nu, void *ctx)
{
    qmc_sigma_ctx_t *q = ctx;
    double rho, rho_s, eps, press;

    qmc_fermi_sea(q->temp, nu, q->m_star, &rho, &rho_s, &eps, &press);
    q->rho_s = rho_s;
    return rho - q->rho;
}

/* Root of a sign change on [a, b]. Same-sign endpoints are rejected. */
static double qmc_bisect(double (*f)(double, void *), void *ctx, double a, double b)
{
    double fa = f(a, ctx);
    double fb = f(b, ctx);
    int i;

    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0) {
        fprintf(stderr, "QMC root is not bracketed on [%.6f, %.6f].\n", a, b);
        exit(1);
    }
    for (i = 0; i < 80; i++) {
        double mid = 0.5 * (a + b);
        double fm = f(mid, ctx);

        if (fa * fm <= 0.0) {
            b = mid;
            fb = fm;
        }
        else {
            a = mid;
            fa = fm;
        }
        if (fabs(b - a) < 1.0e-10) break;
    }
    return 0.5 * (a + b);
}

/* Source at this σ. g_n is minus the derivative of M* with respect to σ,
 * 3 g_σ S_q − (4π/3) R³ dB/dσ, by the envelope theorem.
 * There is no second rearrangement term. mn − g_σ σ has to stay non-negative. */
static void qmc_source_response(qmc_sigma_ctx_t *q, double sigma)
{
    double dB = 0.0;
    double B = qmc_bag_constant(q->input, q->model, sigma, q->m_free, &dB);
    double mn = q->model->mn - q->model->g_sigma * sigma;
    double x;

    if (mn < 0.0) {
        fprintf(stderr, "In-medium light-quark mass is negative at sigma = %.6f GeV.\n", sigma);
        exit(1);
    }
    qmc_source(q->model, mn, B, &q->m_star, &q->radius, &x);
    q->g_n = 3.0 * q->model->g_sigma * qmc_scalar_charge(mn, x, q->radius) - (4.0 * M_PI / 3.0) * q->radius * q->radius * q->radius * dB;
}

/* m_σ² σ + (λ₃/2) σ² − g_N ρ_s. T = 0 uses the Fermi sphere. At T > 0 and
 * ρ_B = 0 the sea is integrated at ν = 0. At T > 0 and ρ_B > 0, ν is solved
 * on [0, M* + 2] before ρ_s is read. */
static double qmc_sigma_residual(double sigma, void *ctx)
{
    qmc_sigma_ctx_t *q = ctx;
    double nu;
    double rho_s;
    double eps;
    double press;
    double e_f;
    double slope = q->model->m_sigma * q->model->m_sigma * sigma + 0.5 * q->model->lambda3 * sigma * sigma;

    qmc_source_response(q, sigma);
    if (q->temp == 0.0) {
        qmc_fermi_zero(q->rho, q->m_star, &rho_s, &eps, &press, &e_f);
        q->rho_s = rho_s;
    }
    else if (q->rho == 0.0) {
        qmc_fermi_sea(q->temp, 0.0, q->m_star, &nu, &q->rho_s, &eps, &press);
    }
    else {
        nu = qmc_bisect(qmc_nu_residual, q, 0.0, q->m_star + 2.0);
        (void)nu;
    }
    return slope - q->g_n * q->rho_s;
}

/* Spectrum bag: in-medium light mass and B(σ). ms, mc, mb and Z stay free.
 * Pair bindings stay zero. Quark content is supplied by the caller. */
static argsBAGModel_t qmc_spectrum_bag(const argsQMCModel_t *model, double mn, double B)
{
    argsBAGModel_t bag = {0};

    bag.mn = mn;
    bag.ms = model->ms;
    bag.mc = model->mc;
    bag.mb = model->mb;
    bag.B = B;
    bag.Z = model->Z;
    return bag;
}

void qmc_solve(const argsInput_t *input, const argsQMCModel_t *model, qmc_matter_t *matter, argsBAGModel_t *bag)
{
    qmc_sigma_ctx_t ctx;
    double rho = qmc_rho_gev(input->rho_b);
    double sigma = 0.0;
    double sigma_hi;
    double dB = 0.0;
    double B;
    double mn;
    double e_f = 0.0;
    double eps_k = 0.0;
    double press_k = 0.0;
    double nu = 0.0;
    double omega;
    double u_sigma;

    if (model->g_omega != 0.0 && !(model->m_omega > 0.0)) {
        fprintf(stderr, "QMC needs m_omega > 0 when g_omega is nonzero.\n");
        exit(1);
    }
    if (model->g_sigma != 0.0 && !(model->m_sigma > 0.0)) {
        fprintf(stderr, "QMC needs m_sigma > 0 when g_sigma is nonzero.\n");
        exit(1);
    }

    ctx.input = input;
    ctx.model = model;
    ctx.rho = rho;
    ctx.temp = input->temp;
    /* Free source at σ = 0 and B = B0. QMCMOD uses this mass as M_N, and the
     * binding energy is measured from it. The boundary root is discarded. */
    qmc_source(model, model->mn, model->B, &ctx.m_free, &ctx.radius, &e_f);
    matter->m_free = ctx.m_free;

    /* σ stays zero when g_σ = 0, and when both density and temperature vanish.
     * The search window keeps m* > 0 and, for QMCMOD, B > 0. A root requires
     * a non-positive residual at zero and a non-negative residual at sigma_hi. */
    if (model->g_sigma != 0.0 && (rho > 0.0 || input->temp > 0.0)) {
        sigma_hi = 0.98 * model->mn / model->g_sigma;
        if (input->model == MODEL_QMCMOD) {
            double sigma_b = 0.98 * model->delta * ctx.m_free / (4.0 * model->g_sigma_B);

            if (sigma_b < sigma_hi) sigma_hi = sigma_b;
        }
        if (!(sigma_hi > 0.0)) {
            fprintf(stderr, "QMC cannot lower the light-quark mass. mn and g_sigma leave no room for sigma.\n");
            exit(1);
        }
        if (qmc_sigma_residual(0.0, &ctx) > 0.0 || qmc_sigma_residual(sigma_hi, &ctx) < 0.0) {
            fprintf(stderr, "QMC scalar field is not bracketed. Couplings do not balance at this density.\n");
            exit(1);
        }
        sigma = qmc_bisect(qmc_sigma_residual, &ctx, 0.0, sigma_hi);
    }
    /* Kinetics at the solved σ. ω is fixed by the net density and does not
     * enter the bag or the scalar equation. */
    qmc_source_response(&ctx, sigma);
    if (input->temp == 0.0) qmc_fermi_zero(rho, ctx.m_star, &ctx.rho_s, &eps_k, &press_k, &e_f);
    else if (rho == 0.0) {
        double rho_net;

        nu = 0.0;
        qmc_fermi_sea(input->temp, 0.0, ctx.m_star, &rho_net, &ctx.rho_s, &eps_k, &press_k);
    }
    else {
        double rho_net;

        nu = qmc_bisect(qmc_nu_residual, &ctx, 0.0, ctx.m_star + 2.0);
        qmc_fermi_sea(input->temp, nu, ctx.m_star, &rho_net, &ctx.rho_s, &eps_k, &press_k);
    }

    B = qmc_bag_constant(input, model, sigma, matter->m_free, &dB);
    mn = model->mn - model->g_sigma * sigma;
    omega = (model->m_omega > 0.0) ? 3.0 * model->g_omega * rho / (model->m_omega * model->m_omega) : 0.0;
    /* U = ½ m_σ² σ² + (λ₃/6) σ³, so ∂U/∂σ is the slope in the σ residual.
     * ε = ε_kin + U + ½ m_ω² ω² and P = P_kin − U + ½ m_ω² ω².
     * μ adds 3 g_ω ω to E_F* at T = 0, and to ν at T > 0. ρ_3 stays zero. */
    u_sigma = 0.5 * model->m_sigma * model->m_sigma * sigma * sigma + model->lambda3 * sigma * sigma * sigma / 6.0;

    matter->field.model = input->model;
    matter->field.sigma = sigma;
    matter->field.omega = omega;
    matter->field.rho3 = 0.0;
    matter->field.B = B;
    matter->field.mn_star = mn;
    matter->m_nucleon = ctx.m_star;
    matter->R_nucleon = ctx.radius;
    matter->rho_s = ctx.rho_s;
    matter->mu = ((input->temp == 0.0) ? e_f : nu) + 3.0 * model->g_omega * omega;
    matter->eps = eps_k + u_sigma + 0.5 * model->m_omega * model->m_omega * omega * omega;
    matter->press = press_k - u_sigma + 0.5 * model->m_omega * model->m_omega * omega * omega;
    matter->ebind = (rho > 0.0) ? matter->eps / rho - matter->m_free : 0.0;
    /* Same mn* and B as the source. The spectrum quark content comes from input. */
    *bag = qmc_spectrum_bag(model, mn, B);
}
