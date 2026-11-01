/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Baryon GI radial moments. K(x) is the gimodel potential; the
 * quadrature is integral_exp_rn. Overlap and <x^2> stay closed forms.
 */

#include <gemstore/model/giradial.h>
#include <gemstore/basis/channel.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/math/integral.h>

#include <math.h>

void scdk_vargs_from_model(scdk_vargs_t *varg, const argsGIModel_t *model, model_type_t mtype)
{
    varg->mn = model->mn;
    varg->ms = model->ms;
    varg->mc = model->mc;
    varg->mb = model->mb;
    varg->alpha[0] = GI_ALPHA_K[0];
    varg->alpha[1] = GI_ALPHA_K[1];
    varg->alpha[2] = GI_ALPHA_K[2];
    varg->gamma[0] = GI_GAMMA_K[0];
    varg->gamma[1] = GI_GAMMA_K[1];
    varg->gamma[2] = GI_GAMMA_K[2];
    varg->b = model->b;
    varg->c = model->c;
    varg->sigma0 = model->sigma_0;
    varg->s = model->s;
    varg->f = 0.5;
    varg->mu = model->mu;
    varg->econt = model->epsilon_cont;
    varg->etens = model->epsilon_tens;
    varg->esov = model->epsilon_sov;
    varg->esos = model->epsilon_sos;
    varg->eCoul = 0.0;
    varg->model = mtype;
}

/* Spin-angular weight is ±1. spinang supplies the real coefficient.
 * quark != 0 selects √(p²+m²) for quark c; otherwise mi, mj are the pair. */
static double gi_moment(potential_t pot, double b11, int n, basis_3b qf, int c,
    const scdk_vargs_t *varg, int quark, system_type_t system,
    double olsi, double olsj, double osds, double otens)
{
    argsGIModel_t model = {0};
    argsGIModelDy_t dyn = {0};
    gi_pot_ctx_t ctx;
    double mi, mj, mk;

    model.b = varg->b;
    model.mu = varg->mu;
    model.c = varg->c;
    model.sigma_0 = varg->sigma0;
    model.s = varg->s;
    model.epsilon_cont = varg->econt;
    model.epsilon_sov = varg->esov;
    model.epsilon_sos = varg->esos;
    model.epsilon_tens = varg->etens;

    if (quark) {
        if (c == 1) mi = qf.m1;
        else if (c == 2) mi = qf.m2;
        else if (c == 3) mi = qf.m3;
        else mi = 0.0;
        mj = mi;
    } else {
        channel_3b_ijk(qf.m1, qf.m2, qf.m3, &mi, &mj, &mk, c);
        (void)mk;
    }

    dyn.model = varg->model;
    dyn.system = system;
    dyn.mi = mi;
    dyn.mj = mj;
    dyn.Cij = -2.0 / 3.0;
    dyn.OCent = 1.0;
    dyn.OSdS = osds;
    dyn.OLSi = olsi;
    dyn.OLSj = olsj;
    dyn.OTens = otens;
    dyn.Sigij = sigma_ij(mi, mj, varg->sigma0, varg->s);
    sigma_k_ij(dyn.Sigij, dyn.Sigkij);
    ctx.model = &model;
    ctx.dyn = &dyn;
    return integral_exp_rn(pot, b11, n, &ctx);
}

/* ∫_0^∞ r^n exp(-a r^2) dr = (1/2) a^{-(n+1)/2} Γ((n+1)/2).
 * GammaSemiInt[n+1] stores that gamma value. */
static double InteExprn(double a,int n)
{
	double GammaSemiInt[100]={
 0.00000000000000000000e+000, 1.77245385090551602729e+000, 1.00000000000000000000e+000, 8.86226925452758013649e-001, 1.00000000000000000000e+000, 1.32934038817913702047e+000, 2.00000000000000000000e+000, 3.32335097044784255118e+000, 6.00000000000000000000e+000, 1.16317283965674489291e+001, 2.40000000000000000000e+001, 5.23427777845535201811e+001, 1.20000000000000000000e+002, 2.87885277815044360996e+002, 7.20000000000000000000e+002, 1.87125430579778834647e+003, 5.04000000000000000000e+003, 1.40344072934834125985e+004, 4.03200000000000000000e+004, 1.19292461994609007087e+005,
 3.62880000000000000000e+005, 1.13327838894878556733e+006, 3.62880000000000000000e+006, 1.18994230839622484570e+007, 3.99168000000000000000e+007, 1.36843365465565857255e+008, 4.79001600000000000000e+008, 1.71054206831957321569e+009, 6.22702080000000000000e+009, 2.30923179223142384118e+010, 8.71782912000000000000e+010, 3.34838609873556456972e+011, 1.30767436800000000000e+012, 5.18999845304012508307e+012, 2.09227898880000000000e+013, 8.56349744751620638707e+013, 3.55687428096000000000e+014, 1.49861205331533611773e+015, 6.40237370572800000000e+015, 2.77243229863337181781e+016,
 1.21645100408832000000e+017, 5.40624298233507504473e+017, 2.43290200817664000000e+018, 1.10827981137869038417e+019, 5.10909421717094400000e+019, 2.38280159446418432596e+020, 1.12400072777760768000e+021, 5.36130358754441473342e+021, 2.58520167388849766400e+022, 1.25990634307293746235e+023, 6.20448401733239439360e+023, 3.08677054052869678277e+024, 1.55112100433309859840e+025, 7.87126487834817679606e+025, 4.03291461126605635584e+026, 2.08588519276226685095e+027, 1.08888694504183521607e+028, 5.73618428009623384013e+028, 3.04888344611713860501e+029, 1.63481251982742664443e+030,
 8.84176199373970195454e+030, 4.82269693349090860109e+031, 2.65252859812191058636e+032, 1.47092256471472712333e+033, 8.22283865417792281772e+033, 4.63340607885139043849e+034, 2.63130836933693530167e+035, 1.50585697562670189251e+036, 8.68331761881188649551e+036, 5.04462086834945133991e+037, 2.95232799039604140847e+038, 1.74039419958056071227e+039, 1.03331479663861449296e+040, 6.17839940851099052856e+040, 3.71993326789901217468e+041, 2.25511578410651154292e+042, 1.37637530912263450463e+043, 8.45668419039941828596e+043, 5.23022617466601111760e+044, 3.25582341330377604009e+045,
 2.03978820811974433586e+046, 1.28605024825499153583e+047, 8.15915283247897734345e+047, 5.20850350543271572014e+048, 3.34525266131638071081e+049, 2.16152895475457702386e+050, 1.40500611775287989854e+051, 9.18649805770695235140e+051, 6.04152630633738356373e+052, 3.99612665510252427286e+053, 2.65827157478844876804e+054, 1.77827636152062330142e+055, 1.19622220865480194561e+056, 8.09115744491883602147e+056, 5.50262215981208894985e+057, 3.76238821188725874998e+058, 2.58623241511168180642e+059, 1.78713440064644790624e+060, 1.24139155925360726708e+061, 8.66760184313527234528e+061};
	return 0.5e0*pow(a,-0.5*(1+n))*GammaSemiInt[n+1];
}

double inteNfi(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qf;
    (void)qi;
    (void)ctx;
    (void)c;
    return InteExprn(b11, n);
}

double inteVogeG(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVcoul, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double inteVcont(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVcont, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 1, 0);
}

double inteVtens(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVtens, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 1);
}

double inteVsovii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsovi, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 1, 0, 0, 0);
}

double inteVsovjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsovj, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 1, 0, 0);
}

/* GIVsovij folds both legs. Split them with the meson sign and OLSi = ±1. */
double inteVsovij(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsovij, b11, n, qf, c, ctx, 0, SYSTEM_MESON, 1, 0, 0, 0);
}

double inteVsovji(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsovij, b11, n, qf, c, ctx, 0, SYSTEM_MESON, -1, 0, 0, 0);
}

double inteVstring(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVconf, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double inteVsosii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsosi, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 1, 0, 0, 0);
}

double inteVsosjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVsosj, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 1, 0, 0);
}

double intepogeG(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVbetaijcoul, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepcont(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaijcont, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double inteptens(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaijtens, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsovii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaiisov, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsovjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltajjsov, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsovji(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaijsov, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsovij(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaijsov, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsosii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltaiisos, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double intepsosjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVdeltajjsos, b11, n, qf, c, ctx, 0, SYSTEM_BARYON, 0, 0, 0, 0);
}

double inteTi(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qi;
    return gi_moment(GIVt_quark, b11, n, qf, c, ctx, 1, SYSTEM_BARYON, 0, 0, 0, 0);
}

double inteRMS(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    (void)qf;
    (void)qi;
    (void)ctx;
    (void)c;
    return InteExprn(b11, n + 2);
}
