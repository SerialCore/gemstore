/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Elastic static OBE pair potential. Radial kernels are the n = 1, 2
 * monopole Fourier pieces in OBELib.wl (poleCentral, poleSS, poleT).
 * Operator coefficients come from OBEPotential with StaticBaryon, classified
 * by the two leg types rather than by a named hadron system:
 *   Dirac–Dirac           BaryonBaryon
 *   meson–meson           MesonMeson
 *   one Dirac, one meson  MesonBaryon
 * Elastic 0-, 0+ and 1- meson legs are covered (0+ uses the 0- vertices).
 * P ↔ P*, AxialVector, finite-mass Gordon, dipole / Λ1 ≠ Λ2, and extra
 * bosons are not.
 * π, σ, ρ, ω masses and G-parity live in obe_boson[] below.
 */

#include <gemstore/model/obemodel.h>

#include <math.h>
#include <stddef.h>

#define OBE_UNPACK(ctx) \
    const argsOBEStatModel_t *args_model = ((const obe_pot_ctx_t *)(ctx))->model; \
    const argsOBEStatModelDy_t *args_dynmc = ((const obe_pot_ctx_t *)(ctx))->dyn; \
    (void)args_model; \
    (void)args_dynmc

enum {
    OBE_EX_SCALAR = 0,
    OBE_EX_PSEUDOSCALAR,
    OBE_EX_VECTOR
};

enum {
    OBE_PI = 0,
    OBE_SIGMA,
    OBE_RHO,
    OBE_OMEGA,
    OBE_NBOSON
};

/* Liu/OBE.nb values. Not input. */
#define OBE_MPI  0.138
#define OBE_FPI  0.132
#define OBE_MSIG 0.600
#define OBE_MRHO 0.770
#define OBE_MOMG 0.780

typedef struct obe_boson {
    int type;
    double mass;
    double heavy_sign;
    double gparity;
    int isovector;
} obe_boson_t;

static const obe_boson_t obe_boson[OBE_NBOSON] = {
    { OBE_EX_PSEUDOSCALAR, OBE_MPI,  1.0, -1.0, 1 }, /* π: HeavySign unused */
    { OBE_EX_SCALAR,       OBE_MSIG, 1.0,  1.0, 0 },
    { OBE_EX_VECTOR,       OBE_MRHO, 1.0,  1.0, 1 }, /* ρ */
    { OBE_EX_VECTOR,       OBE_MOMG, -1.0, -1.0, 0 }  /* ω */
};

/* Dirac and 1- carry spin; elastic 0- does not. */
static int obe_spinning(obe_leg_type_t type)
{
    return type == OBE_LEG_DIRAC || type == OBE_LEG_VECTOR;
}

/*
 * OBELib poleCentral[2], poleSS[1,2], poleT[1,2]: Fourier transforms of
 * 1/(q²+μ²)^n for central, (qi qj) spin-spin, and tensor. n = 1 is the
 * meson pole, n = 2 is the double pole at the cutoff.
 */
static double obe_pole_c2(double mu, double r)
{
    return exp(-mu * r) / (8.0 * M_PI * mu);
}

static double obe_pole_ss1(double mu, double r)
{
    return -mu * mu * exp(-mu * r) / (12.0 * M_PI * r);
}

static double obe_pole_ss2(double mu, double r)
{
    return (2.0 - mu * r) * exp(-mu * r) / (24.0 * M_PI * r);
}

static double obe_pole_t1(double mu, double r)
{
    double mur = mu * r;
    return -(3.0 + 3.0 * mur + mur * mur) * exp(-mur) / (12.0 * M_PI * r * r * r);
}

static double obe_pole_t2(double mu, double r)
{
    return -(1.0 + mu * r) * exp(-mu * r) / (24.0 * M_PI * r);
}

/*
 * Coordinate-space kernel of F_Λ(q)² / (q² + m²) with a common monopole
 * F_Λ = (Λ² − m²)/(Λ² + q²) on each vertex. Partial fractions:
 *   1/(q²+m²) − 1/(q²+Λ²) − (Λ²−m²)/(q²+Λ²)².
 * The (e^{-mr} − e^{-Λr})/r piece is written with sinh so r → 0 is stable.
 * r ≤ 0 returns the analytic value (Λ − m)² / (8πΛ).
 */
static double obe_kern_c(double m, double lam, double r)
{
    double dm, z, sinhc, yuk, sub;

    if (m <= 0.0 || lam <= 0.0) return 0.0;
    if (r <= 0.0) return (lam - m) * (lam - m) / (8.0 * M_PI * lam);

    dm = lam - m;
    z = 0.5 * dm * r;
    if (fabs(z) < 1e-8) sinhc = 1.0 + z * z / 6.0;
    else sinhc = sinh(z) / z;
    yuk = exp(-0.5 * (m + lam) * r) * dm * sinhc / (4.0 * M_PI);
    sub = (lam * lam - m * m) * exp(-lam * r) / (8.0 * M_PI * lam);
    return yuk - sub;
}

/*
 * Same partial fractions through poleSS: kernel of fQ = F_Λ² / (q² + m²)
 * in the spin-spin channel. Small-r series avoids 1/r cancellation.
 * r ≤ 0 returns (Λ − m)² (Λ + 2m) / (24π).
 */
static double obe_kern_ss(double m, double lam, double r)
{
    double d2, num;

    if (m <= 0.0 || lam <= 0.0) return 0.0;
    if (r <= 0.0) return (lam - m) * (lam - m) * (lam + 2.0 * m) / (24.0 * M_PI);

    if (r * lam < 1e-3) {
        double l2 = lam * lam, l3 = l2 * lam, l4 = l2 * l2, l5 = l4 * lam;
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, m5 = m4 * m;
        num = -36.0 * lam * m2
            + 24.0 * l2 * m2 * r
            + 6.0 * l5 * r * r
            - 12.0 * l4 * r
            + 12.0 * l3
            - 10.0 * l3 * m2 * r * r
            + 24.0 * m3
            - 12.0 * m4 * r
            + 4.0 * m5 * r * r;
        return num / (288.0 * M_PI);
    }

    d2 = lam * lam - m * m;
    return obe_pole_ss1(m, r) - obe_pole_ss1(lam, r) - d2 * obe_pole_ss2(lam, r);
}

/*
 * Same partial fractions through poleT. The kernel vanishes at r = 0.
 * The 1/r³ poles cancel only to O(r); use the series for rΛ < 10^{-3}.
 */
static double obe_kern_t(double m, double lam, double r)
{
    double d2, num, m2, m4, l2, l3, l4, l5;

    if (m <= 0.0 || lam <= 0.0 || r <= 0.0) return 0.0;

    if (r * lam < 1e-3) {
        m2 = m * m;
        m4 = m2 * m2;
        l2 = lam * lam;
        l3 = l2 * lam;
        l4 = l2 * l2;
        l5 = l4 * lam;
        num = r * (60.0 * l2 * m2
            + 24.0 * l5 * r
            - 40.0 * l3 * m2 * r
            + m4 * (-30.0 + 16.0 * m * r)
            + 15.0 * l4 * (-2.0 + m2 * r * r));
        return num / (2880.0 * M_PI);
    }

    d2 = lam * lam - m * m;
    return obe_pole_t1(m, r) - obe_pole_t1(lam, r) - d2 * obe_pole_t2(lam, r);
}

/*
 * Fourier transform of q² F_Λ² / (q² + m²) = F_Λ² − m² F_Λ² / (q² + m²).
 * Needed for the magnetic vector spin-spin piece (fSS ~ q²), together with
 * obe_kern_ss acting on fQ.
 */
static double obe_kern_q2c(double m, double lam, double r)
{
    double d2, pc2;

    if (m <= 0.0 || lam <= 0.0) return 0.0;
    d2 = lam * lam - m * m;
    pc2 = (r <= 0.0) ? 1.0 / (8.0 * M_PI * lam) : obe_pole_c2(lam, r);
    return d2 * d2 * pc2 - m * m * obe_kern_c(m, lam, r);
}

/* Liu table g → OBE vertex Coupling / PseudovectorCoupling. */
static double obe_pi_vertex(obe_leg_type_t type, double g_liu)
{
    double den = sqrt(2.0) * OBE_FPI;
    if (type == OBE_LEG_DIRAC) return g_liu * OBE_MPI / den;
    return g_liu / den;
}

static void obe_raw_gf(const argsOBEStatModel_t *model, int boson, int which,
    double *g, double *f)
{
    int first = (which == 1);
    obe_leg_type_t type = first ? model->type1 : model->type2;

    *f = 0.0;
    switch (boson) {
        case OBE_PI:
            *g = obe_pi_vertex(type, first ? model->g_pi1 : model->g_pi2);
            break;
        case OBE_SIGMA:
            *g = first ? model->g_sigma1 : model->g_sigma2;
            break;
        case OBE_RHO:
            *g = first ? model->g_rho1 : model->g_rho2;
            *f = first ? model->f_rho1 : model->f_rho2;
            break;
        case OBE_OMEGA:
            *g = first ? model->g_omega1 : model->g_omega2;
            *f = first ? model->f_omega1 : model->f_omega2;
            break;
        default:
            *g = 0.0;
            break;
    }
}

/*
 * Apply OBELib legSign to one vertex: HeavySign on a meson vector coupling,
 * and −1 when ChargeConjugation is set on a G = −1 exchange (π, ω). Dirac
 * legs are left unchanged. which is 1 or 2.
 */
static void obe_leg_coup(const argsOBEStatModel_t *model, int boson, int which,
    double *g, double *f)
{
    const obe_boson_t *b = &obe_boson[boson];
    obe_leg_type_t type = (which == 1) ? model->type1 : model->type2;
    int cc = (which == 1) ? model->cc1 : model->cc2;
    double s = 1.0;

    obe_raw_gf(model, boson, which, g, f);
    if (type != OBE_LEG_DIRAC && b->type == OBE_EX_VECTOR) s *= b->heavy_sign;
    if (type != OBE_LEG_DIRAC && cc && b->gparity < 0.0) s *= -1.0;
    *g *= s;
    *f *= s;
}

static double obe_iso(const argsOBEStatModel_t *model, int boson)
{
    double iso;

    if (!obe_boson[boson].isovector) return 1.0;
    iso = (boson == OBE_PI) ? model->iso_pi : model->iso_rho;
    return (iso == 0.0) ? 1.0 : iso;
}

/*
 * One boson. Coefficients match OBEPotential in the static elastic limit.
 *
 * Scalar:     Vc = −g1 g2 Y                    for every leg pair
 * Vector g g: Vc = +g1 g2 Y                    electric time-component; 0- included
 * Vector f f: Vss, Vt from β = −f1 f2 / (4M²)  static magneton; no (g+f)
 * Pion:       Vss = Vt = κ × kernel, κ below;  0- elastic pion is zero
 *
 * Pion κ uses the PV vertex on a Dirac leg (divide by m_ex) and the
 * superfield vertex on a meson. The overall sign is −1 for BB and MM
 * (zero or two Dirac legs) and +1 for MB (exactly one Dirac leg), which
 * is the contraction of {0, g q} against {0, −(f/m) q_leg} with q2 = −q.
 * eta multiplies MM pion and MM magnetic vector only.
 */
static void obe_add_exchange(const argsOBEStatModel_t *model, int boson,
    double r, double *vc, double *vss, double *vt)
{
    const obe_boson_t *b = &obe_boson[boson];
    double g1, g2, f1, f2, m, lam, iso, yc, kss, kt, eta;
    int t1 = model->type1;
    int t2 = model->type2;
    int mm = (t1 != OBE_LEG_DIRAC && t2 != OBE_LEG_DIRAC);
    int ndirac;
    double p1, p2, kappa, beta;

    obe_leg_coup(model, boson, 1, &g1, &f1);
    obe_leg_coup(model, boson, 2, &g2, &f2);

    m = b->mass;
    lam = model->lambda;
    iso = obe_iso(model, boson);
    eta = mm ? model->eta : 1.0;
    if (eta == 0.0) eta = 1.0;

    yc = obe_kern_c(m, lam, r);
    kss = obe_kern_ss(m, lam, r);
    kt = obe_kern_t(m, lam, r);

    if (b->type == OBE_EX_SCALAR) {
        *vc += -g1 * g2 * iso * yc;
        return;
    }

    if (b->type == OBE_EX_PSEUDOSCALAR) {
        if (!obe_spinning(t1) || !obe_spinning(t2)) return;
        ndirac = (t1 == OBE_LEG_DIRAC) + (t2 == OBE_LEG_DIRAC);
        p1 = (t1 == OBE_LEG_DIRAC) ? g1 / m : g1;
        p2 = (t2 == OBE_LEG_DIRAC) ? g2 / m : g2;
        kappa = ((ndirac == 1) ? 1.0 : -1.0) * p1 * p2 * iso * eta;
        *vss += kappa * kss;
        *vt += kappa * kt;
        return;
    }

    /* Static vector: electric g from v0 only; spatial nrGi is zero so the
     * magneton is f/(2M) with no Gordon g. */
    *vc += g1 * g2 * iso * yc;
    if (!obe_spinning(t1) || !obe_spinning(t2)) return;
    if (model->tensor_mass1 <= 0.0 || model->tensor_mass2 <= 0.0) return;
    if (f1 == 0.0 || f2 == 0.0) return;

    /* OBELib: (f1/(2 M1))(f2/(2 M2)); M1 and M2 are per-leg TensorMass. */
    beta = -f1 * f2 / (4.0 * model->tensor_mass1 * model->tensor_mass2) * iso * eta;
    *vss += beta * (obe_kern_q2c(m, lam, r) - kss);
    *vt += -beta * kt;
}

void obe_radial(const argsOBEStatModel_t *model, double r, double *vc, double *vss, double *vt)
{
    int i;
    double c = 0.0, ss = 0.0, t = 0.0;

    if (model != NULL && model->lambda > 0.0) {
        for (i = 0; i < OBE_NBOSON; i++) {
            obe_add_exchange(model, i, r, &c, &ss, &t);
        }
    }
    if (vc) *vc = c;
    if (vss) *vss = ss;
    if (vt) *vt = t;
}

/* Fold a spin–orbit matrix element into a radial kernel. NULL dyn leaves the kernel. */
static double obe_op(const argsOBEStatModelDy_t *dyn, double op, double radial)
{
    if (dyn == NULL) return radial;
    return op * radial;
}

double OBEVc(double r, void *ctx)
{
    double vc;
    OBE_UNPACK(ctx);
    obe_radial(args_model, r, &vc, NULL, NULL);
    return obe_op(args_dynmc, args_dynmc ? args_dynmc->OCent : 1.0, vc);
}

double OBEVss(double r, void *ctx)
{
    double vss;
    OBE_UNPACK(ctx);
    obe_radial(args_model, r, NULL, &vss, NULL);
    return obe_op(args_dynmc, args_dynmc ? args_dynmc->OSdS : 1.0, vss);
}

double OBEVtens(double r, void *ctx)
{
    double vt;
    OBE_UNPACK(ctx);
    obe_radial(args_model, r, NULL, NULL, &vt);
    return obe_op(args_dynmc, args_dynmc ? args_dynmc->OTens : 1.0, vt);
}

double OBEVls1(double r, void *ctx)
{
    (void)r;
    (void)ctx;
    return 0.0;
}

double OBEVls2(double r, void *ctx)
{
    (void)r;
    (void)ctx;
    return 0.0;
}
