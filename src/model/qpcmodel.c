/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Quark-pair creation in the 3P0 model, one SHO function per meson.
 * The radial quantum number n starts at 1. beta is the SHO scale in GeV:
 * the momentum-space ground state is exp(-k^2 / (2 beta^2)).
 *
 * The partial-wave amplitude matches Ackleh, Barnes and Swanson,
 * Phys. Rev. D 54, 6811 (1996), up to one overall sign.
 * The constant sqrt(8) in the rate is part of that normalization.
 * Energies and the charge multiplicity belong to the decay solver.
 *
 * qpc_amplitude = rate(gamma) * mslp.
 * mslp sums the spin-orbit recoupling of the two quark-line diagrams.
 * epsilon_function is the spatial overlap of one diagram.
 * sum_d is the radial integral inside that overlap.
 * Quark masses enter only through the momentum sharing in coefficient_x.
 */

#include <gemstore/model/qpcmodel.h>
#include <gemstore/math/soc.h>

#include <math.h>
#include <stdlib.h>

/* Integer power. The exponent may be negative, and the base may be negative. */
static double ipow(double base, int n)
{
    double result = 1.0;
    int k = n;

    if (n < 0) {
        base = 1.0 / base;
        k = -n;
    }
    while (k > 0) {
        if (k & 1) result *= base;
        base *= base;
        k >>= 1;
    }
    return result;
}

/* (-1)^n for a non-negative integer. */
static double phase_int(int n)
{
    return (n & 1) ? -1.0 : 1.0;
}

/* n!. A negative argument drops the term. */
static double fact_n(int n)
{
    if (n < 0) return 0.0;
    return tgamma((double)n + 1.0);
}

/* 9j of rows (j1 j2 j12), (j3 j4 j34), (j13 j24 J),
 * times sqrt((2 j12+1)(2 j34+1)(2 j13+1)(2 j24+1)(2 J+1)). */
static double j9ch(double j1, double j2, double j12, double j3, double j4, double j34,
    double j13, double j24, double J)
{
    return sqrt(2.0 * j12 + 1.0) * sqrt(2.0 * j34 + 1.0) * sqrt(2.0 * j13 + 1.0)
        * sqrt(2.0 * j24 + 1.0) * sqrt(2.0 * J + 1.0)
        * nineJ_symbol(j1, j2, j12, j3, j4, j34, j13, j24, J);
}

/* Solid-harmonic factor that splits rank l into l1 and l-l1. */
static double coefficient_c(int l1, int l)
{
    return sqrt((4.0 * M_PI * fact_n(2 * l + 1)) / (fact_n(2 * l1 + 1) * fact_n(2 * (l - l1) + 1)));
}

/* Couple two spherical harmonics of ranks l1 and l2 to rank l at m = 0. */
static double coefficient_b(int l1, int l2, int l)
{
    return phase_int(l) / sqrt(4.0 * M_PI) * sqrt(2.0 * l1 + 1.0) * sqrt(2.0 * l2 + 1.0)
        * threeJ_symbol(l1, 0.0, l2, 0.0, l, 0.0);
}

/* Momentum fraction on the parent. Daughter shifts are
 * x - m1/(m1+m3) and x - m2/(m2+m3). m1, m2 are the parent quarks
 * and m3 is the created quark. */
static double coefficient_x(double m1, double m2, double m3, double ba, double bb, double bc)
{
    double num = (m1 / (m1 + m3) * bc * bc + m2 / (m2 + m3) * bb * bb) * ba * ba;
    double den = ba * ba * bb * bb + ba * ba * bc * bc + bb * bb * bc * bc;
    return num / den;
}

/* Angular integral of one monomial of degree px.
 * Even px sums ranks 2 sigma. Odd px sums ranks 2 sigma + 1. */
static double coefficient_i(int px, int l5, int l6, int l7, int l8, int L2)
{
    double s = sqrt(2.0 * l5 + 1.0) * sqrt(2.0 * l6 + 1.0) * sqrt(2.0 * l7 + 1.0) * sqrt(2.0 * l8 + 1.0);
    double res = 0.0;
    int sigma;
    int p;

    if ((px % 2) == 0) {
        p = px / 2;
        for (sigma = 0; sigma <= p; sigma++) {
            res += (ipow(4.0, sigma) * (4.0 * sigma + 1.0) * fact_n(p + sigma))
                / (fact_n(2 * p + 2 * sigma + 1) * fact_n(p - sigma))
                * threeJ_symbol(2.0 * sigma, 0.0, l5, 0.0, l7, 0.0)
                * threeJ_symbol(2.0 * sigma, 0.0, l6, 0.0, l8, 0.0)
                * sixJ_symbol(l5, l6, L2, l8, l7, 2.0 * sigma);
        }
        return phase_int(L2) * fact_n(2 * p) * s * res;
    }
    p = (px - 1) / 2;
    for (sigma = 0; sigma <= p; sigma++) {
        res += (ipow(4.0, sigma) * (4.0 * sigma + 3.0) * fact_n(p + sigma + 1))
            / (fact_n(2 * p + 2 * sigma + 3) * fact_n(p - sigma))
            * threeJ_symbol(2.0 * sigma + 1.0, 0.0, l5, 0.0, l7, 0.0)
            * threeJ_symbol(2.0 * sigma + 1.0, 0.0, l6, 0.0, l8, 0.0)
            * sixJ_symbol(l5, l6, L2, l8, l7, 2.0 * sigma + 1.0);
    }
    return 2.0 * phase_int(L2 + 1) * fact_n(2 * p + 1) * s * res;
}

/* Radial integral. Expand the three SHO functions and the momentum shifts,
 * then fold the leftover power of p into coefficient_i.
 * Gamma(z+1) at a non-positive integer is a pole and that term is dropped. */
static double sum_d(double dP, double dm1, double dm2, double dm3, int dna, int dnb, int dnc,
    int dla, int dlb, int dlc, double dba, double dbb, double dbc,
    int dl1, int dl2, int dl3, int dl4, int dl5, int dl6, int dl7, int dl8, int dlll)
{
    double da = coefficient_x(dm1, dm2, dm3, dba, dbb, dbc);
    double db = da - dm1 / (dm1 + dm3);
    double dc = da - dm2 / (dm2 + dm3);
    double product;
    double acc = 0.0;
    double cg = sqrt(0.5 * (1.0 / (dba * dba) + 1.0 / (dbb * dbb) + 1.0 / (dbc * dbc)));
    int di, dp, dl, dd, de, df, du, dv, dw;

    product = phase_int(dna + dnb + dnc - 3)
        * sqrt(8.0 * fact_n(dna - 1) * fact_n(dnb - 1) * fact_n(dnc - 1))
        * sqrt(tgamma(dna - 1 + dla + 1.5) * tgamma(dnb - 1 + dlb + 1.5) * tgamma(dnc - 1 + dlc + 1.5));
    for (di = 0; di < dna; di++) {
        for (dp = 0; dp < dnb; dp++) {
            for (dl = 0; dl < dnc; dl++) {
                for (dd = 0; dd <= di; dd++) {
                    for (de = 0; de <= dp; de++) {
                        for (df = 0; df <= dl; df++) {
                            for (du = 0; du <= di - dd; du++) {
                                for (dv = 0; dv <= dp - de; dv++) {
                                    for (dw = 0; dw <= dl - df; dw++) {
                                        int px = di + dp + dl - dd - de - df - du - dv - dw;
                                        int gpow = 2 * (du + dv + dw) + px - dl1 - dl2 - dl3 - dl4;
                                        double z = (dla + dlb + dlc + 1 - dl1 - dl2 - dl3 - dl4
                                            + 2 * (du + dv + dw) + px + 1) / 2.0;
                                        double term = phase_int(di + dp + dl)
                                            / (tgamma(dla + 1.5 + di) * tgamma(dlb + 1.5 + dp) * tgamma(dlc + 1.5 + dl));
                                        double fz;
                                        int twice;

                                        term /= ipow(dba, 2 * di) * ipow(dbb, 2 * dp) * ipow(dbc, 2 * dl);
                                        term /= fact_n(dna - 1 - di) * fact_n(dnb - 1 - dp) * fact_n(dnc - 1 - dl);
                                        term *= ipow(2.0, di - dd - du + dp - de - dv + dl - df - dw);
                                        term *= ipow(da, di + dd - du) * ipow(db, dp + de - dv) * ipow(dc, dl + df - dw);
                                        term /= fact_n(dd) * fact_n(du) * fact_n(di - dd - du);
                                        term /= fact_n(de) * fact_n(dv) * fact_n(dp - de - dv);
                                        term /= fact_n(df) * fact_n(dw) * fact_n(dl - df - dw);
                                        term *= coefficient_i(px, dl5, dl6, dl7, dl8, dlll);
                                        term *= ipow(dP, dl1 + dl2 + dl3 + dl4 + 2 * (dd + de + df) + px);
                                        term /= ipow(cg, gpow);
                                        /* Factorial[z] = Gamma(z+1). A non-positive integer is a pole and drops out. */
                                        twice = (int)lround(2.0 * z);
                                        if (fabs(2.0 * z - twice) > 1e-8) fz = tgamma(z + 1.0);
                                        else if (twice < 0 && (twice % 2) == 0) fz = 0.0;
                                        else fz = tgamma(z + 1.0);
                                        term *= fz;
                                        acc += term;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return product * acc;
}

/* Spatial overlap of one quark-line diagram.
 * The Gaussian is exp(-(p f)^2). The loops split L_B, L_C, the created
 * pair's L = 1, and L_A, couple those pieces, and call sum_d. */
static double epsilon_function(double efp, double efm1, double efm2, double efm3,
    int efna, int efnb, int efnc, int efla, int eflb, int eflc,
    double efba, double efbb, double efbc, int efl, int efll, int eflll)
{
    double x = coefficient_x(efm1, efm2, efm3, efba, efbb, efbc);
    double shift_b = x - efm1 / (efm1 + efm3);
    double shift_c = x - efm2 / (efm2 + efm3);
    double f = sqrt(0.5 * (x * x / (efba * efba) + shift_b * shift_b / (efbb * efbb) + shift_c * shift_c / (efbc * efbc)));
    double cg = sqrt(0.5 * (1.0 / (efba * efba) + 1.0 / (efbb * efbb) + 1.0 / (efbc * efbc)));
    double product = exp(-(efp * f) * (efp * f))
        / pow(cg, efla + eflb + eflc + 4.0)
        / 2.0
        / pow(efba, efla + 1.5) / pow(efbb, eflb + 1.5) / pow(efbc, eflc + 1.5);
    double acc = 0.0;
    int efll1, efll2, efll3, efll4, efll12, efll5, efll6, efll7, efll8;
    for (efll1 = 0; efll1 <= eflb; efll1++) {
        for (efll2 = 0; efll2 <= eflc; efll2++) {
            for (efll3 = 0; efll3 <= 1; efll3++) {
                for (efll4 = 0; efll4 <= efla; efll4++) {
                    int lo12 = abs(efll1 - efll2);
                    int hi12 = efll1 + efll2;
                    for (efll12 = lo12; efll12 <= hi12; efll12++) {
                        int lo5 = abs(efl - efll12);
                        int hi5 = efl + efll12;
                        for (efll5 = lo5; efll5 <= hi5; efll5++) {
                            int lo6 = abs((eflb - efll1) - (eflc - efll2));
                            int hi6 = abs((eflb - efll1) + (eflc - efll2));
                            for (efll6 = lo6; efll6 <= hi6; efll6++) {
                                int lo7 = abs(efll3 - efll4);
                                int hi7 = efll3 + efll4;
                                for (efll7 = lo7; efll7 <= hi7; efll7++) {
                                    int lo8 = abs((1 - efll3) - (efla - efll4));
                                    int hi8 = abs((1 - efll3) + (efla - efll4));
                                    for (efll8 = lo8; efll8 <= hi8; efll8++) {
                                        double cprod = coefficient_c(efll1, eflb) * coefficient_c(efll2, eflc)
                                            * coefficient_c(efll3, 1) * coefficient_c(efll4, efla);
                                        double b1, b2, b3, b4, b5, s6, j1, j2, term;
                                        if (cprod == 0.0) continue;
                                        b1 = coefficient_b(efll1, efll2, efll12);
                                        if (b1 == 0.0) continue;
                                        b2 = coefficient_b(efl, efll12, efll5);
                                        if (b2 == 0.0) continue;
                                        b3 = coefficient_b(eflb - efll1, eflc - efll2, efll6);
                                        if (b3 == 0.0) continue;
                                        b4 = coefficient_b(efll3, efll4, efll7);
                                        if (b4 == 0.0) continue;
                                        b5 = coefficient_b(1 - efll3, efla - efll4, efll8);
                                        if (b5 == 0.0) continue;
                                        s6 = sixJ_symbol(efl, efll12, efll5, efll6, eflll, efll);
                                        if (s6 == 0.0) continue;
                                        j1 = j9ch(efll1, eflb - efll1, eflb, efll2, eflc - efll2, eflc, efll12, efll6, efll);
                                        if (j1 == 0.0) continue;
                                        j2 = j9ch(efll3, 1 - efll3, 1, efll4, efla - efll4, efla, efll7, efll8, eflll);
                                        if (j2 == 0.0) continue;
                                        term = cprod * b1 * b2 * b3 * b4 * b5 * s6 * j1 * j2;
                                        term *= ipow(shift_b, efll1) * ipow(shift_c, efll2) * ipow(x - 1.0, efll3) * ipow(x, efll4);
                                        term *= phase_int(efll12 + efll6) * sqrt(2.0 * efll5 + 1.0);
                                        term *= sum_d(efp, efm1, efm2, efm3, efna, efnb, efnc, efla, eflb, eflc,
                                            efba, efbb, efbc, efll1, efll2, efll3, efll4, efll5, efll6, efll7, efll8, eflll);
                                        acc += term;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return product * acc;
}

/* Spin-space amplitude of outgoing spin ms and orbital ml, before the rate.
 * Sum the intermediate pair spin, the daughter orbital sum, and the orbital
 * that couples L_A with the created L = 1. A vanishing 6j or 9j is skipped.
 * Diagram 3-2-1-4 swaps the parent quark masses and multiplies by
 * (-1)^(1+S_A+S_B+S_C+L). Betas and daughter quantum numbers stay put. */
static double mslp(int ms, int ml, double mp, double mm1, double mm2, double mm3,
    int mn1, int mn2, int mn3, int msa, int msb, int msc, int mla, int mlb, int mlc,
    int mja, int mjb, int mjc, double mb1, double mb2, double mb3, double mf1432, double mf3214)
{
    double acc = 0.0;
    double pre = sqrt(8.0 / 3.0) * phase_int(1 + msb + ms + mjb + mjc);
    int mss_lo = abs(msa - 1) < abs(msb - msc) ? abs(msa - 1) : abs(msb - msc);
    int mss_hi = abs(msa + 1) > abs(msb + msc) ? abs(msa + 1) : abs(msb + msc);
    int mss, mll, mlll;

    for (mss = mss_lo; mss <= mss_hi; mss++) {
        for (mll = abs(mlb - mlc); mll <= mlb + mlc; mll++) {
            int mlll_lo = abs(mla - 1) < abs(mll - ml) ? abs(mla - 1) : abs(mll - ml);
            int mlll_hi = abs(mla + 1) > (mll + ml) ? abs(mla + 1) : (mll + ml);
            for (mlll = mlll_lo; mlll <= mlll_hi; mlll++) {
                double ang = sqrt(2.0 * mlll + 1.0)
                    * sixJ_symbol(msa, mla, mja, mlll, mss, 1.0)
                    * sixJ_symbol(mss, mll, ms, ml, mja, mlll);
                double rec;
                double term = 0.0;
                if (ang == 0.0) continue;
                rec = j9ch(msb, mlb, mjb, msc, mlc, mjc, mss, mll, ms)
                    * j9ch(0.5, 0.5, msa, 0.5, 0.5, 1.0, msb, msc, mss);
                if (rec == 0.0) continue;
                if (mf1432 != 0.0) {
                    term += mf1432 * epsilon_function(mp, mm1, mm2, mm3, mn1, mn2, mn3,
                        mla, mlb, mlc, mb1, mb2, mb3, ml, mll, mlll);
                }
                if (mf3214 != 0.0) {
                    term += mf3214 * phase_int(1 + msa + msb + msc + ml)
                        * epsilon_function(mp, mm2, mm1, mm3, mn1, mn2, mn3,
                            mla, mlb, mlc, mb1, mb2, mb3, ml, mll, mlll);
                }
                acc += pre * ang * phase_int(mss) * rec * term;
            }
        }
    }
    return acc;
}

double qpc_amplitude(const argsQPCModel_t *model, double p,
    const qpc_meson_t *a, const qpc_meson_t *b, const qpc_meson_t *c,
    int L, int S, double fl1, double fl2)
{
    /* 8 sqrt(pi) gamma / sqrt(8) / sqrt(2). sqrt(8) is the Ackleh constant. */
    double rate = (8.0 * sqrt(M_PI) * model->gamma) / sqrt(8.0) / sqrt(2.0);

    return rate * mslp(S, L, p, model->mq1, model->mq2, model->mq3,
        a->n, b->n, c->n,
        (int)lround(a->S), (int)lround(b->S), (int)lround(c->S),
        (int)lround(a->L), (int)lround(b->L), (int)lround(c->L),
        (int)lround(a->J), (int)lround(b->J), (int)lround(c->J),
        a->beta, b->beta, c->beta, fl1, fl2);
}
