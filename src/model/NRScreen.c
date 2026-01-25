/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <model/NRScreen.h>

#include <math.h>

double NRScreen_T(double p, double m1, double m2)
{
    return m1 + m2 + p * p / (2.0 * m1) + p * p / (2.0 * m2);
}

double NRScreen_Vconf(double r, double alpha_s, double b, double mu, double c)
{
    if (r == 0.0) {
        return 0.0;
    }
    
    double Cij = -4.0 / 3.0;
    double coul = alpha_s / r;
    double screen = b * (1.0 - exp(-mu * r)) / mu;
    
    return Cij * coul - (3.0 / 4.0) * Cij * screen - (3.0 / 4.0) * Cij * c;
}

double NRScreen_Vcont(double r, double alpha_s, double sigma, double m1, double m2, double sds)
{
    double magnet = 32.0 * alpha_s * sds / (9.0 * m1 * m2);
    double factor = pow(sigma, 3.0) / sqrt(M_PI);
    double exp_term = exp(-sigma * sigma * r * r);

    return magnet * factor * exp_term;
}

double NRScreen_Vsocm(double r, double alpha_s, double m1, double m2, double lds1, double lds2)
{
    if (r == 0.0) {
        return 0.0;
    }
    
    double Cij = -4.0 / 3.0;
    double tensor = alpha_s / pow(r, 3.0);
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2) + (lds1 + lds2) / (m1 * m2));

    return -Cij * tensor * magnet;
}

double NRScreen_Vsotp(double r, double alpha_s, double b, double mu, double m1, double m2, double lds1, double lds2)
{
    if (r == 0.0) {
        return 0.0;
    }

    double Cij = -4.0 / 3.0;
    double dcoul_dr = -alpha_s / pow(r, 2.0);
    double dscreen_dr = b * exp(-mu * r);
    double dV_dr = Cij * dcoul_dr - (3.0 / 4.0) * Cij * dscreen_dr;
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2));
    
    return -0.5 / r * dV_dr * magnet;
}

double NRScreen_Vtens(double r, double alpha_s, double m1, double m2, double tens)
{
    if (r == 0.0) {
        return 0.0;
    }

    double tensor = alpha_s / pow(r, 3.0);
    double magnet = 4.0 * tens / (3.0 * m1 * m2);

    return tensor * magnet;
}