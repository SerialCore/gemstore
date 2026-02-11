/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_NRSCREEN
#define GEMSTORE_MODEL_NRSCREEN

#include <gemstore/model/model.h>

/* Default meson parameters for model NRScreen */
extern const argsModel_t argsNRScreen_meson;

/* Kinetic energy for NRScreen */
double NRScreen_T(double p, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

/* Confining potential for NRScreen */
double NRScreen_Vconf(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

/* Contact potential for NRScreen */
double NRScreen_Vcont(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

/* Spin-orbit coulping for NRScreen */
double NRScreen_Vsocm(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

/* Thomas precession for NRScreen */
double NRScreen_Vsotp(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

/* Tenser potential for NRScreen */
double NRScreen_Vtens(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

#endif