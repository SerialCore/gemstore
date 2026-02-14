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
double NRScreen_T(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Confining potential for NRScreen */
double NRScreen_Vconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Contact potential for NRScreen */
double NRScreen_Vcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for NRScreen */
double NRScreen_Vsocm(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for NRScreen */
double NRScreen_Vsotp(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Tenser potential for NRScreen */
double NRScreen_Vtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

#endif