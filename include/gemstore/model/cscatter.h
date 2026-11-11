/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CSCATTER
#define GEMSTORE_MODEL_CSCATTER

#include <gemstore/input.h>
#include <gemstore/model/eftmodel.h>
#include <gemstore/model/halmodel.h>
#include <gemstore/model/obemodel.h>
#include <gemstore/model/ykwmodel.h>
#include <gemstore/types.h>

/* Projected two-body potential passed to the effective-range solver.
 * The context pointers stay valid for the call. */
typedef struct scatter_job {
    model_type_t model;
    eft_pot_ctx_t eft;
    hal_pot_ctx_t hal;
    obe_pot_ctx_t obe;
    ykw_pot_ctx_t ykw;
} scatter_job_t;

/* Two-body S-wave scattering length for one MOLECUL2 input.
 * Projects the diagonal l=0 potential, then calls the effective-range solver.
 * On success writes the reduced mass in GeV and a0, r0 in fm. */
void scatter_ere_mol2(const argsInput_t *input, double *mu_gev, double *a0_fm, double *r0_fm);

#endif
