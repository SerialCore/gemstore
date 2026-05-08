/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CBARYON_HELPERS
#define GEMSTORE_MODEL_CBARYON_HELPERS

#include <gemstore/basis/basis.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

typedef struct baryon_legacy_model {
    double mq;
    double ms;
    double mc;
    double mb;
    double b;
    double K;
    double alpha_ss;
    double alpha_so;
    double alpha_ten;
    double C;
    double Lambda;
} baryon_legacy_model_t;

baryon_legacy_model_t baryon_legacy_model_from(const argsGIModel_t *args_model);

void baryon_build_hamiltonian(const basis_list *qnlist_spfy,
    const basis_list *qnlist_full,
    const baryon_legacy_model_t *model,
    matrix_t *Hfi,
    matrix_t *Nfi);

#endif
