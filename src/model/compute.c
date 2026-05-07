/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/cmeson.h>

#include <gemstore/math/matrix.h>
#include <gemstore/math/interplt.h>

#include <gemstore/param/argset.h>

#include <gemstore/types.h>
#include <gemstore/parse.h>
#include <gemstore/print.h>

#include <stdio.h>
#include <stdlib.h>

void compute_spectra_meson(const argsInput_t *input)
{
    int nmax = input->nmax;

    /* compute mass eigenvalues and eigenvectors */
    array_t eigenvalue = array_init(nmax);
    array_t rmsradius = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    argsGIModel_t args_model = argsGIModel_from(input);
    argsGIModelDy_t args_dynmc = {0};

    if (input->model == MODEL_GISTRING) {
        args_dynmc.model = MODEL_GISTRING;
    }
    else if (input->model == MODEL_GISCREEN) {
        args_dynmc.model = MODEL_GISCREEN;
    }

    /* Route to appropriate basis and system dispatcher */
    if (input->system == SYSTEM_MESON) {
        args_dynmc.system = SYSTEM_MESON;
        if (input->orbit == ORBIT_GEM) {
            spectra_meson_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
            radius_meson_GEM(input, &eigenvector, &rmsradius, nmax);
        }
        else if (input->orbit == ORBIT_CRG) {
            spectra_meson_CRG(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
            radius_meson_CRG(input, &eigenvector, &rmsradius, nmax);
        }
        else {
            fprintf(stderr, "Error: Unsupported orbit type for meson system.\n");
            exit(1);
        }
    }

    /* interpolate anomalies in RMS radius */
    print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);
    interpolate_quadratic(rmsradius.value, rmsradius.len);
    print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);

    /* write output into file */
    write_meson_spectra(input, &eigenvalue, &rmsradius, &eigenvector, nmax);

    /* choose to write potential or wavefunction */
    if (input->print_pot) {
        write_potential_GI(input, &args_model, &args_dynmc);
    }
    if (input->print_wfn) {
        write_meson_wfn(input, &eigenvector);
    }

    array_free(&eigenvalue);
    array_free(&rmsradius);
    matrix_free(&eigenvector);
}
