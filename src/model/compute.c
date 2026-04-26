/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/spectra.h>
#include <gemstore/model/radius.h>
#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/interplt.h>

#include <gemstore/types.h>
#include <gemstore/parse.h>
#include <gemstore/print.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void compute_spectra_meson(const argsInput_t *input)
{
    int nmax = input->nmax;

    /* compute mass eigenvalues and eigenvectors */
    array_t eigenvalue = array_init(nmax);
    array_t rmsradius = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    argsGIModel_t args_model = argsGIModel_from(input);
    argsGIModelDy_t args_dynmc = {0};

    if (input->param == PARAM_GISTRING_MESON) args_model = argsGIString_meson;
    else if (input->param == PARAM_GISCREEN_MESON) args_model = argsGIScreen_meson;
    else if (input->param == PARAM_GISCREEN_BBBAR) args_model = argsGIScreen_bbbar;
    else if (input->param == PARAM_GISCREEN_CCBAR) args_model = argsGIScreen_ccbar;
    else if (input->param == PARAM_GISTRING_CUSTOM) {
        parse_param_GISTRING(input->param_file, &args_model);
    }
    else if (input->param == PARAM_GISCREEN_CUSTOM) {
        parse_param_GISCREEN(input->param_file, &args_model);
    }
    else {
        array_free(&eigenvalue);
        array_free(&rmsradius);
        matrix_free(&eigenvector);
        return;
    }

    if (input->model == MODEL_GISTRING) {
        args_dynmc.model = MODEL_GISTRING;
        args_dynmc.system = SYSTEM_MESON;
        
        /* Route to appropriate basis dispatcher */
        if (input->orbit == ORBIT_GEM) {
            spectra_meson_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else if (input->orbit == ORBIT_CRG) {
            spectra_meson_CRG(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else {
            fprintf(stderr, "Unknown orbit basis type: %d\n", input->orbit);
        }
    }
    else if (input->model == MODEL_GISCREEN) {
        args_dynmc.model = MODEL_GISCREEN;
        args_dynmc.system = SYSTEM_MESON;
        
        /* Route to appropriate basis dispatcher */
        if (input->orbit == ORBIT_GEM) {
            spectra_meson_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else if (input->orbit == ORBIT_CRG) {
            spectra_meson_CRG(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else {
            fprintf(stderr, "Unknown orbit basis type: %d\n", input->orbit);
        }
    }
    else {
        array_free(&eigenvalue);
        array_free(&rmsradius);
        matrix_free(&eigenvector);
        return;
    }

    /* compute RMS radius */
    radius_meson_GEM(input, &eigenvector, &rmsradius, nmax);
    print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);

    /* interpolate anomalies in RMS radius */
    interpolate_quadratic(rmsradius.value, rmsradius.len);
    print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);
    
    write_meson_spectra(input, &eigenvalue, &rmsradius, &eigenvector, nmax);

    array_free(&eigenvalue);
    array_free(&rmsradius);
    matrix_free(&eigenvector);
}
