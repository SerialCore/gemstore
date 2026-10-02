/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/cmeson.h>
#include <gemstore/model/cbaryon.h>
#include <gemstore/model/wfntrans.h>

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

    array_t eigenvalue = array_init(nmax);
    array_t rmsradius = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);

    if (input->model == MODEL_NRSTRING || input->model == MODEL_NRSCREEN) {
        argsNRModel_t args_model = argsNRModel_from(input);
        argsNRModelDy_t args_dynmc = {0};
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->orbit == ORBIT_GEM) {
            meson_nrmodel_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else if (input->orbit == ORBIT_SHO) {
            meson_nrmodel_SHO(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else {
            fprintf(stderr, "Error: Unsupported orbit type for meson system.\n");
            exit(1);
        }

        radius_meson_rms(input, &eigenvector, &rmsradius, nmax);
        interpolate_divergence(rmsradius.value, rmsradius.len);
        print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);
        write_meson_spectra(input, &eigenvalue, &rmsradius, &eigenvector, nmax);

        if (input->print_pot) {
            write_meson_nr_pot(input, &args_model, &args_dynmc);
        }
        if (input->print_wfn) {
            write_meson_wfn(input, &eigenvector);
        }

        array_free(&eigenvalue);
        array_free(&rmsradius);
        matrix_free(&eigenvector);
    }
    else if (input->model == MODEL_GISTRING || input->model == MODEL_GISCREEN) {
        argsGIModel_t args_model = argsGIModel_from(input);
        argsGIModelDy_t args_dynmc = {0};
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->orbit == ORBIT_GEM) {
            meson_gimodel_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else if (input->orbit == ORBIT_SHO) {
            meson_gimodel_SHO(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, nmax);
        }
        else {
            fprintf(stderr, "Error: Unsupported orbit type for meson system.\n");
            exit(1);
        }

        radius_meson_rms(input, &eigenvector, &rmsradius, nmax);
        interpolate_divergence(rmsradius.value, rmsradius.len);
        print_meson_spectra(&eigenvalue, &rmsradius, &eigenvector, nmax);
        write_meson_spectra(input, &eigenvalue, &rmsradius, &eigenvector, nmax);

        if (input->print_pot) {
            write_meson_pot(input, &args_model, &args_dynmc);
        }
        if (input->print_wfn) {
            write_meson_wfn(input, &eigenvector);
        }

        array_free(&eigenvalue);
        array_free(&rmsradius);
        matrix_free(&eigenvector);
    }
    else {
        fprintf(stderr, "Error: Unsupported model for meson spectra.\n");
        exit(1);
    }
}

void compute_spectra_baryon(const argsInput_t *input)
{
    argsGIModel_t args_model = argsGIModel_from(input);
    argsGIModelDy_t args_dynmc = {0};
    array_t eigenvalue = {0};
    matrix_t eigenvector = {0};
    matrix_t overlap = {0};
    array_t rms12 = {0};
    array_t rms13 = {0};
    array_t rms23 = {0};

    if (input->model == MODEL_NRSTRING || input->model == MODEL_NRSCREEN) {
        fprintf(stderr, "Error: Non-relativistic baryon SPECTRA is not implemented.\n");
        exit(1);
    }
    else if (input->model == MODEL_GISTRING || input->model == MODEL_GISCREEN) {
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->orbit != ORBIT_GEM) {
            fprintf(stderr, "Error: Baryon SPECTRA currently supports GEM basis only.\n");
            exit(1);
        }
        spectra_baryon_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12, &rms13, &rms23);

        print_baryon_spectra(&eigenvalue, &eigenvector, &overlap, &rms12, &rms13, &rms23, eigenvalue.len);
        write_baryon_spectra(input, &eigenvalue, &eigenvector, &rms12, &rms13, &rms23, eigenvalue.len);

        if (input->print_pot) {
            write_baryon_pot(input, &args_model, &args_dynmc);
        }

        array_free(&eigenvalue);
        array_free(&rms12);
        array_free(&rms13);
        array_free(&rms23);
        matrix_free(&eigenvector);
        matrix_free(&overlap);
    }
    else {
        fprintf(stderr, "Error: Unsupported model for baryon spectra.\n");
        exit(1);
    }
}
