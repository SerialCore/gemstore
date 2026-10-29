/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/bagmodel.h>
#include <gemstore/model/cmeson.h>
#include <gemstore/model/cbaryon.h>
#include <gemstore/model/cmol2.h>

#include <gemstore/math/matrix.h>
#include <gemstore/math/interplt.h>

#include <gemstore/param/argset.h>

#include <gemstore/input.h>
#include <gemstore/types.h>
#include <gemstore/print.h>
#include <gemstore/write.h>

#include <stdio.h>
#include <stdlib.h>

void compute_spectra_bag(const argsInput_t *input)
{
    argsBAGModel_t param;
    bag_mit_t state[CMI_COLOR_MAX * CMI_SPIN_MAX];
    char order[7] = {0};
    int nstate = 0;

    param = argsBAGModel_from(input);
    bag_compute(input, &param, state, &nstate, order);
    print_bag_spectra(state, nstate, order[0] ? order : NULL);
    write_bag_spectra(input, state, nstate, order[0] ? order : NULL);
}

void compute_spectra_meson(const argsInput_t *input)
{
    int nmax = input->nmax;

    array_t eigenvalue = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    matrix_t overlap = {0};
    array_t rms12 = {0};

    if (input->model == MODEL_NRSTRING || input->model == MODEL_NRSCREEN) {
        argsNRModel_t args_model = argsNRModel_from(input);
        argsNRModelDy_t args_dynmc = {0};
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->basis == BASIS_GEM) {
            meson_nrmodel_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12);
        }
        else if (input->basis == BASIS_SHO) {
            meson_nrmodel_SHO(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12);
        }
        else {
            fprintf(stderr, "Error: Unsupported orbit type for meson system.\n");
            exit(1);
        }

        interpolate_divergence(rms12.value, rms12.len);
        print_meson_spectra(&eigenvalue, &eigenvector, &overlap, &rms12, eigenvalue.len);
        write_meson_spectra(input, &eigenvalue, &rms12, &eigenvector, eigenvalue.len);

        if (input->print_pot) {
            write_meson_nr_pot(input, &args_model, &args_dynmc);
        }
        if (input->print_wfn) {
            write_meson_wfn(input, &eigenvector);
        }

        array_free(&eigenvalue);
        array_free(&rms12);
        matrix_free(&eigenvector);
        matrix_free(&overlap);
    }
    else if (input->model == MODEL_GISTRING || input->model == MODEL_GISCREEN) {
        argsGIModel_t args_model = argsGIModel_from(input);
        argsGIModelDy_t args_dynmc = {0};
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->basis == BASIS_GEM) {
            meson_gimodel_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12);
        }
        else if (input->basis == BASIS_SHO) {
            meson_gimodel_SHO(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12);
        }
        else {
            fprintf(stderr, "Error: Unsupported orbit type for meson system.\n");
            exit(1);
        }

        interpolate_divergence(rms12.value, rms12.len);
        print_meson_spectra(&eigenvalue, &eigenvector, &overlap, &rms12, eigenvalue.len);
        write_meson_spectra(input, &eigenvalue, &rms12, &eigenvector, eigenvalue.len);

        if (input->print_pot) {
            write_meson_gi_pot(input, &args_model, &args_dynmc);
        }
        if (input->print_wfn) {
            write_meson_wfn(input, &eigenvector);
        }

        array_free(&eigenvalue);
        array_free(&rms12);
        matrix_free(&eigenvector);
        matrix_free(&overlap);
    }
    else {
        fprintf(stderr, "Error: Unsupported model for meson spectra.\n");
        exit(1);
    }
}

void compute_spectra_baryon(const argsInput_t *input)
{
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
        argsGIModel_t args_model = argsGIModel_from(input);
        argsGIModelDy_t args_dynmc = {0};
        args_dynmc.model = input->model;
        args_dynmc.system = input->system;

        if (input->basis != BASIS_GEM) {
            fprintf(stderr, "Error: Baryon SPECTRA currently supports GEM basis only.\n");
            exit(1);
        }
        baryon_gimodel_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12, &rms13, &rms23);

        print_baryon_spectra(&eigenvalue, &eigenvector, &overlap, &rms12, &rms13, &rms23, eigenvalue.len);
        write_baryon_spectra(input, &eigenvalue, &eigenvector, &rms12, &rms13, &rms23, eigenvalue.len);

        if (input->print_pot) {
            write_baryon_gi_pot(input, &args_model, &args_dynmc);
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

void compute_spectra_mol2(const argsInput_t *input)
{
    int nmax = input->nmax;
    int nch = mol2_channel_count(input);
    int ntot = nch * nmax;
    array_t eigenvalue;
    matrix_t eigenvector;
    matrix_t overlap = {0};
    array_t rms12 = {0};
    argsOBEModel_t args_model;
    argsOBEModelDy_t args_dynmc = {0};

    if (input->model != MODEL_OBESTAT) {
        fprintf(stderr, "Error: MOLECUL2 SPECTRA supports OBESTAT only.\n");
        exit(1);
    }
    if (input->basis != BASIS_GEM) {
        fprintf(stderr, "Error: MOLECUL2 SPECTRA currently supports GEM basis only.\n");
        exit(1);
    }
    if (input->m1 <= 0.0 || input->m2 <= 0.0) {
        fprintf(stderr, "Error: MOLECUL2 needs positive physical masses m1 and m2.\n");
        exit(1);
    }

    if (nch < 1 || ntot < 1) {
        fprintf(stderr, "Error: MOLECUL2 has no partial-wave channels for this J, L.\n");
        exit(1);
    }

    eigenvalue = array_init(ntot);
    eigenvector = matrix_init(ntot, ntot);

    args_model = argsOBEModel_from(input);
    mol2_obestat_GEM(input, &args_model, &args_dynmc, &eigenvalue, &eigenvector, &overlap, &rms12);

    print_mol2_spectra(&eigenvalue, &eigenvector, &overlap, &rms12, input->m1 + input->m2, eigenvalue.len);
    write_meson_spectra(input, &eigenvalue, &rms12, &eigenvector, eigenvalue.len);

    if (input->print_pot) {
        write_obe_pot(input, &args_model);
    }
    if (input->print_wfn) {
        write_mol2_wfn(input, &eigenvector);
    }

    array_free(&eigenvalue);
    array_free(&rms12);
    matrix_free(&eigenvector);
    matrix_free(&overlap);
}
