/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/compute.h>
#include <gemstore/model/cbaryon.h>
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

    /* fix anomalies in mass and RMS radius */
    interpolate_fix_divergence(rmsradius.value, rmsradius.len);
    interpolate_fix_divergence(rmsradius.value, rmsradius.len);
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

void compute_spectra_baryon(const argsInput_t *input)
{
    int nmax = input->nmax * input->nmax * 3;
    int basis_len = 0;
    int out_len;

    if (input->orbit != ORBIT_GEM) {
        fprintf(stderr, "Error: baryon spectra currently supports only GEM basis.\n");
        exit(1);
    }

    array_t eigenvalue = array_init(nmax);
    matrix_t eigenvector = matrix_init(nmax, nmax);
    matrix_t overlap = matrix_init(nmax, nmax);

    spectra_baryon_GEM(input, &eigenvalue, &eigenvector, &overlap, nmax, &basis_len);
    out_len = (eigenvalue.len < input->nmax) ? eigenvalue.len : input->nmax;

    printf("\n");
    printf("================================================================================\n");
    printf("                         BARYON BASIS SUMMARY (GEM)                            \n");
    printf("================================================================================\n");
    printf("\n");
    printf("nmax per Jacobi axis:   %d\n", input->nmax);
    printf("Jacobi truncation Lmax: %d\n", input->Lmax);
    printf("Requested J^P:          %.1f^%c\n", input->J, (input->P > 0) ? '+' : '-');
    printf("Pair color factor:      -2/3 for (12), (13), (23)\n");
    printf("Solver scope:           native direct integration, S-wave only\n");
    printf("Internal matrix size:   %d x %d\n", basis_len, basis_len);
    printf("Retained subspace:      %d\n", eigenvalue.len);
    printf("Reported states:        %d\n", out_len);
    printf("\n");

    print_baryon_spectra(&eigenvalue, &eigenvector, &overlap, out_len);
    write_baryon_spectra(input, &eigenvalue, &eigenvector, out_len);

    if (input->print_pot || input->print_wfn) {
        fprintf(stderr, "Warning: baryon potential and wavefunction export are not implemented yet for the native solver.\n");
    }

    array_free(&eigenvalue);
    matrix_free(&eigenvector);
    matrix_free(&overlap);
}
