/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/write.h>
#include <gemstore/math/soc.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/model/nrmodel.h>
#include <gemstore/model/wfntrans.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>

int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len)
{
    if (input == NULL || mass == NULL || radius == NULL || vector == NULL || len <= 0) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_spectra()\n");
        return 0;
    }

    int nmax = input->nmax;
    FILE *pf;
    cJSON *root = NULL;
    cJSON *states = NULL;
    char *json_text = NULL;
    int state = 0;

    char path[267];
    sprintf(path, "%s%s", input->project, ".state.json");
    pf = fopen(path, "w");

    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    /* Create root object with "states" array (no other metadata) */
    root = cJSON_CreateObject();
    if (root == NULL) {
        fprintf(stderr, "Error: Cannot allocate JSON root for %s\n", path);
        fclose(pf);
        return 0;
    }

    states = cJSON_AddArrayToObject(root, "states");
    if (states == NULL) {
        fprintf(stderr, "Error: Cannot allocate states array for %s\n", path);
        cJSON_Delete(root);
        fclose(pf);
        return 0;
    }

    /* Populate states array with state objects */
    for (int n = 0; n < len; n++) {
        cJSON *state_obj = cJSON_CreateObject();
        cJSON *eigenvector = cJSON_CreateArray();

        if (state_obj == NULL || eigenvector == NULL) {
            cJSON_Delete(eigenvector);
            cJSON_Delete(state_obj);
            fprintf(stderr, "Error: Cannot allocate state JSON for %s\n", path);
            cJSON_Delete(root);
            fclose(pf);
            return 0;
        }

        cJSON_AddNumberToObject(state_obj, "index", n + 1);
        cJSON_AddNumberToObject(state_obj, "mass", mass->value[n]);
        cJSON_AddNumberToObject(state_obj, "rms_radius", radius->value[n]);

        for (int d = 0; d < nmax; d++) {
            cJSON_AddItemToArray(eigenvector, cJSON_CreateNumber(vector->value[n][d]));
        }

        cJSON_AddItemToObject(state_obj, "eigenvector", eigenvector);
        cJSON_AddItemToArray(states, state_obj);
    }

    json_text = cJSON_Print(root);
    if (json_text == NULL) {
        fprintf(stderr, "Error: Cannot serialize JSON output for %s\n", path);
        cJSON_Delete(root);
        fclose(pf);
        return 0;
    }

    if (fputs(json_text, pf) == EOF) {
        fprintf(stderr, "Error: Failed writing JSON output to %s\n", path);
        cJSON_Delete(root);
        free(json_text);
        fclose(pf);
        return 0;
    }

    state = fclose(pf) == 0 ? 1 : 0;
    cJSON_Delete(root);
    free(json_text);

    return state;
}

int write_baryon_spectra(const argsInput_t *input, const array_t *mass, const matrix_t *vector,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, int len)
{
    if (input == NULL || mass == NULL || vector == NULL || len <= 0) {
        fprintf(stderr, "Error: Invalid input parameters to write_baryon_spectra()\n");
        return 0;
    }

    FILE *pf;
    cJSON *root = NULL;
    cJSON *states = NULL;
    char *json_text = NULL;
    int state = 0;

    char path[267];
    sprintf(path, "%s%s", input->project, ".state.json");
    pf = fopen(path, "w");

    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    root = cJSON_CreateObject();
    if (root == NULL) {
        fprintf(stderr, "Error: Cannot allocate JSON root for %s\n", path);
        fclose(pf);
        return 0;
    }

    states = cJSON_AddArrayToObject(root, "states");
    if (states == NULL) {
        fprintf(stderr, "Error: Cannot allocate states array for %s\n", path);
        cJSON_Delete(root);
        fclose(pf);
        return 0;
    }

    for (int n = 0; n < len; n++) {
        cJSON *state_obj = cJSON_CreateObject();
        cJSON *eigenvector = cJSON_CreateArray();

        if (state_obj == NULL || eigenvector == NULL) {
            cJSON_Delete(eigenvector);
            cJSON_Delete(state_obj);
            fprintf(stderr, "Error: Cannot allocate state JSON for %s\n", path);
            cJSON_Delete(root);
            fclose(pf);
            return 0;
        }

        cJSON_AddNumberToObject(state_obj, "index", n + 1);
        cJSON_AddNumberToObject(state_obj, "mass", mass->value[n]);
        if (rms12 && n < rms12->len) {
            cJSON_AddNumberToObject(state_obj, "rms_r12", rms12->value[n]);
        }
        if (rms13 && n < rms13->len) {
            cJSON_AddNumberToObject(state_obj, "rms_r13", rms13->value[n]);
        }
        if (rms23 && n < rms23->len) {
            cJSON_AddNumberToObject(state_obj, "rms_r23", rms23->value[n]);
        }

        for (int d = 0; d < vector->col; d++) {
            cJSON_AddItemToArray(eigenvector, cJSON_CreateNumber(vector->value[n][d]));
        }

        cJSON_AddItemToObject(state_obj, "eigenvector", eigenvector);
        cJSON_AddItemToArray(states, state_obj);
    }

    json_text = cJSON_Print(root);
    if (json_text == NULL) {
        fprintf(stderr, "Error: Cannot serialize JSON output for %s\n", path);
        cJSON_Delete(root);
        fclose(pf);
        return 0;
    }

    if (fputs(json_text, pf) == EOF) {
        fprintf(stderr, "Error: Failed writing JSON output to %s\n", path);
        cJSON_Delete(root);
        free(json_text);
        fclose(pf);
        return 0;
    }

    state = fclose(pf) == 0 ? 1 : 0;
    cJSON_Delete(root);
    free(json_text);

    return state;
}

int write_meson_wfn(const argsInput_t *input, const matrix_t *vector)
{
    if (input == NULL || vector == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_wfn()\n");
        return 0;
    }

    double rmin = 0.01;
    double rmax = 10.0;
    double dr = 0.01;
    int nstate = vector->row;
    int state = 1;

    for (int n = 0; n < nstate; n++) {
        FILE *pf;
        char path[275];
        double normalized;

        sprintf(path, "%s.wfn.%d.dat", input->project, n + 1);
        pf = fopen(path, "w");
        if (pf == NULL) {
            fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
            state = 0;
            continue;
        }

        normalized = get_normalized_factor(input, vector->value[n]);
        fprintf(pf, "# r[fm]    R(r)[fm^(-3/2)]    u(r)=rR(r)[fm^(-1/2)]    r^2|R(r)|^2[fm^(-1)]\n");

        for (double r = rmin; r <= rmax; r += dr) {
            double wfn = get_state_wfn_value(input, vector->value[n], normalized, r);
            double uwfn = r * wfn;
            double density = r * r * wfn * wfn;
            fprintf(pf, "%.8f    %.8e    %.8e    %.8e\n", r, wfn, uwfn, density);
        }

        if (fclose(pf) != 0) {
            fprintf(stderr, "Error: Failed to close file %s\n", path);
            state = 0;
        }
    }

    return state;
}

int write_meson_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc)
{
    if (input == NULL || args_model == NULL || args_dynmc == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_pot()\n");
        return 0;
    }

    /* prepare variables */
    int f1 = input->f1, f2 = input->f2;
    double s1 = 0.5, s2 = 0.5;
    double S = input->S, L = input->L, J = input->J;
    double m1 = getmq(f1, args_model);
    double m2 = getmq(f2, args_model);
    double sigmaij = sigma_ij(m1, m2, args_model->sigma_0, args_model->s);
    args_dynmc->mi = m1;
    args_dynmc->mj = m2;
    args_dynmc->Cij = -4.0 / 3.0;
    args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->Sigij = sigmaij;
    sigma_k_ij(sigmaij, args_dynmc->Sigkij);

    FILE *pf;
    char path[264];
    sprintf(path, "%s%s", input->project, ".pot.dat");
    pf = fopen(path, "w");
    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    double fm = 5.06773093854369882649;
    double rmin = 0.01;
    double rmax = 10.0;
    double dr = 0.01;
    gi_pot_ctx_t ctx = { args_model, args_dynmc };
    for (double r = rmin; r <= rmax; r += dr) {
        double rGeV = r * fm; /* convert fm to GeV^-1 */
        double potential = GIVconf(rGeV, &ctx)
            + GIVcoul(rGeV, &ctx)
            + GIVcont(rGeV, &ctx)
            + GIVsovi(rGeV, &ctx)
            + GIVsovj(rGeV, &ctx)
            + GIVsovij(rGeV, &ctx)
            + GIVsosi(rGeV, &ctx)
            + GIVsosj(rGeV, &ctx)
            + GIVtens(rGeV, &ctx);

        fprintf(pf, "%.8f    %.8e\n", r, potential);
    }

    if (fclose(pf) != 0) {
        fprintf(stderr, "Error: Failed to close file %s\n", path);
        return 0;
    }

    return 1;
}

int write_meson_nr_pot(const argsInput_t *input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc)
{
    if (input == NULL || args_model == NULL || args_dynmc == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_nr_pot()\n");
        return 0;
    }

    int f1 = input->f1, f2 = input->f2;
    double s1 = 0.5, s2 = 0.5;
    double S = input->S, L = input->L, J = input->J;
    args_dynmc->mi = getmq(f1, args_model);
    args_dynmc->mj = getmq(f2, args_model);
    args_dynmc->Cij = -4.0 / 3.0;
    args_dynmc->OCent = operator_center_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OSdS = operator_sdots_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSi = operator_ldotsi_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OLSj = operator_ldotsj_sl(s1, s2, S, L, s1, s2, S, L, J);
    args_dynmc->OTens = operator_tensor_sl(s1, s2, S, L, s1, s2, S, L, J);

    FILE *pf;
    char path[264];
    sprintf(path, "%s%s", input->project, ".pot.dat");
    pf = fopen(path, "w");
    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    double fm = 5.06773093854369882649;
    double rmin = 0.01;
    double rmax = 10.0;
    double dr = 0.01;
    nr_pot_ctx_t ctx = { args_model, args_dynmc };
    for (double r = rmin; r <= rmax; r += dr) {
        double rGeV = r * fm;
        double potential = NRVconf(rGeV, &ctx)
            + NRVcoul(rGeV, &ctx)
            + NRVcont(rGeV, &ctx)
            + NRVsov(rGeV, &ctx)
            + NRVsos(rGeV, &ctx)
            + NRVtens(rGeV, &ctx);

        fprintf(pf, "%.8f    %.8e\n", r, potential);
    }

    if (fclose(pf) != 0) {
        fprintf(stderr, "Error: Failed to close file %s\n", path);
        return 0;
    }

    return 1;
}

int write_baryon_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc)
{
    if (input == NULL || args_model == NULL || args_dynmc == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_baryon_pot()\n");
        return 0;
    }

    FILE *pf;
    char path[264];
    sprintf(path, "%s%s", input->project, ".pot.dat");
    pf = fopen(path, "w");
    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    double fm = 5.06773093854369882649;
    double rmin = 0.01;
    double rmax = 10.0;
    double dr = 0.01;
    gi_pot_ctx_t ctx = { args_model, args_dynmc };
    args_dynmc->Cij = -2.0 / 3.0;
    args_dynmc->OCent = 1.0;
    args_dynmc->OSdS = 0.0;
    args_dynmc->OLSi = 0.0;
    args_dynmc->OLSj = 0.0;
    args_dynmc->OTens = 0.0;
    args_dynmc->system = SYSTEM_BARYON;

    fprintf(pf, "# r[fm]    V12[GeV]    V13[GeV]    V23[GeV]    Vsum[GeV]  (central GI, colour -2/3)\n");
    for (double r = rmin; r <= rmax; r += dr) {
        double rGeV = r * fm;
        double Vpair[3];
        int flav[3][2] = {{input->f1, input->f2}, {input->f3, input->f1}, {input->f2, input->f3}};
        for (int p = 0; p < 3; p++) {
            double mi = getmq(flav[p][0], args_model);
            double mj = getmq(flav[p][1], args_model);
            args_dynmc->mi = mi;
            args_dynmc->mj = mj;
            args_dynmc->Sigij = sigma_ij(mi, mj, args_model->sigma_0, args_model->s);
            sigma_k_ij(args_dynmc->Sigij, args_dynmc->Sigkij);
            Vpair[p] = GIVconf(rGeV, &ctx) + GIVcoul(rGeV, &ctx);
        }
        fprintf(pf, "%.8f    %.8e    %.8e    %.8e    %.8e\n",
            r, Vpair[0], Vpair[1], Vpair[2], Vpair[0] + Vpair[1] + Vpair[2]);
    }

    if (fclose(pf) != 0) {
        fprintf(stderr, "Error: Failed to close file %s\n", path);
        return 0;
    }

    return 1;
}