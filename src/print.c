/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/print.h>
#include <gemstore/parse.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/param/argset.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

void print_logo()
{
    printf("\n");
    printf("  ██████╗  ███████╗  ███╗   ███╗  ███████╗  ████████╗   ██████╗  ██████╗  ███████╗\n");
    printf(" ██╔════╝  ██╔════╝  ████╗ ████║  ██╔════╝  ╚══██╔══╝  ██╔═══██║ ██╔══██╗ ██╔════╝\n");
    printf(" ██║  ███╗ █████╗    ██╔████╔██║  ███████╗     ██║     ██║   ██║ ██████╔╝ █████╗  \n");
    printf(" ██║   ██║ ██╔══╝    ██║╚██╔╝██║  ╚════██║     ██║     ██║   ██║ ██╔══██╗ ██╔══╝  \n");
    printf(" ╚██████╔╝ ███████╗  ██║ ╚═╝ ██║  ███████║     ██║     ╚██████╔╝ ██║  ██║ ███████╗\n");
    printf("  ╚═════╝  ╚══════╝  ╚═╝     ╚═╝  ╚══════╝     ╚═╝      ╚═════╝  ╚═╝  ╚═╝ ╚══════╝\n");
    printf("\n");
}

void print_help()
{
    printf("gemstore: hadron spectroscopy tools using Gaussian Expanding Method, Godfrey-Isgur models and more.\n\n");
    printf("Usage: gemstore [--compute FILE] [--fitting TARGET] [--debug UNIT]\n\n");
    printf("Arguments:\n");
    printf("  -c, --compute         perform computation with input FILE that contains full instructions\n");
    printf("  -f, --fitting         fit TARGET such as GIScreen_meson, GIScreen_ccbar, GIScreen_bbbar, GIQuadra_light\n");
    printf("  -d, --debug           debug UNIT such as su3_product, soc_operator, casimir_operator, \n");
    printf("                        color_wfn, spin_wfn, isospin_wfn, orbit_wfn, eigen_system\n");
    printf("  -h,--help             show this help\n");
    printf("  -v,--version          show version\n\n");
}

void print_copyright()
{
    printf("\n");
    printf("================================================================================\n");
    printf("                              COPYRIGHT NOTICE                                \n");
    printf("================================================================================\n");
    printf("\n");
    printf("GEMSTORE - Hadron Spectroscopy Simulation Tools\n");
    printf("\n");
    printf("Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>\n");
    printf("\n");
    printf("This program is free software: you can redistribute it and/or modify\n");
    printf("it under the terms of the GNU General Public License as published by\n");
    printf("the Free Software Foundation, either version 3 of the License, or\n");
    printf("(at your option) any later version.\n");
    printf("\n");
    printf("This program is distributed in the hope that it will be useful,\n");
    printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n");
    printf("GNU General Public License for more details.\n");
    printf("\n");
    printf("You should have received a copy of the GNU General Public License\n");
    printf("along with this program. If not, see <https://www.gnu.org/licenses/>.\n");
    printf("\n");
    printf("SPDX-License-Identifier: GPL-3.0-or-later\n");
    printf("\n");
    printf("================================================================================\n");
    printf("\n");
}

void print_input_parameters(const argsInput_t *input)
{
    if (input == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                          INPUT PARAMETERS SUMMARY                            \n");
    printf("================================================================================\n");
    printf("\n");

    /* Task and Model Configuration */
    printf("Configuration Settings:\n");
    printf("  Project Name:         %-50s\n", input->project);
    printf("  Task Type:            %-50s\n", task_type_str[input->task]);
    printf("  Model Type:           %-50s\n", model_type_str[input->model]);
    printf("  Parameter Set:        %-50s\n", param_type_str[input->param]);
    printf("  System Type:          %-50s\n", system_type_str[input->system]);
    printf("  Basis Type:           %-50s\n", orbit_type_str[input->orbit]);
    if (input->param == PARAM_GISTRING_CUSTOM || input->param == PARAM_GISCREEN_CUSTOM) {
        printf("  Parameter File:       %-50s\n", input->param_file);
    }
    printf("\n");

    /* Quark Flavor Configuration */
    printf("Quark Flavor Configuration:\n");
    if (input->system == SYSTEM_MESON) {
        printf("  Flavor 1 (f1):        %-50d\n", input->f1);
        printf("  Flavor 2 (f2):        %-50d\n", input->f2);
    }
    else if (input->system == SYSTEM_BARYON) {
        printf("  Flavor 1 (f1):        %-50d\n", input->f1);
        printf("  Flavor 2 (f2):        %-50d\n", input->f2);
        printf("  Flavor 3 (f3):        %-50d\n", input->f3);
    }
    printf("\n");

    /* Angular Momentum Quantum Numbers */
    printf("Angular Momentum Quantum Numbers:\n");
    if (input->system == SYSTEM_MESON) {
        printf("  Spin Momentum (S):    %-50.6f\n", input->S);
        printf("  Orbital Momentum (L): %-50.6f\n", input->L);
    }
    else if (input->system == SYSTEM_BARYON) {
        printf("  1<->2 symmetry:       %-50d\n", input->f12);
        printf("  Parity (P):           %-50d\n", input->P);
        printf("  Jacobi Lmax:          %-50d\n", input->Lmax);
    }
    printf("  Total Momentum (J):   %-50.6f\n", input->J);
    printf("\n");

    /* Basis Parameters */
    printf("Basis Parameters:\n");
    if (input->orbit == ORBIT_GEM || input->orbit == ORBIT_CRG) {
        printf("  Number of Gaussians (nmax): %-46d\n", input->nmax);
        printf("  Minimum Range (rmin):       %-46.6f fm\n", input->rmin);
        printf("  Maximum Range (rmax):       %-46.6f fm\n", input->rmax);
    }
    if (input->orbit == ORBIT_CRG) {
        printf("  Oscillation Scale (omega):  %-46.6f\n", input->omega);
    }
    if (input->orbit == ORBIT_SHO) {
        printf("  Harmonic Scale (beta):      %-46.6f\n", input->beta);
    }
    printf("\n");

    /* Model Parameters */
    argsGIModel_t args_model = argsGIModel_from(input);
    printf("Model Parameters:\n");
    printf("  Quark Masses:\n");
    printf("    Light Quark (n)   (mn):  %-46.6f GeV\n", args_model.mn);
    printf("    Strange Quark (s) (ms):  %-46.6f GeV\n", args_model.ms);
    printf("    Charm Quark (c)   (mc):  %-46.6f GeV\n", args_model.mc);
    printf("    Bottom Quark (b)  (mb):  %-46.6f GeV\n", args_model.mb);
    printf("\n");

    printf("  Potential Parameters:\n");
    printf("    String Tension (b):      %-46.6f\n", args_model.b);
    printf("    Screen Length (mu):      %-46.6f\n", args_model.mu);
    printf("    Constant Potential (c):  %-46.6f\n", args_model.c);
    printf("\n");

    printf("  Gaunov-Isgur Smearing Parameters:\n");
    printf("    sigma_0 (Center):        %-46.6f\n", args_model.sigma_0);
    printf("    sigma (Spin-Spin):       %-46.6f\n", args_model.s);
    printf("    epsilon_cont (Contact):  %-46.6f\n", args_model.epsilon_cont);
    printf("    epsilon_sov (Spin-Orbit):%-46.6f\n", args_model.epsilon_sov);
    printf("    epsilon_sos (Thomas):    %-46.6f\n", args_model.epsilon_sos);
    printf("    epsilon_tens (Tensor):   %-46.6f\n", args_model.epsilon_tens);
    printf("\n");
}

void print_meson_spectra(const array_t *eigenvalue, const array_t *rmsradius, const matrix_t *eigenvector, int nmax)
{
    if (eigenvalue == NULL || rmsradius == NULL || eigenvector == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                            MESON RESULTS SUMMARY                             \n");
    printf("================================================================================\n");
    printf("\n");

    /* Compute global statistics */
    double min_mass = eigenvalue->value[0];
    double max_mass = eigenvalue->value[0];
    double min_rms = rmsradius->value[0];
    double max_rms = rmsradius->value[0];
    double total_norm = 0.0;
    int num_states = 0;

    for (int n = 0; n < nmax; n++) {
        if (n < eigenvalue->len) {
            if (eigenvalue->value[n] < min_mass) min_mass = eigenvalue->value[n];
            if (eigenvalue->value[n] > max_mass) max_mass = eigenvalue->value[n];
        }
        if (n < rmsradius->len) {
            if (rmsradius->value[n] < min_rms) min_rms = rmsradius->value[n];
            if (rmsradius->value[n] > max_rms) max_rms = rmsradius->value[n];
        }
        num_states++;
    }

    /* Print header row */
    printf("%-6s%-12s%-12s%-12s%-15s%-15s\n", 
        "State", "Mass(GeV)", "RMS(fm)", "max|c|", "||c||^2", "norm_check");
    printf("------+-----------+-----------+-----------+--------------+--------------\n");

    /* Print each state's information */
    for (int n = 0; n < nmax; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        
        for (int i = 0; i < nmax; i++) {
            double c = fabs(eigenvector->value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }
        
        total_norm += norm;
        
        /* Determine normalization status */
        const char *norm_status = "";
        if (fabs(norm - 1.0) < 1e-6) {
            norm_status = "✓ OK";
        } else if (norm > 1.0) {
            norm_status = "⚠ OVER";
        } else {
            norm_status = "⚠ UNDER";
        }
        
        printf("%-6d%-12.6f%-12.6f%-12.6f%-15.10f%-15s\n", 
            n+1, eigenvalue->value[n], rmsradius->value[n], maxc, norm, norm_status);
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", num_states);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n", 
        min_mass, max_mass, max_mass - min_mass);
    printf("  RMS radius range:    %.6f - %.6f fm (Δ=%.6f fm)\n", 
        min_rms, max_rms, max_rms - min_rms);
    printf("  Total norm sum:      %.10f (should ≈ %d)\n", total_norm, nmax);
    printf("  Average norm per state: %.10f\n", total_norm / num_states);
    printf("\n");
}

void print_baryon_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap, int len)
{
    if (eigenvalue == NULL || eigenvector == NULL || overlap == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                           BARYON RESULTS SUMMARY                              \n");
    printf("================================================================================\n");
    printf("\n");

    double min_mass = eigenvalue->value[0];
    double max_mass = eigenvalue->value[0];
    double total_norm = 0.0;

    printf("%-6s%-12s%-12s%-15s%-15s\n",
        "State", "Mass(GeV)", "max|c|", "c^T N c", "norm_check");
    printf("------+-----------+-----------+--------------+--------------\n");

    for (int n = 0; n < len; n++) {
        double norm = 0.0;
        double maxc = 0.0;

        if (eigenvalue->value[n] < min_mass) min_mass = eigenvalue->value[n];
        if (eigenvalue->value[n] > max_mass) max_mass = eigenvalue->value[n];

        for (int i = 0; i < eigenvector->col; i++) {
            double c = fabs(eigenvector->value[n][i]);
            if (c > maxc) maxc = c;
        }

        for (int i = 0; i < eigenvector->col; i++) {
            for (int j = 0; j < eigenvector->col; j++) {
                norm += eigenvector->value[n][i] * overlap->value[i][j] * eigenvector->value[n][j];
            }
        }

        total_norm += norm;

        printf("%-6d%-12.6f%-12.6f%-15.10f%-15s\n",
            n + 1, eigenvalue->value[n], maxc, norm,
            (fabs(norm - 1.0) < 1e-6) ? "OK" : "CHECK");
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", len);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n",
        min_mass, max_mass, max_mass - min_mass);
    printf("  Total norm sum:      %.10f\n", total_norm);
    printf("\n");
}

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

int write_baryon_spectra(const argsInput_t *input, const array_t *mass, const matrix_t *vector, int len)
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

static void write_state_wfn(const argsInput_t *input, const double *eigenvector, FILE *file)
{
    if (input == NULL || eigenvector == NULL || file == NULL) {
        return;
    }

    /* variables for printing */
    int L = (int)input->L;
    double rmin = 0.01;
    double rmax = 10.0;
    double dr = 0.01;
    double fm = 5.06773093854369882649; /* fm to GeV^-1 conversion */

    int nmax = input->nmax;
    double normalized = 1.0;
    double overlap_sum = 0.0;
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, input->rmax, input->rmin);
        basis[i].param = input->omega;
    }

    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);
            double overlap_ij = 0.0;

            if (input->orbit == ORBIT_GEM) {
                overlap_ij = integral_nlr_overlap(GRnlr, factor, &basis[i], &basis[j]);
            }
            else if (input->orbit == ORBIT_CRG) {
                overlap_ij = integral_crg_overlap(CGRnlr, factor, &basis[i], &basis[j]);
            }
            else if (input->orbit == ORBIT_SHO) {
                overlap_ij = (i == j) ? 1.0 : 0.0;
            }

            overlap_sum += eigenvector[i] * overlap_ij * eigenvector[j];
        }
    }

    if (fabs(overlap_sum) > 1e-12) {
        normalized = sqrt(1.0 / overlap_sum);
    }

    /* Evaluate wave function at each radial point */
    for (double r = rmin; r <= rmax; r += dr) {
        double rGeV = r * fm; /* convert fm to GeV^-1 for consistency with potential */
        double psi_r = 0.0;

        /* Sum over basis functions: psi(r) = sum_n c[n] * phi_n(r) */
        for (int n = 0; n < nmax; n++) {
            double N = n + 1;
            double c_n = eigenvector[n];
            double basis_func = 0.0;

            /* Compute basis function phi_n(r) depending on basis type */
            if (input->orbit == ORBIT_GEM) {
                double nu = getnu(N, nmax, input->rmax, input->rmin);
                basis_func = GRnlr(rGeV, N, L, nu) * exp(-nu * rGeV * rGeV);
            }
            else if (input->orbit == ORBIT_CRG) {
                double nu = getnu(N, nmax, input->rmax, input->rmin);
                double omega = input->omega;
                complex basis_func_complex = CGRnlr(rGeV, N, L, nu, omega) * exp(-nu * rGeV * rGeV);
                basis_func = creal(basis_func_complex);
            }
            else if (input->orbit == ORBIT_SHO) {
                double beta = input->beta;
                basis_func = SRnlr(rGeV, n, L, beta) * exp(-0.5 * beta * beta * rGeV * rGeV);
            }

            psi_r += c_n * basis_func;
        }

        /* Output the overlap-normalized radial wave function. */
        fprintf(file, "%.8f    %.8e\n", r, psi_r * normalized);
    }

    free(basis);
}

int write_meson_wfn(const argsInput_t *input, const matrix_t *vector)
{
    if (input == NULL || vector == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_wfn()\n");
        return 0;
    }

    int nmax = input->nmax;
    int state = 1;

    for (int n = 0; n < nmax; n++) {
        FILE *pf;
        char path[275];

        sprintf(path, "%s.wfn.%d.dat", input->project, n + 1);
        pf = fopen(path, "w");
        if (pf == NULL) {
            fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
            state = 0;
            continue;
        }

        write_state_wfn(input, vector->value[n], pf);

        if (fclose(pf) != 0) {
            fprintf(stderr, "Error: Failed to close file %s\n", path);
            state = 0;
        }
    }

    return state;
}

int write_potential_GI(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc)
{
    if (input == NULL || args_model == NULL || args_dynmc == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to write_potential_GI()\n");
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
    for (double r = rmin; r <= rmax; r += dr) {
        double rGeV = r * fm; /* convert fm to GeV^-1 */
        double potential = GIVconf(rGeV, args_model, args_dynmc)
            + GIVcoul(rGeV, args_model, args_dynmc)
            + GIVcont(rGeV, args_model, args_dynmc)
            + GIVsovi(rGeV, args_model, args_dynmc)
            + GIVsovj(rGeV, args_model, args_dynmc)
            + GIVsovij(rGeV, args_model, args_dynmc)
            + GIVsosi(rGeV, args_model, args_dynmc)
            + GIVsosj(rGeV, args_model, args_dynmc)
            + GIVtens(rGeV, args_model, args_dynmc);

        fprintf(pf, "%.8f    %.8e\n", r, potential);
    }

    if (fclose(pf) != 0) {
        fprintf(stderr, "Error: Failed to close file %s\n", path);
        return 0;
    }

    return 1;
}
