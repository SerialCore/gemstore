/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/print.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
    printf("Usage: gemstore [--input FILE] [--fitting TARGET] [--print ITEM] [--debug UNIT]\n\n");
    printf("Arguments:\n");
    printf("  -i, --input           input FILE that constains full instructions\n");
    printf("  -f, --fitting         fit TARGET such as GIScreen_meson, GIScreen_ccbar, GIScreen_bbbar, GIQuadra_light\n");
    printf("  -d, --debug           debug UNIT such as su3_product, soc_operator, casimir_operator, \n");
    printf("                        color_wfn, spin_wfn, isospin_wfn, orbit_wfn, eigen_system\n");
    printf("  -p, --print           print ITEM such as potential, wavefunction\n");
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
    printf("  Task Type:            %-50d\n", input->task);
    printf("  Model Type:           %-50d\n", input->model);
    printf("  System Type:          %-50d\n", input->system);
    printf("  Basis Type:           %-50d\n", input->orbit);
    printf("\n");

    /* Quark Flavor Configuration */
    printf("Quark Flavor Configuration:\n");
    printf("  Flavor 1 (f1):        %-50d\n", input->f1);
    printf("  Flavor 2 (f2):        %-50d\n", input->f2);
    printf("  Flavor 3 (f3):        %-50d\n", input->f3);
    printf("  Flavor 4 (f4):        %-50d\n", input->f4);
    printf("\n");

    /* Angular Momentum Quantum Numbers */
    printf("Angular Momentum Quantum Numbers:\n");
    printf("  Spin Momentum (S):    %-50.6f\n", input->S);
    printf("  Orbital Momentum (L): %-50.6f\n", input->L);
    printf("  jl (Orbital jl):      %-50.6f\n", input->jl);
    printf("  Total Momentum (J):   %-50.6f\n", input->J);
    printf("\n");

    /* Basis Parameters */
    printf("Basis Parameters:\n");
    if (input->orbit == ORBIT_GEM || input->orbit == ORBIT_CRG || input->orbit == ORBIT_CSM) {
        printf("  Number of Gaussians (nmax): %-46d\n", input->nmax);
        printf("  Minimum Range (rmin):       %-46.6f fm\n", input->rmin);
        printf("  Maximum Range (rmax):       %-46.6f fm\n", input->rmax);
    }
    if (input->orbit == ORBIT_CRG) {
        printf("  Oscillation Scale (omega):  %-46.6f\n", input->omega);
    }
    if (input->orbit == ORBIT_CSM) {
        printf("  Rotation Angle (theta):     %-46.6f\n", input->theta);
    }
    if (input->orbit == ORBIT_SHO) {
        printf("  Harmonic Scale (beta):      %-46.6f\n", input->beta);
    }
    printf("\n");

    /* Model Parameters */
    printf("Model Parameters:\n");
    printf("  Quark Masses:\n");
    printf("    Light Quark (n)   (mn):  %-46.6f GeV\n", input->params.mn);
    printf("    Strange Quark (s) (ms):  %-46.6f GeV\n", input->params.ms);
    printf("    Charm Quark (c)   (mc):  %-46.6f GeV\n", input->params.mc);
    printf("    Bottom Quark (b)  (mb):  %-46.6f GeV\n", input->params.mb);
    printf("\n");

    printf("  Potential Parameters:\n");
    printf("    String Tension (b):      %-46.6f\n", input->params.b);
    printf("    Screen Length (mu):      %-46.6f\n", input->params.mu);
    printf("    Constant Potential (c):  %-46.6f\n", input->params.c);
    printf("\n");

    printf("  Gaunov-Isgur Smearing Parameters:\n");
    printf("    sigma_0 (Center):        %-46.6f\n", input->params.sigma_0);
    printf("    sigma (Spin-Spin):       %-46.6f\n", input->params.s);
    printf("    epsilon_cont (Contact):  %-46.6f\n", input->params.epsilon_cont);
    printf("    epsilon_sov (Spin-Orbit):%-46.6f\n", input->params.epsilon_sov);
    printf("    epsilon_sos (Thomas):    %-46.6f\n", input->params.epsilon_sos);
    printf("    epsilon_tens (Tensor):   %-46.6f\n", input->params.epsilon_tens);
    printf("\n");
}

void print_debug_results(const array_t *eigenvalue, const array_t *rmsradius, const matrix_t *eigenvector, int nmax)
{
    if (eigenvalue == NULL || rmsradius == NULL || eigenvector == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                            DEBUG RESULTS SUMMARY                             \n");
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

    char path[265];
    sprintf(path, "%s%s", input->project, ".out.json");
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

    /* Timestamp and metadata */
    time_t now = time(NULL);
    struct tm *timeinfo = localtime(&now);
    char time_str[80];
    strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", timeinfo);
    cJSON_AddStringToObject(root, "generated", time_str);
    cJSON_AddStringToObject(root, "project", input->project);
    cJSON_AddNumberToObject(root, "task", input->task);
    cJSON_AddNumberToObject(root, "model", input->model);
    cJSON_AddNumberToObject(root, "system", input->system);
    cJSON_AddNumberToObject(root, "basis", input->orbit);

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
