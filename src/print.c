/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/print.h>
#include <gemstore/parse.h>
#include <gemstore/input.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>
#include <gemstore/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    printf("Usage: gemstore [--compute FILE] [--debug UNIT]\n\n");
    printf("Arguments:\n");
    printf("  -c, --compute         run FILE. task is SPECTRA or FITTING\n");
    printf("  -d, --debug           debug UNIT such as su3_product, soc_operator, casimir_operator, \n");
    printf("                        color_wfn, spin_wfn, isospin_wfn, orbit_wfn, eigen_system,\n");
    printf("                        eigen_system_complex\n");
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
    if (input->task == TASK_FITTING) {
        printf("  Fit Target:           %-50s\n", fitting_type_str[input->target]);
    }
    else {
        printf("  Parameter Set:        %-50s\n", param_type_str[input->param]);
    }
    printf("  System Type:          %-50s\n", system_type_str[input->system]);
    printf("  Basis Type:           %-50s\n", basis_type_str[input->basis]);
    if (input->param == PARAM_GISTRING_CUSTOM || input->param == PARAM_GISCREEN_CUSTOM
        || input->param == PARAM_NRSTRING_CUSTOM || input->param == PARAM_NRSCREEN_CUSTOM
        || input->param == PARAM_BAGMIT_CUSTOM || input->param == PARAM_OBESTAT_CUSTOM) {
        printf("  Parameter File:       %-50s\n", input->param_file);
    }
    printf("\n");

    /* Basis Parameters */
    printf("Basis Parameters:\n");
    if (input->basis == BASIS_GEM) {
        printf("  Number of Gaussians (nmax): %-46d\n", input->nmax);
        printf("  Minimum Range (rmin):       %-46.6f fm\n", input->rmin);
        printf("  Maximum Range (rmax):       %-46.6f fm\n", input->rmax);
    }
    if (input->basis == BASIS_SHO) {
        printf("  Number of Oscillators (nmax): %-44d\n", input->nmax);
        printf("  Harmonic Scale (beta):        %-44.6f GeV\n", input->beta);
    }
    if (input->basis == BASIS_NONE) {
        printf("  No color-magnetic basis\n");
    }
    if (input->basis == BASIS_CMI) {
        printf("  Color channels:");
        for (int i = 0; i < input->ncolor; i++) printf(" %s", cmicolor_type_str[input->color_id[i]]);
        printf("\n");
        printf("  Spin channels: %d\n", input->nspin);
        for (int i = 0; i < input->nspin; i++) {
            const cmi_spin_t *ch = &input->spin_ch[i];
            printf("    %d  s12=%.1f s34=%.1f s123=%.1f s45=%.1f s456=%.1f st=%.1f st3=%.1f\n",
                i + 1, ch->s12, ch->s34, ch->s123, ch->s45, ch->s456, ch->st, ch->st3);
        }
    }
    printf("\n");

    if (input->task == TASK_FITTING) {
        printf("Fitting task running ...\n");
        return;
    }

    if (input->system == SYSTEM_MOLECUL2 || input->system == SYSTEM_MOLECUL3) {
        printf("Hadron Configuration:\n");
    }
    else {
        printf("Quark Flavor Configuration:\n");
    }
    if (input->system == SYSTEM_MESON) {
        printf("  Flavor 1 (f1):        %-50d\n", input->f1);
        printf("  Flavor 2 (f2):        %-50d\n", input->f2);
    }
    else if (input->system == SYSTEM_BARYON) {
        printf("  Flavor 1 (f1):        %-50d\n", input->f1);
        printf("  Flavor 2 (f2):        %-50d\n", input->f2);
        printf("  Flavor 3 (f3):        %-50d\n", input->f3);
    }
    else if (input->system == SYSTEM_BAG) {
        printf("  n quarks (nn):        %-50d\n", input->nn);
        printf("  s quarks (ns):        %-50d\n", input->ns);
        printf("  c quarks (nc):        %-50d\n", input->nc);
        printf("  b quarks (nb):        %-50d\n", input->nb);
    }
    else if (input->system == SYSTEM_MOLECUL2) {
        printf("  Leg 1 type:           %-50s\n", obe_leg_type_str[input->type1]);
        printf("  Leg 2 type:           %-50s\n", obe_leg_type_str[input->type2]);
        printf("  Physical mass m1:     %-50.6f GeV\n", input->m1);
        printf("  Physical mass m2:     %-50.6f GeV\n", input->m2);
        printf("  cc1, cc2:             %d, %d\n", input->cc1, input->cc2);
        printf("  eta:                  %-50.6f\n", input->eta);
        printf("  iso_pi (pion):        %-50.6f\n", input->iso_pi);
        printf("  iso_rho (rho):        %-50.6f\n", input->iso_rho);
    }
    printf("\n");

    /* Angular Momentum Quantum Numbers */
    printf("Angular Momentum Quantum Numbers:\n");
    if (input->system == SYSTEM_MESON || input->system == SYSTEM_MOLECUL2) {
        printf("  Spin Momentum (S):    %-50.6f\n", input->S);
        printf("  Orbital Momentum (L): %-50.6f\n", input->L);
        printf("  Total Momentum (J):   %-50.6f\n", input->J);
    }
    else if (input->system == SYSTEM_BARYON) {
        printf("  1<->2 symmetry:       %-50d\n", input->sym12);
        printf("  Parity (P):           %-50d\n", input->P);
        printf("  Jacobi Lmax:          %-50d\n", input->Lmax);
        printf("  Total Momentum (J):   %-50.6f\n", input->J);
    }
    else if (input->system == SYSTEM_BAG) {
        printf("  Defined in color-magnetic basis parameters\n");
    }
    printf("\n");

    /* Model Parameters */
    printf("Model Parameters:\n");
    if (input->model == MODEL_GISTRING || input->model == MODEL_GISCREEN) {
        argsGIModel_t args_model = argsGIModel_from(input);
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
    else if (input->model == MODEL_NRSTRING || input->model == MODEL_NRSCREEN) {
        argsNRModel_t args_model = argsNRModel_from(input);
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
        printf("    Coupling Constant (as):  %-46.6f\n", args_model.alpha_s);
        printf("    Smearing (sigma):        %-46.6f\n", args_model.sigma);
        printf("\n");
    }
    else if (input->model == MODEL_BAGMIT) {
        argsBAGModel_t bag = argsBAGModel_from(input);
        printf("  Quark Masses:\n");
        printf("    Light Quark (n)   (mn):  %-46.6f GeV\n", bag.mn);
        printf("    Strange Quark (s) (ms):  %-46.6f GeV\n", bag.ms);
        printf("    Charm Quark (c)   (mc):  %-46.6f GeV\n", bag.mc);
        printf("    Bottom Quark (b)  (mb):  %-46.6f GeV\n", bag.mb);
        printf("\n");
    
        printf("  Bag Parameters:\n");
        printf("    Bag Constant (B):        %-44.6e GeV^4\n", bag.B);
        printf("    Zero Point (Z):          %-46.6f\n", bag.Z);
        printf("\n");
    }
    else if (input->model == MODEL_OBESTAT) {
        argsOBEModel_t obe = argsOBEModel_from(input);
        printf("  Parameter file (cutoff, tensor masses, couplings):\n");
        printf("    Lambda:                  %-46.6f GeV\n", obe.lambda);
        printf("    Tensor mass M1:          %-46.6f GeV\n", obe.tensor_mass1);
        printf("    Tensor mass M2:          %-46.6f GeV\n", obe.tensor_mass2);
        printf("    g_pi:                    %.6f, %.6f\n", obe.g_pi1, obe.g_pi2);
        printf("    g_sigma:                 %.6f, %.6f\n", obe.g_sigma1, obe.g_sigma2);
        printf("    g_rho:                   %.6f, %.6f\n", obe.g_rho1, obe.g_rho2);
        printf("    f_rho:                   %.6f, %.6f\n", obe.f_rho1, obe.f_rho2);
        printf("    g_omega:                 %.6f, %.6f\n", obe.g_omega1, obe.g_omega2);
        printf("    f_omega:                 %.6f, %.6f\n", obe.f_omega1, obe.f_omega2);
        printf("\n");
        printf("  Legs, iso, eta, and cc are taken from the system JSON.\n");
        printf("\n");
    }
}

void print_meson_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap, 
    const array_t *rms12, int len)
{
    if (eigenvalue == NULL || eigenvector == NULL || overlap == NULL || rms12 == NULL) {
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
    double min_rms = rms12->value[0];
    double max_rms = rms12->value[0];
    double total_norm = 0.0;
    int num_states = 0;

    for (int n = 0; n < len; n++) {
        if (n < eigenvalue->len) {
            if (eigenvalue->value[n] < min_mass) min_mass = eigenvalue->value[n];
            if (eigenvalue->value[n] > max_mass) max_mass = eigenvalue->value[n];
        }
        if (n < rms12->len) {
            if (rms12->value[n] < min_rms) min_rms = rms12->value[n];
            if (rms12->value[n] > max_rms) max_rms = rms12->value[n];
        }
        num_states++;
    }

    /* Print header row */
    printf("%-6s%-12s%-12s%-12s%-15s\n",
        "State", "Mass(GeV)", "RMS(fm)", "max|c|", "c^T N c");
    printf("------+-----------+-----------+-----------+--------------\n");

    /* Print each state's information */
    for (int n = 0; n < len; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        
        int nbas = eigenvector->col;
        for (int i = 0; i < nbas; i++) {
            double c = fabs(eigenvector->value[n][i]);
            if (c > maxc) maxc = c;
        }
        for (int i = 0; i < nbas; i++) {
            for (int j = 0; j < nbas; j++) {
                norm += eigenvector->value[n][i] * overlap->value[i][j] * eigenvector->value[n][j];
            }
        }
        
        total_norm += norm;
        
        printf("%-6d%-12.6f%-12.6f%-12.6f%-15.10f\n",
            n+1, eigenvalue->value[n], rms12->value[n], maxc, norm);
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", num_states);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n", 
        min_mass, max_mass, max_mass - min_mass);
    printf("  RMS radius range:    %.6f - %.6f fm (Δ=%.6f fm)\n", 
        min_rms, max_rms, max_rms - min_rms);
    printf("  Total c^T N c:       %.10f\n", total_norm);
    printf("\n");
}

void print_mol2_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap,
    const array_t *rms12, double threshold, int len)
{
    if (eigenvalue == NULL || eigenvector == NULL || overlap == NULL || rms12 == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                          MOLECULE RESULTS SUMMARY                            \n");
    printf("================================================================================\n");
    printf("\n");
    printf("  Threshold m1+m2:     %.6f GeV\n", threshold);
    printf("  Ebind = Mass - (m1+m2). Negative is bound.\n");
    printf("\n");

    double min_mass = eigenvalue->value[0];
    double max_mass = eigenvalue->value[0];
    double min_rms = rms12->value[0];
    double max_rms = rms12->value[0];
    double total_norm = 0.0;
    int num_states = 0;

    for (int n = 0; n < len; n++) {
        if (n < eigenvalue->len) {
            if (eigenvalue->value[n] < min_mass) min_mass = eigenvalue->value[n];
            if (eigenvalue->value[n] > max_mass) max_mass = eigenvalue->value[n];
        }
        if (n < rms12->len) {
            if (rms12->value[n] < min_rms) min_rms = rms12->value[n];
            if (rms12->value[n] > max_rms) max_rms = rms12->value[n];
        }
        num_states++;
    }

    printf("%-6s%-12s%-12s%-12s%-12s%-15s\n",
        "State", "Mass(GeV)", "Ebind(GeV)", "RMS(fm)", "max|c|", "c^T N c");
    printf("------+-----------+----------+-----------+-----------+--------------\n");

    for (int n = 0; n < len; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        int nbas = eigenvector->col;
        for (int i = 0; i < nbas; i++) {
            double c = fabs(eigenvector->value[n][i]);
            if (c > maxc) maxc = c;
        }
        for (int i = 0; i < nbas; i++) {
            for (int j = 0; j < nbas; j++) {
                norm += eigenvector->value[n][i] * overlap->value[i][j] * eigenvector->value[n][j];
            }
        }
        total_norm += norm;
        printf("%-6d%-12.6f%-12.6f%-12.6f%-12.6f%-15.10f\n",
            n + 1, eigenvalue->value[n], eigenvalue->value[n] - threshold,
            rms12->value[n], maxc, norm);
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", num_states);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n",
        min_mass, max_mass, max_mass - min_mass);
    printf("  Binding range:       %.6f - %.6f GeV\n",
        min_mass - threshold, max_mass - threshold);
    printf("  RMS radius range:    %.6f - %.6f fm (Δ=%.6f fm)\n",
        min_rms, max_rms, max_rms - min_rms);
    printf("  Total c^T N c:       %.10f\n", total_norm);
    printf("\n");
}

void print_baryon_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, int len)
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

    printf("%-6s%-12s%-10s%-10s%-10s%-12s%-15s\n",
        "State", "Mass(GeV)", "r12(fm)", "r13(fm)", "r23(fm)", "max|c|", "c^T N c");
    printf("------+-----------+---------+---------+---------+-----------+--------------\n");

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

        double r12 = (rms12 && n < rms12->len) ? rms12->value[n] : 0.0;
        double r13 = (rms13 && n < rms13->len) ? rms13->value[n] : 0.0;
        double r23 = (rms23 && n < rms23->len) ? rms23->value[n] : 0.0;

        printf(" %-5d%-12.6f%-10.4f%-10.4f%-10.4f%-12.6f%-15.10f\n",
            n + 1, eigenvalue->value[n], r12, r13, r23, maxc, norm);
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", len);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n",
        min_mass, max_mass, max_mass - min_mass);
    printf("  Total norm sum:      %.10f\n", total_norm);
    printf("\n");
}

void print_bag_spectra(const bag_mit_t *state, int nstate, const char *order)
{
    double min_mass;
    double max_mass;

    if (state == NULL || nstate < 1) return;

    printf("\n");
    printf("================================================================================\n");
    printf("                             BAG RESULTS SUMMARY                              \n");
    printf("================================================================================\n");
    printf("\n");
    printf("  Quark content:       n=%d  s=%d  c=%d  b=%d\n",
        state[0].n[0], state[0].n[1], state[0].n[2], state[0].n[3]);
    if (order == NULL) printf("  Configuration:       no need\n");
    else printf("  Configuration:       %s\n", order);
    printf("\n");

    printf("%-6s%-8s%-12s%-11s%-11s%-11s%-11s%-11s\n",
        "State", "st", "Mass(GeV)", "R(1/GeV)", "x_n", "x_s", "x_c", "x_b");
    printf("------+-------+-----------+----------+----------+----------+----------+---------\n");
    min_mass = state[0].mass;
    max_mass = state[0].mass;
    for (int n = 0; n < nstate; n++) {
        if (state[n].mass < min_mass) min_mass = state[n].mass;
        if (state[n].mass > max_mass) max_mass = state[n].mass;
        if (isfinite(state[n].st)) {
            printf("%-6d%-8.1f%-12.6f%-11.6f%-11.6f%-11.6f%-11.6f%-11.6f\n",
                n, state[n].st, state[n].mass, state[n].R,
                state[n].x[0], state[n].x[1], state[n].x[2], state[n].x[3]);
        }
        else {
            printf("%-6d%-8.1s%-12.6f%-11.6f%-11.6f%-11.6f%-11.6f%-11.6f\n",
                n, "---", state[n].mass, state[n].R,
                state[n].x[0], state[n].x[1], state[n].x[2], state[n].x[3]);
        }
    }
    printf("\n");
    printf("  Number of states:    %d\n", nstate);
    printf("  Mass range:          %.6f - %.6f GeV\n", min_mass, max_mass);
    printf("\n");
}