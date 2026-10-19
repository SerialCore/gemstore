/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/cmisym.h>
#include <gemstore/basis/color.h>
#include <gemstore/basis/spin.h>

#include <math.h>
#include <stdlib.h>

/* nq; pair 12/45: 1 = 6/6̄, 0 = 3̄/3;
   cluster 123/456: 2 fully symmetric, 1 mixed, 0 fully antisymmetric; -1 unused. */
static const struct {
    int nq;
    int c12, c123, c45, c456;
    intrin_wfn_t (*wfn)(void);
} color_ch[] = {
    [CMI_COLOR_MESON]    = {2, -1, -1, -1, -1, color_wfn_meson},
    [CMI_COLOR_BARYON]   = {3, -1, -1, -1, -1, color_wfn_baryon},
    [CMI_COLOR_TETRA11]  = {4, -1, -1, -1, -1, color_wfn_tetra11},
    [CMI_COLOR_TETRA33]  = {4, -1, -1, -1, -1, color_wfn_tetra33},
    [CMI_COLOR_TETRA66]  = {4, -1, -1, -1, -1, color_wfn_tetra66},
    [CMI_COLOR_TETRA88]  = {4, -1, -1, -1, -1, color_wfn_tetra88},
    [CMI_COLOR_PENTA688] = {5,  1,  1, -1, -1, color_wfn_penta688},
    [CMI_COLOR_PENTA388] = {5,  0,  1, -1, -1, color_wfn_penta388},
    [CMI_COLOR_PENTA11]  = {5,  0,  0, -1, -1, color_wfn_penta11},
    [CMI_COLOR_HEXA1010] = {6,  1,  2,  1,  2, color_wfn_hexa1010},
    [CMI_COLOR_HEXA6868] = {6,  1,  1,  1,  1, color_wfn_hexa6868},
    [CMI_COLOR_HEXA6838] = {6,  1,  1,  0,  1, color_wfn_hexa6838},
    [CMI_COLOR_HEXA3868] = {6,  0,  1,  1,  1, color_wfn_hexa3868},
    [CMI_COLOR_HEXA3838] = {6,  0,  1,  0,  1, color_wfn_hexa3838},
    [CMI_COLOR_HEXA11]   = {6,  0,  0,  0,  0, color_wfn_hexa11}
};

static int color_ok(int id)
{
    return id >= 0 && id <= CMI_COLOR_HEXA11 && color_ch[id].wfn != NULL;
}

static int eq_spin(double a, double b)
{
    return fabs(a - b) < 1e-12;
}

/* Identical flavors: color 6 with s=0, color 3̄ with s=1. */
static int pauli_pair(int fa, int fb, int csym, double s)
{
    int ssym;

    if (fa != fb) return 1;
    if (eq_spin(s, 1.0)) ssym = 1;
    else if (eq_spin(s, 0.0)) ssym = 0;
    else return 0;
    if (csym < 0) return 0;
    return csym != ssym;
}

/* Identical flavors: color 10 forbidden; color 1 with s=3/2; color 8 with s=1/2. */
static int pauli_three(int fa, int fb, int fc, int csym, double s12, double s123)
{
    int ssym;

    if (!(fa == fb && fb == fc)) return 1;
    if (csym < 0) return 0;
    if (eq_spin(s123, 1.5)) ssym = eq_spin(s12, 1.0) ? 1 : -1;
    else if (eq_spin(s123, 0.5)) ssym = 0;
    else ssym = -1;
    if (ssym < 0) return 0;
    if (csym == 2) return 0;
    if (csym == 0) return ssym == 1;
    return ssym == 0;
}

int cmi_color_wfn(int id, intrin_wfn_t *wfn)
{
    if (wfn == NULL || !color_ok(id)) return -1;
    *wfn = color_ch[id].wfn();
    return 0;
}

int cmi_spin_wfn(int nquark, const cmi_spin_t *channel, intrin_wfn_t *wfn)
{
    if (channel == NULL || wfn == NULL) return -1;
    switch (nquark) {
        case 2: *wfn = spin_wfn_meson(channel->st, channel->st3); break;
        case 3: *wfn = spin_wfn_baryon(channel->s12, channel->st, channel->st3); break;
        case 4: *wfn = spin_wfn_tetra(channel->s12, channel->s34, channel->st, channel->st3); break;
        case 5: *wfn = spin_wfn_penta(channel->s12, channel->s123, channel->s45, channel->st, channel->st3); break;
        case 6: *wfn = spin_wfn_hexa(channel->s12, channel->s123, channel->s45, channel->s456, channel->st, channel->st3); break;
        default: return -1;
    }
    return 0;
}

int cmi_color_spin_allowed(int color_id, const cmi_spin_t *channel, int nf, const int flavor[])
{
    if (channel == NULL || flavor == NULL || !color_ok(color_id)) return 0;
    if (color_ch[color_id].nq != nf) return 0;
    if (nf < 5) return 1;

    if (!pauli_pair(flavor[0], flavor[1], color_ch[color_id].c12, channel->s12)) return 0;
    if (!pauli_three(flavor[0], flavor[1], flavor[2], color_ch[color_id].c123, channel->s12, channel->s123)) return 0;
    if (nf == 6) {
        if (!pauli_pair(flavor[3], flavor[4], color_ch[color_id].c45, channel->s45)) return 0;
        if (!pauli_three(flavor[3], flavor[4], flavor[5], color_ch[color_id].c456, channel->s45, channel->s456)) return 0;
    }
    return 1;
}

int cmi_load_basis(const argsInput_t *input, intrin_wfn_t *color, intrin_wfn_t *spin)
{
    int nquark;

    if (input == NULL || color == NULL || spin == NULL) return -1;
    if (input->basis != BASIS_CMI || input->ncolor < 1 || input->nspin < 1) return -1;
    nquark = input->nn + input->ns + input->nc + input->nb;
    if (nquark < 2 || nquark > 6) return -1;

    for (int i = 0; i < input->ncolor; i++) {
        int id = input->color_id[i];

        if (!color_ok(id) || color_ch[id].nq != nquark) return -1;
        if (cmi_color_wfn(id, &color[i]) != 0) return -1;
    }
    for (int i = 0; i < input->nspin; i++) {
        if (cmi_spin_wfn(nquark, &input->spin_ch[i], &spin[i]) != 0) return -1;
    }
    return 0;
}
