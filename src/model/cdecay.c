/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Open-flavor two-body decay. The QPC amplitude is qpc_amplitude.
 *   Gamma = 2 pi (EB EC / MA) p F sum |M_LS|^2.
 */

#include <gemstore/model/cdecay.h>
#include <gemstore/model/qpcmodel.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void decay_qpc_meson(const argsInput_t *input, qpc_result_t *result)
{
    const qpc_meson_t *a = &input->qpc_initial;
    const qpc_meson_t *b = &input->qpc_final[0];
    const qpc_meson_t *c = &input->qpc_final[1];
    int la = (int)lround(a->L);
    int ja = (int)lround(a->J);
    int lb = (int)lround(b->L);
    int jb = (int)lround(b->J);
    int lc = (int)lround(c->L);
    int jc = (int)lround(c->J);
    int pa = ((la + 1) & 1) ? -1 : 1;
    int pb = ((lb + 1) & 1) ? -1 : 1;
    int pc = ((lc + 1) & 1) ? -1 : 1;
    int Ls[QPC_WAVE_MAX];
    int Ss[QPC_WAVE_MAX];
    int nwave = 0;
    int y;
    int i;
    double disc;
    double p;
    double m2;
    argsQPCModel_t qpc;

    result->open = 0;
    result->nwave = 0;
    result->p = 0.0;
    result->EB = 0.0;
    result->EC = 0.0;
    result->width = 0.0;
    qpc = argsQPCModel_from(input);
    disc = (a->mass * a->mass - (b->mass + c->mass) * (b->mass + c->mass))
        * (a->mass * a->mass - (b->mass - c->mass) * (b->mass - c->mass));
    if (disc <= 0.0 || a->mass <= 0.0) return;
    p = sqrt(disc) / (2.0 * a->mass);
    for (y = abs(jb - jc); y <= jb + jc; y++) {
        int x;

        for (x = 0; x <= ja + y + 2; x++) {
            int tri = abs(x - y) <= ja && ja <= x + y;
            int parity = ((x & 1) ? -1 : 1) * pb * pc == pa;
            if (!tri || !parity) continue;
            if (nwave >= QPC_WAVE_MAX) {
                fprintf(stderr, "3P0DECAY: more than %d partial waves\n", QPC_WAVE_MAX);
                exit(1);
            }
            Ls[nwave] = x;
            Ss[nwave] = y;
            nwave++;
        }
    }
    m2 = 0.0;
    for (i = 0; i < nwave; i++) {
        double amp = qpc_amplitude(&qpc, p, a, b, c, Ls[i], Ss[i], input->qpc_fl1, input->qpc_fl2);

        result->wave[i].L = Ls[i];
        result->wave[i].S = Ss[i];
        result->wave[i].M = amp;
        m2 += amp * amp;
    }
    result->open = 1;
    result->nwave = nwave;
    result->p = p;
    result->EB = sqrt(b->mass * b->mass + p * p);
    result->EC = sqrt(c->mass * c->mass + p * p);
    result->width = 2.0 * M_PI * (result->EB * result->EC / a->mass) * p * input->qpc_F * m2 * 1000.0;
}
