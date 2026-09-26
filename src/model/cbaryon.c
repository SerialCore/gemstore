/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Baryon SPECTRA with SCDK matrix elements (recycle/mfi.h, debug.h, res.h).
 * Basis is the three Jacobi frames c=1,2,3 (pair 12 / 31 / 23). SCDK
 * converts those Gaussians onto the pair a potential acts on.
 * Overcomplete N is pruned by its eigenvalues before H is solved.
 */

#include <gemstore/model/cbaryon.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/threebody.h>
#include <gemstore/basis/orbit.h>

#include <gemstore/math/scdkme.h>
#include <gemstore/math/sumckdk.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>

#include <gemstore/param/argset.h>
#include <gemstore/thread.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    basis_list qnlist_spfy;
    basis_list qnlist_full;
    matrix_t **mlsj;
    sumckdk_scdk ****scdk[21];
    scdk_vargs_t varg;

    matrix_t Nfi;
    matrix_t VogeG[3];
    matrix_t Vcont[3];
    matrix_t Vtens[3];
    matrix_t Vsovii[3];
    matrix_t Vsovjj[3];
    matrix_t Vsovji[3];
    matrix_t Vsovij[3];
    matrix_t Vstring[3];
    matrix_t Vsosii[3];
    matrix_t Vsosjj[3];
    matrix_t pogeG[3];
    matrix_t pcont[3];
    matrix_t ptens[3];
    matrix_t psovii[3];
    matrix_t psovjj[3];
    matrix_t psovji[3];
    matrix_t psovij[3];
    matrix_t psosii[3];
    matrix_t psosjj[3];
    matrix_t T[3];
    matrix_t rmsr[3];
    matrix_t rmsl[3];
} baryon_job_t;

static const vtype scdk_vfn[21] = {
    vcent, vcent, vcent,
    vcont, vcont, vcont,
    vtens, vtens, vtens,
    vsoii, vsoii, vsoii,
    vsojj, vsojj, vsojj,
    vsoji, vsoji, vsoji,
    vsoij, vsoij, vsoij
};

static const int scdk_fg[21] = {
    1, 1, 1,
    1, 1, 1,
    2, 2, 2,
    3, 3, 3,
    4, 4, 4,
    5, 5, 5,
    6, 6, 6
};

static const int scdk_pair[21] = {
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3
};

static int baryon_collect_mats(baryon_job_t *job, matrix_t **list)
{
    int n = 0;
    list[n++] = &job->Nfi;
    for (int p = 0; p < 3; p++) {
        list[n++] = &job->VogeG[p];
        list[n++] = &job->Vcont[p];
        list[n++] = &job->Vtens[p];
        list[n++] = &job->Vsovii[p];
        list[n++] = &job->Vsovjj[p];
        list[n++] = &job->Vsovji[p];
        list[n++] = &job->Vsovij[p];
        list[n++] = &job->Vstring[p];
        list[n++] = &job->Vsosii[p];
        list[n++] = &job->Vsosjj[p];
        list[n++] = &job->pogeG[p];
        list[n++] = &job->pcont[p];
        list[n++] = &job->ptens[p];
        list[n++] = &job->psovii[p];
        list[n++] = &job->psovjj[p];
        list[n++] = &job->psovji[p];
        list[n++] = &job->psovij[p];
        list[n++] = &job->psosii[p];
        list[n++] = &job->psosjj[p];
        list[n++] = &job->T[p];
        list[n++] = &job->rmsr[p];
        list[n++] = &job->rmsl[p];
    }
    return n;
}

static void add_reduced(matrix_t *H, const matrix_t *vt, const matrix_t *M, matrix_t *tmp)
{
    matrix_productT(vt, M, tmp);
    matrix_sum(H, tmp, H);
}

static void add_sandwiched(matrix_t *H, const matrix_t *vt, const matrix_t *V, const matrix_t *P,
    matrix_t *tV, matrix_t *tP, matrix_t *tmp)
{
    matrix_productT(vt, V, tV);
    matrix_productT(vt, P, tP);
    matrix_sandwich(tV, tP, tmp);
    matrix_sum(H, tV, H);
}

static void mlsj_push(matrix_t *m, double coe, double ms1, double ms2, double ms3, double mrho, double mlam)
{
    int n = m->row;
    m->value = (double **)realloc(m->value, sizeof(double *) * (size_t)(n + 1));
    m->value[n] = (double *)malloc(sizeof(double) * 6);
    m->value[n][0] = coe;
    m->value[n][1] = ms1;
    m->value[n][2] = ms2;
    m->value[n][3] = ms3;
    m->value[n][4] = mrho;
    m->value[n][5] = mlam;
    m->row = n + 1;
    m->col = 6;
}

static void baryon_mlsj_jl(baryon_job_t *job)
{
    basis_list *qnlist = &job->qnlist_spfy;

    job->mlsj = (matrix_t **)malloc(sizeof(matrix_t *) * (size_t)qnlist->len_list);
    for (int i = 0; i < qnlist->len_list; i++) {
        job->mlsj[i] = (matrix_t *)malloc(sizeof(matrix_t) * (size_t)qnlist->len_part[i]);
        for (int j = 0; j < qnlist->len_part[i]; j++) {
            job->mlsj[i][j].value = (double **)malloc(sizeof(double *) * 0);
            job->mlsj[i][j].row = 0;
            job->mlsj[i][j].col = 6;

            double coe = qnlist->qnum[i][j].coe;
            double s1 = qnlist->qnum[i][j].s1;
            double s2 = qnlist->qnum[i][j].s2;
            double s3 = qnlist->qnum[i][j].s3;
            int c = qnlist->qnum[i][j].c;
            int lrho = qnlist->qnum[i][j].lrho;
            int llam = qnlist->qnum[i][j].llam;
            int L = qnlist->qnum[i][j].L;
            double sij = qnlist->qnum[i][j].sij;
            double jl = qnlist->qnum[i][j].jl;
            double J = qnlist->qnum[i][j].J;
            double MJ = J;
            double si, sj, sk;

            getijk(s1, s2, s3, &si, &sj, &sk, c);
            for (double msi = -si; msi <= si + 1e-12; msi += 1.0) {
                for (double msj = -sj; msj <= sj + 1e-12; msj += 1.0) {
                    for (double msk = -sk; msk <= sk + 1e-12; msk += 1.0) {
                        for (int mrho = -lrho; mrho <= lrho; mrho++) {
                            for (int mlam = -llam; mlam <= llam; mlam++) {
                                double cgf = coe
                                    * clebsch_gordan(si, msi, sj, msj, sij, msi + msj)
                                    * clebsch_gordan(lrho, mrho, llam, mlam, L, mrho + mlam)
                                    * clebsch_gordan(sij, msi + msj, L, mrho + mlam, jl, msi + msj + mrho + mlam)
                                    * clebsch_gordan(jl, msi + msj + mrho + mlam, sk, msk, J, MJ);
                                if (cgf != 0.0) {
                                    double ms1, ms2, ms3;
                                    get123(&ms1, &ms2, &ms3, msi, msj, msk, c);
                                    mlsj_push(&job->mlsj[i][j], cgf, ms1, ms2, ms3, (double)mrho, (double)mlam);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void baryon_basis_build(const argsInput_t *input, const argsGIModel_t *model,
    basis_list *spfy, basis_list *full)
{
    double m1 = getmq(input->f1, model);
    double m2 = getmq(input->f2, model);
    double m3 = getmq(input->f3, model);
    double J = input->J;
    int P = input->P;
    int f12 = input->f12;
    int Lmax = input->Lmax;
    double s1 = 0.5, s2 = 0.5, s3 = 0.5;

    basis_list_init(spfy);
    basis_list_init(full);

    /* Three pair frames. Identical quarks: recycle |c_a⟩ + η|c_b⟩
     * (1↔2 maps c=2 onto c=3). Distinguishable nsc: three independent charts. */
    int nchan[4] = {0};
    int id12 = threebody_pair_identical(input->f1, input->f2, input->f3, 1);
    int id31 = threebody_pair_identical(input->f1, input->f2, input->f3, 2);
    int id23 = threebody_pair_identical(input->f1, input->f2, input->f3, 3);

    for (int lrho = 0; lrho <= Lmax; lrho++) {
        for (int llam = 0; llam <= Lmax - lrho; llam++) {
            if (((lrho + llam) % 2 == 0) ? (P != 1) : (P != -1)) {
                continue;
            }
            for (int L = abs(lrho - llam); L <= lrho + llam; L++) {
                for (double sij = 0.0; sij <= 1.0 + 1e-9; sij += 1.0) {
                    double eta = threebody_exchange_eta(f12, sij, lrho);
                    int pauli_ok = (eta == 1.0);
                    int keep1 = !id12 || pauli_ok;
                    int keep2 = !id31 || pauli_ok;
                    int keep3 = !id23 || pauli_ok;

                    for (double jl = fabs(L - sij); jl <= L + sij + 1e-9; jl += 1.0) {
                        if (fabs(jl - s3) - 1e-9 > J || jl + s3 + 1e-9 < J) {
                            continue;
                        }

                        if (id12) {
                            if (keep1) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 1, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[1]++;
                            }
                            if (keep2) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 2, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[2]++;
                                basis_list_push(spfy, 0, -1, -1, eta, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 3, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[3]++;
                            }
                        }
                        else if (id23) {
                            if (keep3) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 3, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[3]++;
                            }
                            if (keep1) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 1, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[1]++;
                                basis_list_push(spfy, 0, -1, -1, eta, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 2, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[2]++;
                            }
                        }
                        else if (id31) {
                            if (keep2) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 2, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[2]++;
                            }
                            if (keep1) {
                                basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 1, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[1]++;
                                basis_list_push(spfy, 0, -1, -1, eta, m1, m2, m3, s1, s2, s3,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 3, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                                nchan[3]++;
                            }
                        }
                        else {
                            basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                0.0, 0.0, 0.0, 0.0, 0.0, 1, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                            nchan[1]++;
                            basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                0.0, 0.0, 0.0, 0.0, 0.0, 2, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                            nchan[2]++;
                            basis_list_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                0.0, 0.0, 0.0, 0.0, 0.0, 3, lrho, llam, L, sij, jl, J, 0, 0, 0.0, 0.0);
                            nchan[3]++;
                        }
                    }
                }
            }
        }
    }

    basis_list_push_full(spfy, full, input->rmin, input->rmax, input->nmax);
    printf("Baryon Jacobi GEM: c=1(12) / c=2(31) / c=3(23)\n");
    printf("  angular channels: c1=%d c2=%d c3=%d  (n_spfy=%d n_full=%d)\n",
        nchan[1], nchan[2], nchan[3], spfy->len_list, full->len_list);
    fflush(stdout);
}

static void calc_scdk_mt(void *args)
{
    argsThread_t *mtargs = (argsThread_t *)args;
    baryon_job_t *job = (baryon_job_t *)mtargs->args;
    int ith = mtargs->index;
    int nf = ith / 21;
    int op = ith % 21;
    int *len_part = job->qnlist_spfy.len_part;
    int len_list = job->qnlist_spfy.len_list;

    for (int nfp = 0; nfp < len_part[nf]; nfp++) {
        for (int ni = 0; ni < len_list; ni++) {
            for (int nip = 0; nip < len_part[ni]; nip++) {
                sumckdk_scdk_vtype(&job->scdk[op][nf][nfp][ni][nip],
                    scdk_vfn[op], scdk_fg[op], job->qnlist_spfy, job->mlsj,
                    nf, nfp, ni, nip, mtargs->mutex, mtargs->lock, scdk_pair[op]);
            }
        }
    }
}

static double me_add(txrp_vtype kin, sumckdk_scdk scdk, basis_qnum qf, basis_qnum qi,
    scdk_vargs_t varg, scdk_inte_fn iv, int c)
{
    return inteVcenPartA(kin, scdk, qf, qi, varg, iv, c);
}

static void getmfi(void *args)
{
    argsThread_t *mtarg = (argsThread_t *)args;
    baryon_job_t *job = (baryon_job_t *)mtarg->args;
    int nf = mtarg->index;

    for (int ni = 0; ni < job->qnlist_full.len_list; ni++) {
        double nfi = 0.0;
        double VogeG[3] = {0}, Vcont[3] = {0}, Vtens[3] = {0};
        double Vsovii[3] = {0}, Vsovjj[3] = {0}, Vsovji[3] = {0}, Vsovij[3] = {0};
        double Vstring[3] = {0}, Vsosii[3] = {0}, Vsosjj[3] = {0};
        double pogeG[3] = {0}, pcont[3] = {0}, ptens[3] = {0};
        double psovii[3] = {0}, psovjj[3] = {0}, psovji[3] = {0}, psovij[3] = {0};
        double psosii[3] = {0}, psosjj[3] = {0}, T[3] = {0}, rmsr[3] = {0}, rmsl[3] = {0};

        for (int nfp = 0; nfp < job->qnlist_full.len_part[nf]; nfp++) {
            for (int nip = 0; nip < job->qnlist_full.len_part[ni]; nip++) {
                int mapf1 = job->qnlist_full.qnum[nf][nfp].map1;
                int mapf2 = job->qnlist_full.qnum[nf][nfp].map2;
                int mapi1 = job->qnlist_full.qnum[ni][nip].map1;
                int mapi2 = job->qnlist_full.qnum[ni][nip].map2;
                basis_qnum qf = job->qnlist_full.qnum[nf][nfp];
                basis_qnum qi = job->qnlist_full.qnum[ni][nip];

                nfi += me_add(tir_cent, job->scdk[0][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteNfi, 1);

                for (int p = 0; p < 3; p++) {
                    int c = p + 1;
                    VogeG[p] += me_add(tir_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVogeG, c);
                    Vcont[p] += me_add(tir_cent, job->scdk[3 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVcont, c);
                    Vtens[p] += me_add(tir_tens, job->scdk[6 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVtens, c);
                    Vsovii[p] += me_add(tir_soii, job->scdk[9 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsovii, c);
                    Vsovjj[p] += me_add(tjr_sojj, job->scdk[12 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsovjj, c);
                    Vsovji[p] += me_add(t1r_soji, job->scdk[15 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsovji, c);
                    Vsovij[p] += me_add(tir_soij, job->scdk[18 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsovij, c);
                    Vstring[p] += me_add(tir_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVstring, c);
                    Vsosii[p] += me_add(tir_soii, job->scdk[9 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsosii, c);
                    Vsosjj[p] += me_add(tjr_sojj, job->scdk[12 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteVsosjj, c);
                    pogeG[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepogeG, c);
                    pcont[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepcont, c);
                    ptens[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteptens, c);
                    psovii[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsovii, c);
                    psovjj[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsovjj, c);
                    psovji[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsovji, c);
                    psovij[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsovij, c);
                    psosii[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsosii, c);
                    psosjj[p] += me_add(t1p_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, intepsosjj, c);
                    T[p] += me_add(tpi_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteTi, c);
                    rmsr[p] += me_add(tir_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteRMS, c);
                    rmsl[p] += me_add(t2r_cent, job->scdk[0 + p][mapf1][mapf2][mapi1][mapi2], qf, qi, job->varg, inteRMS, c);
                }
            }
        }

        job->Nfi.value[nf][ni] = nfi;
        for (int p = 0; p < 3; p++) {
            job->VogeG[p].value[nf][ni] = VogeG[p];
            job->Vcont[p].value[nf][ni] = Vcont[p];
            job->Vtens[p].value[nf][ni] = Vtens[p];
            job->Vsovii[p].value[nf][ni] = Vsovii[p];
            job->Vsovjj[p].value[nf][ni] = Vsovjj[p];
            job->Vsovji[p].value[nf][ni] = Vsovji[p];
            job->Vsovij[p].value[nf][ni] = Vsovij[p];
            job->Vstring[p].value[nf][ni] = Vstring[p];
            job->Vsosii[p].value[nf][ni] = Vsosii[p];
            job->Vsosjj[p].value[nf][ni] = Vsosjj[p];
            job->pogeG[p].value[nf][ni] = pogeG[p];
            job->pcont[p].value[nf][ni] = pcont[p];
            job->ptens[p].value[nf][ni] = ptens[p];
            job->psovii[p].value[nf][ni] = psovii[p];
            job->psovjj[p].value[nf][ni] = psovjj[p];
            job->psovji[p].value[nf][ni] = psovji[p];
            job->psovij[p].value[nf][ni] = psovij[p];
            job->psosii[p].value[nf][ni] = psosii[p];
            job->psosjj[p].value[nf][ni] = psosjj[p];
            job->T[p].value[nf][ni] = T[p];
            job->rmsr[p].value[nf][ni] = rmsr[p];
            job->rmsl[p].value[nf][ni] = rmsl[p];
        }
    }
}

void spectra_baryon_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out,
    array_t *rms12, array_t *rms13, array_t *rms23)
{
    baryon_job_t job = {0};

    scdk_vargs_from_model(&job.varg, args_model, args_dynmc->model);
    baryon_basis_build(args_input, args_model, &job.qnlist_spfy, &job.qnlist_full);

    int nbas = job.qnlist_full.len_list;
    if (nbas <= 0) {
        fprintf(stderr, "Error: empty baryon basis. Check J, P, sym12, and Lmax.\n");
        exit(1);
    }

    baryon_mlsj_jl(&job);

    int nsp = job.qnlist_spfy.len_list;
    for (int op = 0; op < 21; op++) {
        threebody_scdk_table_alloc(job.qnlist_spfy.len_part, nsp, &job.scdk[op]);
    }

    thread_load(calc_scdk_mt, &job, 21 * nsp, getNumCores());

    matrix_t *mats[80];
    int nmats = baryon_collect_mats(&job, mats);
    for (int i = 0; i < nmats; i++) {
        *mats[i] = matrix_init(nbas, nbas);
    }

    thread_load(getmfi, &job, nbas, getNumCores());

    for (int i = 0; i < nmats; i++) {
        matrix_symmetrize(mats[i]);
    }

    matrix_t vt;
    int nkeep = threebody_overlap_basis(&job.Nfi, 1e-8, &vt);
    if (nkeep <= 0) {
        fprintf(stderr, "Error: baryon overlap N has no positive eigenvalues.\n");
        exit(1);
    }

    matrix_t temp = matrix_init(nkeep, nkeep);
    matrix_t tV = matrix_init(nkeep, nkeep);
    matrix_t tP = matrix_init(nkeep, nkeep);

    /* Recycle eigsys: each pair is mapped with SCDK, then βVβ per pair. */
    matrix_t Hfi = matrix_init(nkeep, nkeep);
    for (int p = 0; p < 3; p++) {
        add_reduced(&Hfi, &vt, &job.T[p], &temp);
        add_reduced(&Hfi, &vt, &job.Vstring[p], &temp);
        add_sandwiched(&Hfi, &vt, &job.VogeG[p], &job.pogeG[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vcont[p], &job.pcont[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vtens[p], &job.ptens[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsovii[p], &job.psovii[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsovjj[p], &job.psovjj[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsovji[p], &job.psovji[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsovij[p], &job.psovij[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsosii[p], &job.psosii[p], &tV, &tP, &temp);
        add_sandwiched(&Hfi, &vt, &job.Vsosjj[p], &job.psosjj[p], &tV, &tP, &temp);
    }

    *e_out = array_init(nkeep);
    matrix_t ut = matrix_init(nkeep, nkeep);

#ifdef LAPACKE
    lapack_standard(Hfi.value, nkeep, e_out->value, ut.value, nkeep);
#else
    eigen_standard(Hfi.value, nkeep, e_out->value, ut.value, nkeep);
#endif

    *v_out = matrix_init(nkeep, nbas);
    matrix_product(&ut, &vt, v_out);

    *n_out = matrix_init(nbas, nbas);
    matrix_copy(n_out, &job.Nfi);

    *rms12 = array_init(nkeep);
    *rms13 = array_init(nkeep);
    *rms23 = array_init(nkeep);
    for (int n = 0; n < nkeep; n++) {
        rms12->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[0]))) / GEMSTORE_FM;
        rms13->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[1]))) / GEMSTORE_FM;
        rms23->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[2]))) / GEMSTORE_FM;
    }

    if (args_input->print_wfn) {
        char path[280];
        sprintf(path, "%s.basis.dat", args_input->project);
        FILE *bf = fopen(path, "w");
        if (bf) {
            fprintf(bf, "# i  c  nrho nlam lrho llam L  sij  jl  coe\n");
            for (int i = 0; i < nbas; i++) {
                basis_qnum q = job.qnlist_full.qnum[i][0];
                fprintf(bf, "%d  %d  %d %d  %d %d %d  %.1f  %.1f  %.6f\n",
                    i, q.c, q.nrho, q.nlam, q.lrho, q.llam, q.L, q.sij, q.jl, q.coe);
            }
            fclose(bf);
        }
        for (int n = 0; n < nkeep; n++) {
            sprintf(path, "%s.wfn.%d.dat", args_input->project, n + 1);
            FILE *wf = fopen(path, "w");
            if (!wf) {
                continue;
            }
            fprintf(wf, "# i  coefficient   (see %s.basis.dat)\n", args_input->project);
            for (int i = 0; i < nbas; i++) {
                fprintf(wf, "%d  %.16e\n", i, v_out->value[n][i]);
            }
            fclose(wf);
        }
    }

    for (int op = 0; op < 21; op++) {
        threebody_scdk_table_free(job.qnlist_spfy.len_part, nsp, &job.scdk[op]);
    }
    for (int i = 0; i < nsp; i++) {
        for (int j = 0; j < job.qnlist_spfy.len_part[i]; j++) {
            matrix_free(&job.mlsj[i][j]);
        }
        free(job.mlsj[i]);
    }
    free(job.mlsj);
    basis_list_free(&job.qnlist_spfy);
    basis_list_free(&job.qnlist_full);

    matrix_free(&temp);
    matrix_free(&tV);
    matrix_free(&tP);
    matrix_free(&vt);
    matrix_free(&ut);
    matrix_free(&Hfi);
    for (int i = 0; i < nmats; i++) {
        matrix_free(mats[i]);
    }
}
