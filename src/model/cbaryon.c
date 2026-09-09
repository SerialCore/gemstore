/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Baryon SPECTRA on a Jacobi GEM basis, structured like cmeson.c.
 *
 * Wave functions live in a single Jacobi frame (c=1: ρ = r1−r2). Pair
 * operators still loop over pair = 1,2,3. Same-channel spatial MEs use the
 * meson 1D GRnlr/GRnlp integrals; off-channel central MEs map both Gaussians
 * onto the pair frame (complete-the-square + solid-harmonic addition).
 *
 * GIVt multiplies OCent: same-channel OCent must be a full q-number Kronecker,
 * otherwise the 1D radial path mixes (lρ,lλ,jl) blocks and H acquires a kernel.
 */

#include <gemstore/model/cbaryon.h>
#include <gemstore/model/gimodel.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/jacobi.h>
#include <gemstore/basis/orbit.h>

#include <gemstore/math/soc.h>
#include <gemstore/math/solidharm.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>

#include <gemstore/param/argset.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Angular Kronecker in one Jacobi frame. N and same-channel central/T use this;
 * do not replace it with recoupled meson operator_center_sl. */
static int baryon_qn_match(const basis_qnum *a, const basis_qnum *b)
{
    return a->c == b->c
        && a->lrho == b->lrho && a->llam == b->llam && a->L == b->L
        && a->sij == b->sij && a->jl == b->jl && a->J == b->J;
}

/* True if the pair that defines Jacobi channel c is flavour-identical. */
static int pair_identical(int f1, int f2, int f3, int c)
{
    switch (c) {
        case 1: return f1 == f2;
        case 2: return f3 == f1;
        case 3: return f2 == f3;
        default: return 0;
    }
}

/* |(lρ lλ)L, sij; jl⟩ → |(sij lρ)jρ, lλ; jl⟩ so meson pair operators apply. */
static double recouple_L_to_lrho(double sij, int lrho, int llam, int L, double jl, double jrho)
{
    double phase = pow(-1.0, sij + llam + jl + L);
    double hat = sqrt((2.0 * L + 1.0) * (2.0 * jrho + 1.0));

    return phase * hat * sixJ_symbol(lrho, llam, L, jl, sij, jrho);
}

/* |(L sij)jl, s3; J⟩ → |(sij s3)S, L; J⟩. Used for off-channel s_i·s_j. */
static double recouple_jl_to_S(double sij, int L, double jl, double s3, double J, double S)
{
    double phase = pow(-1.0, sij + L + s3 + J);
    double hat = sqrt((2.0 * jl + 1.0) * (2.0 * S + 1.0));

    return phase * hat * sixJ_symbol(sij, L, jl, s3, J, S);
}

/* Spin overlap between pair-spin sija in channel ca and sijb in channel cb
 * at total spin S. Adjacent channels (1↔3, 2↔3) pick up (−1)^{3/2+S}. */
static double recouple_spin_pair(int ca, double sija, int cb, double sijb, double S)
{
    if (ca == cb) {
        return (sija == sijb) ? 1.0 : 0.0;
    }

    double hat = sqrt((2.0 * sija + 1.0) * (2.0 * sijb + 1.0));
    double sixj = sixJ_symbol(0.5, 0.5, sija, 0.5, S, sijb);

    if ((ca == 1 && cb == 3) || (ca == 3 && cb == 1) ||
        (ca == 2 && cb == 3) || (ca == 3 && cb == 2)) {
        return pow(-1.0, 1.5 + S) * hat * sixj;
    }

    return pow(-1.0, sija + sijb) * hat * sixj;
}

/* Apply a meson (si sj s l) operator on the ρ pair, same Jacobi channel. */
static double baryon_op_apply(operator_sl osl, const basis_qnum *a, const basis_qnum *b)
{
    if (a->c != b->c || a->llam != b->llam || a->J != b->J) {
        return 0.0;
    }

    double sum = 0.0;
    double jmin = fabs(a->sij - a->lrho);
    double jmax = a->sij + a->lrho;
    double jminp = fabs(b->sij - b->lrho);
    double jmaxp = b->sij + b->lrho;

    if (jminp > jmin) jmin = jminp;
    if (jmaxp < jmax) jmax = jmaxp;

    for (double jrho = jmin; jrho <= jmax + 1e-9; jrho += 1.0) {
        double rec_a = recouple_L_to_lrho(a->sij, a->lrho, a->llam, a->L, a->jl, jrho);
        double rec_b = recouple_L_to_lrho(b->sij, b->lrho, b->llam, b->L, b->jl, jrho);
        sum += rec_a * rec_b * osl(0.5, 0.5, a->sij, a->lrho, 0.5, 0.5, b->sij, b->lrho, jrho);
    }

    return sum;
}

/* ⟨s_i·s_j⟩ on quark pair `pair`. Off-channel: recouple to total S then to
 * the pair spin of that pair. */
static double baryon_op_sdots_pair(const basis_qnum *a, const basis_qnum *b, int pair)
{
    if (a->L != b->L || a->J != b->J) {
        return 0.0;
    }
    if (pair == a->c) {
        return baryon_op_apply(operator_sdots_sl, a, b);
    }

    double sum = 0.0;
    for (double S = 0.5; S <= 1.5 + 1e-9; S += 1.0) {
        double ra = recouple_jl_to_S(a->sij, a->L, a->jl, 0.5, a->J, S);
        double rb = recouple_jl_to_S(b->sij, b->L, b->jl, 0.5, b->J, S);
        for (double sijp = 0.0; sijp <= 1.0 + 1e-9; sijp += 1.0) {
            double ca = recouple_spin_pair(a->c, a->sij, pair, sijp, S);
            double cb = recouple_spin_pair(b->c, b->sij, pair, sijp, S);
            double sdots = 0.5 * (sijp * (sijp + 1.0) - 1.5);
            sum += ra * rb * ca * cb * sdots;
        }
    }
    return sum;
}

static void baryon_set_operators(argsGIModelDy_t *dyn, const basis_qnum *a, const basis_qnum *b, int pair)
{
    /* Same-channel central/GIVt: Kronecker on the full baryon q-numbers so
     * the 1D GRnlr path cannot mix different (lρ,lλ,L,sij,jl). Off-channel
     * central may connect (1,0)↔(0,1) at fixed L; the spatial ME is RR. */
    if (pair == a->c && pair == b->c) {
        dyn->OCent = baryon_qn_match(a, b) ? 1.0 : 0.0;
        dyn->OSdS  = baryon_op_apply(operator_sdots_sl, a, b);
        dyn->OLSi  = baryon_op_apply(operator_ldotsi_sl, a, b);
        dyn->OLSj  = baryon_op_apply(operator_ldotsj_sl, a, b);
        dyn->OTens = baryon_op_apply(operator_tensor_sl, a, b);
        return;
    }

    if (a->J == b->J && a->L == b->L && a->jl == b->jl && a->sij == b->sij) {
        dyn->OCent = 1.0;
    }
    else {
        dyn->OCent = 0.0;
    }
    dyn->OSdS  = baryon_op_sdots_pair(a, b, pair);
    dyn->OLSi  = 0.0;   /* off-channel L·S / tensor not implemented */
    dyn->OLSj  = 0.0;
    dyn->OTens = 0.0;
}

/* GRnlr without the r^l factor; solid harmonics restore |x|^l Y_lm. */
static double gem_pref(int l, double nu)
{
    return pow(2.0, l / 2.0 + 1.25) * pow(nu, l / 2.0 + 0.75) / sqrt(tgamma(l + 1.5));
}

static double pot_one(double r, void *ctx)
{
    (void)r;
    (void)ctx;
    return 1.0;
}

/* Central ME of pot(|r_pair|) after mapping both Gaussians onto the pair frame.
 * Do not feed b11 into integral_nlr_hamilton: that double-counts GRnlr ν^{3/4}. */
static double me_pair_reduced_r(potential_t pot, gi_pot_ctx_t *ctx,
    const basis_qnum *qa, const basis_qnum *qb, int pair)
{
    jacobi_shift_t sh;
    double pref;
    int M;

    if (!jacobi_gaussian_shift(qa->m1, qa->m2, qa->m3,
            qa->c, qa->nurho, qa->nulam,
            qb->c, qb->nurho, qb->nulam,
            pair, &sh)) {
        return 0.0;
    }

    pref = gem_pref(qa->lrho, qa->nurho) * gem_pref(qa->llam, qa->nulam)
        * gem_pref(qb->lrho, qb->nurho) * gem_pref(qb->llam, qb->nulam);
    M = qa->L;  /* scalar ME is M-independent; stretched M=L is enough */

    return pref * solidharm_central_me(
        qa->lrho, qa->llam, qa->L, M,
        qb->lrho, qb->llam, qb->L,
        sh.al_a, sh.be_a, sh.ga_a, sh.de_a,
        sh.al_b, sh.be_b, sh.ga_b, sh.de_b,
        sh.b11, sh.aRR, pot, ctx);
}

/* Spatial ME of a pair operator. momentum=1 is GI β(p),δ(p): same-channel
 * 1D GRnlp only; off-channel smearing is not implemented. */
static double me_pair_spatial(potential_t pot, gi_pot_ctx_t *ctx,
    const basis_qnum *qa, const basis_qnum *qb,
    const argsOrbit_t *rho_a, const argsOrbit_t *rho_b,
    const argsOrbit_t *lam_a, const argsOrbit_t *lam_b,
    int pair, int momentum)
{
    if (qa->c == pair && qb->c == pair) {
        double overlap_lam = integral_nlr_overlap(GRnlr,
            1.0 / sqrt(lam_a->scale + lam_b->scale), lam_a, lam_b);
        if (momentum) {
            double factor_p = sqrt(4.0 * rho_a->scale * rho_b->scale / (rho_a->scale + rho_b->scale));
            return integral_nlp_hamilton(GRnlp, pot, factor_p, rho_a, rho_b, ctx) * overlap_lam;
        }
        double factor_r = 1.0 / sqrt(rho_a->scale + rho_b->scale);
        return integral_nlr_hamilton(GRnlr, pot, factor_r, rho_a, rho_b, ctx) * overlap_lam;
    }

    if (momentum) {
        return 0.0;
    }
    return me_pair_reduced_r(pot, ctx, qa, qb, pair);
}

/* Gold checks: same-channel Coulomb 1D vs RR, and V=1 frame independence. */
static void baryon_check_reduce_identity(const basis_qnum *qn, const argsOrbit_t *rho, const argsOrbit_t *lam,
    int nbas, gi_pot_ctx_t *ctx)
{
    int ncheck = 0;
    double maxrel = 0.0;

    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            if (qn[i].c != qn[j].c) {
                continue;
            }
            int pair = qn[i].c;
            double direct = me_pair_spatial(GIVcoul, ctx, &qn[i], &qn[j],
                &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
            double reduced = me_pair_reduced_r(GIVcoul, ctx, &qn[i], &qn[j], pair);
            double denom = fabs(direct) + 1e-15;
            double rel = fabs(direct - reduced) / denom;
            if (rel > maxrel) {
                maxrel = rel;
            }

            /* V=1 is the overlap: independent of which pair frame we reduce in. */
            double u1 = me_pair_reduced_r(pot_one, ctx, &qn[i], &qn[j], pair);
            double u2 = me_pair_reduced_r(pot_one, ctx, &qn[i], &qn[j], pair == 1 ? 2 : 1);
            double urel = fabs(u1 - u2) / (fabs(u1) + 1e-15);
            if (urel > maxrel) {
                maxrel = urel;
            }

            ncheck++;
            if (ncheck >= 16) {
                break;
            }
        }
        if (ncheck >= 16) {
            break;
        }
    }

    if (ncheck == 0) {
        return;
    }
    if (maxrel > 1e-8) {
        fprintf(stderr, "Error: Jacobi reduce identity failed, maxrel=%.3e\n", maxrel);
        exit(1);
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

    /* One Jacobi frame is complete. Do not open c=2,3 until N and T have
     * cross-channel blocks (otherwise H is a false direct sum). */
    for (int c = 1; c <= 1; c++) {
        for (int lrho = 0; lrho <= Lmax; lrho++) {
            for (int llam = 0; llam <= Lmax - lrho; llam++) {
                /* P = (−1)^{lρ+lλ} */
                if (((lrho + llam) % 2 == 0) ? (P != 1) : (P != -1)) {
                    continue;
                }

                for (int L = abs(lrho - llam); L <= lrho + llam; L++) {
                    for (double sij = 0.0; sij <= 1.0 + 1e-9; sij += 1.0) {
                        if (pair_identical(input->f1, input->f2, input->f3, c)) {
                            /* (−1)^{sij+lρ} for the identical pair in this channel */
                            int phase = ((1 + (int)sij + lrho) % 2 == 0) ? 1 : -1;
                            if (f12 * phase != 1) {
                                continue;
                            }
                        }

                        for (double jl = fabs(L - sij); jl <= L + sij + 1e-9; jl += 1.0) {
                            if (fabs(jl - s3) - 1e-9 > J || jl + s3 + 1e-9 < J) {
                                continue;
                            }

                            basis_list_push(spfy, 1, -1, -1, 1.0,
                                m1, m2, m3, s1, s2, s3,
                                0.0, 0.0, 0.0, 0.0, 0.0,
                                c, lrho, llam, L, sij, jl, J,
                                0, 0, 0.0, 0.0);
                        }
                    }
                }
            }
        }
    }

    basis_list_push_full(spfy, full, input->rmin, input->rmax, input->nmax);
}

void spectra_baryon_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out)
{
    basis_list qnlist_spfy;
    basis_list qnlist_full;

    baryon_basis_build(args_input, args_model, &qnlist_spfy, &qnlist_full);

    int nbas = qnlist_full.len_list;
    if (nbas <= 0) {
        fprintf(stderr, "Error: empty baryon basis. Check J, P, sym12, and Lmax.\n");
        exit(1);
    }

    basis_qnum *qn = (basis_qnum *)malloc((size_t)nbas * sizeof(basis_qnum));
    argsOrbit_t *rho = (argsOrbit_t *)malloc((size_t)nbas * sizeof(argsOrbit_t));
    argsOrbit_t *lam = (argsOrbit_t *)malloc((size_t)nbas * sizeof(argsOrbit_t));
    for (int i = 0; i < nbas; i++) {
        qn[i] = qnlist_full.qnum[i][0];
        rho[i].n = qn[i].nrho;
        rho[i].l = qn[i].lrho;
        rho[i].scale = qn[i].nurho;
        rho[i].param = 0.0;
        lam[i].n = qn[i].nlam;
        lam[i].l = qn[i].llam;
        lam[i].scale = qn[i].nulam;
        lam[i].param = 0.0;
    }

    /* construct matrices */
    matrix_t mT = matrix_init(nbas, nbas);
    matrix_t mbetaijCoul = matrix_init(nbas, nbas);
    matrix_t mdeltaijCont = matrix_init(nbas, nbas);
    matrix_t mdeltaiiSov = matrix_init(nbas, nbas);
    matrix_t mdeltajjSov = matrix_init(nbas, nbas);
    matrix_t mdeltaijSov = matrix_init(nbas, nbas);
    matrix_t mdeltaiiSos = matrix_init(nbas, nbas);
    matrix_t mdeltajjSos = matrix_init(nbas, nbas);
    matrix_t mdeltaijTens = matrix_init(nbas, nbas);
    matrix_t mVcoul = matrix_init(nbas, nbas);
    matrix_t mVconf = matrix_init(nbas, nbas);
    matrix_t mVcont = matrix_init(nbas, nbas);
    matrix_t mVsovi = matrix_init(nbas, nbas);
    matrix_t mVsovj = matrix_init(nbas, nbas);
    matrix_t mVsovij = matrix_init(nbas, nbas);
    matrix_t mVsosi = matrix_init(nbas, nbas);
    matrix_t mVsosj = matrix_init(nbas, nbas);
    matrix_t mVtens = matrix_init(nbas, nbas);
    matrix_t tmT = matrix_init(nbas, nbas);
    matrix_t tmbetaijCoul = matrix_init(nbas, nbas);
    matrix_t tmdeltaijCont = matrix_init(nbas, nbas);
    matrix_t tmdeltaiiSov = matrix_init(nbas, nbas);
    matrix_t tmdeltajjSov = matrix_init(nbas, nbas);
    matrix_t tmdeltaijSov = matrix_init(nbas, nbas);
    matrix_t tmdeltaiiSos = matrix_init(nbas, nbas);
    matrix_t tmdeltajjSos = matrix_init(nbas, nbas);
    matrix_t tmdeltaijTens = matrix_init(nbas, nbas);
    matrix_t tmVcoul = matrix_init(nbas, nbas);
    matrix_t tmVconf = matrix_init(nbas, nbas);
    matrix_t tmVcont = matrix_init(nbas, nbas);
    matrix_t tmVsovi = matrix_init(nbas, nbas);
    matrix_t tmVsovj = matrix_init(nbas, nbas);
    matrix_t tmVsovij = matrix_init(nbas, nbas);
    matrix_t tmVsosi = matrix_init(nbas, nbas);
    matrix_t tmVsosj = matrix_init(nbas, nbas);
    matrix_t tmVtens = matrix_init(nbas, nbas);
    matrix_t Hfi = matrix_init(nbas, nbas);
    matrix_t Nfi = matrix_init(nbas, nbas);

    double Cij = -2.0 / 3.0;  /* baryon pair colour; meson is −4/3 */
    gi_pot_ctx_t ctx = { args_model, args_dynmc };

    {
        double mi, mj, mk;
        jacobi_pair_mass(qn[0].m1, qn[0].m2, qn[0].m3, qn[0].c, &mi, &mj, &mk);
        args_dynmc->Cij = Cij;
        args_dynmc->OCent = 1.0;
        args_dynmc->mi = mi;
        args_dynmc->mj = mj;
        args_dynmc->Sigij = sigma_ij(mi, mj, args_model->sigma_0, args_model->s);
        sigma_k_ij(args_dynmc->Sigij, args_dynmc->Sigkij);
        baryon_check_reduce_identity(qn, rho, lam, nbas, &ctx);
    }

    /* calculate matrix elements */
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            /* matrix_init uses malloc; pair loop accumulates with += */
            mbetaijCoul.value[i][j] = 0.0;
            mdeltaijCont.value[i][j] = 0.0;
            mdeltaiiSov.value[i][j] = 0.0;
            mdeltajjSov.value[i][j] = 0.0;
            mdeltaijSov.value[i][j] = 0.0;
            mdeltaiiSos.value[i][j] = 0.0;
            mdeltajjSos.value[i][j] = 0.0;
            mdeltaijTens.value[i][j] = 0.0;
            mVcoul.value[i][j] = 0.0;
            mVconf.value[i][j] = 0.0;
            mVcont.value[i][j] = 0.0;
            mVsovi.value[i][j] = 0.0;
            mVsovj.value[i][j] = 0.0;
            mVsovij.value[i][j] = 0.0;
            mVsosi.value[i][j] = 0.0;
            mVsosj.value[i][j] = 0.0;
            mVtens.value[i][j] = 0.0;

            double factor_r = 1.0 / sqrt(rho[i].scale + rho[j].scale);
            double factor_p = sqrt(4.0 * rho[i].scale * rho[j].scale / (rho[i].scale + rho[j].scale));
            double factor_rl = 1.0 / sqrt(lam[i].scale + lam[j].scale);
            double factor_pl = sqrt(4.0 * lam[i].scale * lam[j].scale / (lam[i].scale + lam[j].scale));
            double overlap_lam = integral_nlr_overlap(GRnlr, factor_rl, &lam[i], &lam[j]);
            double overlap_rho = integral_nlr_overlap(GRnlr, factor_r, &rho[i], &rho[j]);
            int c = qn[i].c;
            double mi, mj, mk;

            /* 1D Iρ Iλ is valid only in the same Jacobi frame. */
            Nfi.value[i][j] = overlap_rho * overlap_lam * (baryon_qn_match(&qn[i], &qn[j]) ? 1.0 : 0.0);

            jacobi_pair_mass(qn[i].m1, qn[i].m2, qn[i].m3, c, &mi, &mj, &mk);
            args_dynmc->Cij = Cij;
            args_dynmc->mi = mi;
            args_dynmc->mj = mj;
            args_dynmc->Sigij = sigma_ij(mi, mj, args_model->sigma_0, args_model->s);
            sigma_k_ij(args_dynmc->Sigij, args_dynmc->Sigkij);
            baryon_set_operators(args_dynmc, &qn[i], &qn[j], c);

            /* T = Σ_k √(m_k²+p_k²): pair quarks on ρ, spectator on λ. */
            mT.value[i][j] = integral_nlp_hamilton(GRnlp, GIVt, factor_p, &rho[i], &rho[j], &ctx) * overlap_lam;
            args_dynmc->mi = mk;
            args_dynmc->mj = mk;
            mT.value[i][j] += integral_nlp_hamilton(GRnlp, GIVt_quark, factor_pl, &lam[i], &lam[j], &ctx) * overlap_rho;

            for (int pair = 1; pair <= 3; pair++) {
                jacobi_pair_mass(qn[i].m1, qn[i].m2, qn[i].m3, pair, &mi, &mj, &mk);
                args_dynmc->mi = mi;
                args_dynmc->mj = mj;
                args_dynmc->Sigij = sigma_ij(mi, mj, args_model->sigma_0, args_model->s);
                sigma_k_ij(args_dynmc->Sigij, args_dynmc->Sigkij);
                baryon_set_operators(args_dynmc, &qn[i], &qn[j], pair);

                mbetaijCoul.value[i][j] += me_pair_spatial(GIVbetaijcoul, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltaijCont.value[i][j] += me_pair_spatial(GIVdeltaijcont, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltaiiSov.value[i][j] += me_pair_spatial(GIVdeltaiisov, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltajjSov.value[i][j] += me_pair_spatial(GIVdeltajjsov, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltaijSov.value[i][j] += me_pair_spatial(GIVdeltaijsov, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltaiiSos.value[i][j] += me_pair_spatial(GIVdeltaiisos, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltajjSos.value[i][j] += me_pair_spatial(GIVdeltajjsos, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mdeltaijTens.value[i][j] += me_pair_spatial(GIVdeltaijtens, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 1);
                mVcoul.value[i][j] += me_pair_spatial(GIVcoul, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVconf.value[i][j] += me_pair_spatial(GIVconf, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVcont.value[i][j] += me_pair_spatial(GIVcont, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVsovi.value[i][j] += me_pair_spatial(GIVsovi, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVsovj.value[i][j] += me_pair_spatial(GIVsovj, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVsovij.value[i][j] += me_pair_spatial(GIVsovij, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVsosi.value[i][j] += me_pair_spatial(GIVsosi, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVsosj.value[i][j] += me_pair_spatial(GIVsosj, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
                mVtens.value[i][j] += me_pair_spatial(GIVtens, &ctx, &qn[i], &qn[j], &rho[i], &rho[j], &lam[i], &lam[j], pair, 0);
            }
        }
    }

    /* prepare a random symmetric matrix */
    matrix_t rand = matrix_random(nbas, nbas);
    matrix_t temp = matrix_init(nbas, nbas);
    matrix_transpose(&rand, &temp);
    matrix_sum(&rand, &temp, &rand);

    matrix_t vt = matrix_init(nbas, nbas);

    *e_out = array_init(nbas);

#ifdef LAPACKE
    lapack_general(rand.value, Nfi.value, nbas, e_out->value, vt.value, nbas);
#else
    eigen_general(rand.value, Nfi.value, nbas, e_out->value, vt.value, nbas);
#endif

    /* construct new orthogonal basis */
    matrix_productT(&vt, &Nfi, &temp);
    for (int k = 0; k < nbas; k++) {
        double norm = sqrt(fabs(temp.value[k][k]));
        if (norm > 1e-10) {
            for (int i = 0; i < nbas; i++) {
                vt.value[k][i] /= norm;
            }
        }
        else {
            printf("Warning: singular vector %d, norm=%.2e\n", k, norm);
        }
    }

    /* transform Hamiltonian matrices in orthogonal basis */
    matrix_productT(&vt, &mT, &tmT);
    matrix_productT(&vt, &mbetaijCoul, &tmbetaijCoul);
    matrix_productT(&vt, &mdeltaijCont, &tmdeltaijCont);
    matrix_productT(&vt, &mdeltaiiSov, &tmdeltaiiSov);
    matrix_productT(&vt, &mdeltajjSov, &tmdeltajjSov);
    matrix_productT(&vt, &mdeltaijSov, &tmdeltaijSov);
    matrix_productT(&vt, &mdeltaiiSos, &tmdeltaiiSos);
    matrix_productT(&vt, &mdeltajjSos, &tmdeltajjSos);
    matrix_productT(&vt, &mdeltaijTens, &tmdeltaijTens);
    matrix_productT(&vt, &mVcoul, &tmVcoul);
    matrix_productT(&vt, &mVconf, &tmVconf);
    matrix_productT(&vt, &mVcont, &tmVcont);
    matrix_productT(&vt, &mVsovi, &tmVsovi);
    matrix_productT(&vt, &mVsovj, &tmVsovj);
    matrix_productT(&vt, &mVsovij, &tmVsovij);
    matrix_productT(&vt, &mVsosi, &tmVsosi);
    matrix_productT(&vt, &mVsosj, &tmVsosj);
    matrix_productT(&vt, &mVtens, &tmVtens);

    /* Meson-style PVP: H = T + Vconf + β V β + δ V δ + … .
     * Pair sums are taken first, so β of one pair can sandwich V of another. */
    matrix_sum(&tmT, &tmVconf, &Hfi);
    matrix_productT(&tmbetaijCoul, &tmVcoul, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijCont, &tmVcont, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaiiSov, &tmVsovi, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltajjSov, &tmVsovj, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijSov, &tmVsovij, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaiiSos, &tmVsosi, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltajjSos, &tmVsosj, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);
    matrix_productT(&tmdeltaijTens, &tmVtens, &temp);
    matrix_sum(&Hfi, &temp, &Hfi);

    /* final eigen system */
    matrix_t ut = matrix_init(nbas, nbas);

#ifdef LAPACKE
    lapack_standard(Hfi.value, nbas, e_out->value, ut.value, nbas);
#else
    eigen_standard(Hfi.value, nbas, e_out->value, ut.value, nbas);
#endif

    *v_out = matrix_init(nbas, nbas);
    matrix_product(&ut, &vt, v_out);

    *n_out = matrix_init(nbas, nbas);
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            n_out->value[i][j] = Nfi.value[i][j];
        }
    }

    free(qn);
    free(rho);
    free(lam);
    basis_list_free(&qnlist_spfy);
    basis_list_free(&qnlist_full);
    matrix_free(&temp);
    matrix_free(&rand);
    matrix_free(&vt);
    matrix_free(&ut);
    matrix_free(&mT);
    matrix_free(&mbetaijCoul);
    matrix_free(&mdeltaijCont);
    matrix_free(&mdeltaiiSov);
    matrix_free(&mdeltajjSov);
    matrix_free(&mdeltaijSov);
    matrix_free(&mdeltaiiSos);
    matrix_free(&mdeltajjSos);
    matrix_free(&mdeltaijTens);
    matrix_free(&mVcoul);
    matrix_free(&mVconf);
    matrix_free(&mVcont);
    matrix_free(&mVsovi);
    matrix_free(&mVsovj);
    matrix_free(&mVsovij);
    matrix_free(&mVsosi);
    matrix_free(&mVsosj);
    matrix_free(&mVtens);
    matrix_free(&tmT);
    matrix_free(&tmbetaijCoul);
    matrix_free(&tmdeltaijCont);
    matrix_free(&tmdeltaiiSov);
    matrix_free(&tmdeltajjSov);
    matrix_free(&tmdeltaijSov);
    matrix_free(&tmdeltaiiSos);
    matrix_free(&tmdeltajjSos);
    matrix_free(&tmdeltaijTens);
    matrix_free(&tmVcoul);
    matrix_free(&tmVconf);
    matrix_free(&tmVcont);
    matrix_free(&tmVsovi);
    matrix_free(&tmVsovj);
    matrix_free(&tmVsovij);
    matrix_free(&tmVsosi);
    matrix_free(&tmVsosj);
    matrix_free(&tmVtens);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}
