/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2026, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/basis.h>
#include <gemstore/basis/orbit.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void basis_list_init(basis_list *qnlist)
{
    qnlist->qnum = (argsBasis_t **)malloc(sizeof(argsBasis_t *) * 0);
    qnlist->len_part = (int *)malloc(sizeof(int) * 0);
    qnlist->len_list = 0;
}

void basis_list_logs(basis_list qnlist)
{
    int i, j;
    printf("-------------------basis_list start--------------------\n");
    for (i = 0; i < qnlist.len_list; i++)
    {
        printf("%3d: |", i);
        for (j = 0; j < qnlist.len_part[i]; j++)
        {
            printf("j=%d map1=%3d map2=%3d coe=%10.6f c=%d nrho=%2d nlam=%2d |", j, qnlist.qnum[i][j].map1, qnlist.qnum[i][j].map2, qnlist.qnum[i][j].coe, qnlist.qnum[i][j].c, qnlist.qnum[i][j].nrho, qnlist.qnum[i][j].nlam);
        }
        printf("\n");
    }
    printf("-------------------basis_list end---------------------\n\n");
}

void basis_list_push(basis_list *qnlist, int newqnQ, int map1, int map2, double coe, double m1, double m2, double m3, double s1, double s2, double s3, double t1, double t2, double t3, double tij, double T, int c, int lrho, int llam, int L, double sij, double jl, double J, int nrho, int nlam, double nurho, double nulam)
{
    int i, j;
    if (1 == newqnQ)
    {
        i = qnlist->len_list;
        j = 0;
        qnlist->qnum = (argsBasis_t **)realloc(qnlist->qnum, sizeof(argsBasis_t *) * (i + 1));
        qnlist->qnum[i] = (argsBasis_t *)malloc(sizeof(argsBasis_t) * 1);
        qnlist->len_part = (int *)realloc(qnlist->len_part, sizeof(int) * (i + 1));
        qnlist->len_part[i] = 1;
        qnlist->len_list++;
    }
    else
    {
        i = qnlist->len_list - 1;
        j = qnlist->len_part[i];
        qnlist->qnum[i] = (argsBasis_t *)realloc(qnlist->qnum[i], sizeof(argsBasis_t) * (j + 1));
        qnlist->len_part[i]++;
    }

    if (map1 < 0 || map2 < 0)
    {
        qnlist->qnum[i][j].map1 = i;
        qnlist->qnum[i][j].map2 = j;
    }
    else
    {
        qnlist->qnum[i][j].map1 = map1;
        qnlist->qnum[i][j].map2 = map2;
    }

    qnlist->qnum[i][j].coe = coe;
    qnlist->qnum[i][j].m1 = m1;
    qnlist->qnum[i][j].m2 = m2;
    qnlist->qnum[i][j].m3 = m3;
    qnlist->qnum[i][j].s1 = s1;
    qnlist->qnum[i][j].s2 = s2;
    qnlist->qnum[i][j].s3 = s3;
    qnlist->qnum[i][j].t1 = t1;
    qnlist->qnum[i][j].t2 = t2;
    qnlist->qnum[i][j].t3 = t3;
    qnlist->qnum[i][j].tij = tij;
    qnlist->qnum[i][j].T = T;
    qnlist->qnum[i][j].c = c;
    qnlist->qnum[i][j].lrho = lrho;
    qnlist->qnum[i][j].llam = llam;
    qnlist->qnum[i][j].L = L;
    qnlist->qnum[i][j].sij = sij;
    qnlist->qnum[i][j].jl = jl;
    qnlist->qnum[i][j].J = J;
    qnlist->qnum[i][j].nrho = nrho;
    qnlist->qnum[i][j].nlam = nlam;
    qnlist->qnum[i][j].nurho = nurho;
    qnlist->qnum[i][j].nulam = nulam;
}

void basis_list_push_full(basis_list *qnlist_spfy, basis_list *qnlist_full, double rmin, double rmax, int nmax)
{
    int newqnQ;
    for (int i = 0; i < qnlist_spfy->len_list; i++)
    {
        for (int nrho = 1; nrho <= nmax; nrho++)
        {
            for (int nlam = 1; nlam <= nmax; nlam++)
            {
                for (int j = 0; j < qnlist_spfy->len_part[i]; j++)
                {
                    if (0 == j)
                    {
                        newqnQ = 1;
                    }
                    else
                    {
                        newqnQ = 0;
                    }

                    basis_list_push(qnlist_full, newqnQ,
                                    qnlist_spfy->qnum[i][j].map1,
                                    qnlist_spfy->qnum[i][j].map2,
                                    qnlist_spfy->qnum[i][j].coe,
                                    qnlist_spfy->qnum[i][j].m1,
                                    qnlist_spfy->qnum[i][j].m2,
                                    qnlist_spfy->qnum[i][j].m3,
                                    qnlist_spfy->qnum[i][j].s1,
                                    qnlist_spfy->qnum[i][j].s2,
                                    qnlist_spfy->qnum[i][j].s3,
                                    qnlist_spfy->qnum[i][j].t1,
                                    qnlist_spfy->qnum[i][j].t2,
                                    qnlist_spfy->qnum[i][j].t3,
                                    qnlist_spfy->qnum[i][j].tij,
                                    qnlist_spfy->qnum[i][j].T,
                                    qnlist_spfy->qnum[i][j].c,
                                    qnlist_spfy->qnum[i][j].lrho,
                                    qnlist_spfy->qnum[i][j].llam,
                                    qnlist_spfy->qnum[i][j].L,
                                    qnlist_spfy->qnum[i][j].sij,
                                    qnlist_spfy->qnum[i][j].jl,
                                    qnlist_spfy->qnum[i][j].J,
                                    nrho, nlam, getnu(rmin, rmax, nmax, nrho), getnu(rmin, rmax, nmax, nlam));
                }
            }
        }
    }
}

// SL/jj 基构造函数（支持 S-D 混合：L_min to L_max, e.g., 0 to 2）
void build_basis_sd(basis_list *qnlist, double s1, double s2, double J, int P,  // P: parity (-1)^L
                    int L_min, int L_max, int nmax, double rmin, double rmax) {  // 径向参数
    basis_list_init(qnlist);
    
    // 外层：L 循环 (S-D: L=0,2; parity 约束: (-1)^L = P)
    for (int L = L_min; L <= L_max; L += 2) {  // 步长 2 确保 parity 一致 (S-D 同奇偶)
        if (pow(-1, L) != P) continue;  // parity 过滤
        
        // 内层：SL 基 (Mathematica 风格)
        for (double S = fabs(s1 - s2); S <= s1 + s2; S += 1.0) {
            if (fabs(S - L) > J || S + L < J) continue;  // 三角不等式
            
            // jj 子基 (可选，如果用 jj 描述)
            for (double jl = fabs(L - s1); jl <= L + s1; jl += 1.0) {
                if (fabs(jl - s2) > J || jl + s2 < J) continue;
                
                // 径向循环 (Fortran 风格: nrho, nlam)
                for (int nrho = 1; nrho <= nmax; nrho++) {
                    for (int nlam = 1; nlam <= nmax; nlam++) {
                        // 推入基 (类似 basis_list_push)
                        basis_list_push(qnlist, 1, -1, -1, 1.0,  // newqnQ=1, coe=1.0
                                        s1, s2, 0.0,  // m3=0 (meson 无第三质量)
                                        s1, s2, 0.0,  // s3=0
                                        0.0, 0.0, 0.0, 0.0, 0.0,  // t=0 (meson 无味)
                                        0, L, 0, L,  // lrho=L, llam=0 (Jacobi?)
                                        S, jl, J,    // sij=S, jl, J
                                        nrho, nlam,  // 径向量子数
                                        getnu(rmin, rmax, nmax, nrho),  // nurho
                                        getnu(rmin, rmax, nmax, nlam)); // nulam
                    }
                }
            }
        }
    }
    
    // 日志 (可选，类似 basis_list_logs)
    basis_list_logs(*qnlist);
}

// 辅助：计算 parity (-1)^L
static int parity(int L) { return (L % 2 == 0) ? 1 : -1; }

// 构建单一 L 的 SL 和 JJ 子基（严格按 Mathematica）
static void build_sub_basis_for_L(basis_list *qnlist, double s1, double s2, double J,
                                  int L, int use_jj_mode, int nmax, double rmin, double rmax) {
    // 先构建 SL 子列表（对应 SL = {}）
    int sl_count = 0;
    for (double S = fabs(s1 - s2); S <= s1 + s2 + 1e-10; S += 1.0) {
        if (fabs(S - L) <= J && J <= S + L) sl_count++;
    }

    // 再构建 JJ 子列表（对应 JJ = {}）
    int jj_count = 0;
    if (use_jj_mode) {
        for (double jl = fabs(L - s1); jl <= L + s1 + 1e-10; jl += 1.0) {
            if (fabs(jl - s2) <= J && J <= jl + s2) jj_count++;
        }
    } else {
        jj_count = sl_count;  // SL 模式下，jj_count = sl_count（单位矩阵）
    }

    // 分配临时数组（避免多次 realloc）
    argsBasis_t *sl_temp = malloc(sl_count * sizeof(argsBasis_t));
    argsBasis_t *jj_temp = malloc(jj_count * sizeof(argsBasis_t));
    int sl_idx = 0, jj_idx = 0;

    // 填充 SL
    for (double S = fabs(s1 - s2); S <= s1 + s2 + 1e-10; S += 1.0) {
        if (fabs(S - L) <= J && J <= S + L) {
            argsBasis_t q;
            q.s1 = s1; q.s2 = s2; q.sij = S; q.L = L; q.J = J;
            q.jl = 0.0;  // SL 模式下 jl 无意义，可设 0
            q.lrho = L; q.llam = 0;  // 假设 lrho=L, llam=0（meson Jacobi）
            // 径向（可选，若不需要可注释）
            for (int nr = 1; nr <= nmax; nr++) {
                for (int nl = 1; nl <= nmax; nl++) {
                    q.nrho = nr; q.nlam = nl;
                    q.nurho = getnu(rmin, rmax, nmax, nr);
                    q.nulam = getnu(rmin, rmax, nmax, nl);
                    // 推入（这里简化，实际可按需复制）
                    basis_list_push(qnlist, 1, -1, -1, 1.0,
                                    s1, s2, 0.0, s1, s2, 0.0,
                                    0.0,0.0,0.0,0.0,0.0,
                                    0, L, 0, L, S, 0.0, J,
                                    nr, nl, q.nurho, q.nulam);
                }
            }
            sl_temp[sl_idx++] = q;
        }
    }

    // 填充 JJ 或 SL 模式
    if (use_jj_mode) {
        for (double jl = fabs(L - s1); jl <= L + s1 + 1e-10; jl += 1.0) {
            if (fabs(jl - s2) <= J && J <= jl + s2) {
                argsBasis_t q = {0};  // 清零
                q.s1 = s1; q.s2 = s2; q.L = L; q.J = J; q.jl = jl;
                // sij 和其他字段根据需要设（通常 sij 不直接用）
                jj_temp[jj_idx++] = q;
            }
        }
    } else {
        // SL 模式：直接复制 SL 列表
        for (int i = 0; i < sl_count; i++) {
            jj_temp[i] = sl_temp[i];
            jj_temp[i].jl = 0.0;  // 或其他标记
        }
    }

    free(sl_temp);
    free(jj_temp);
}

// 主函数：支持 S-D 混合（多 L）
void build_basis_mixed(basis_list *qnlist, double s1, double s2, double J, int P,
                       int L_min, int L_max, int step, int use_jj_mode,
                       int nmax, double rmin, double rmax) {
    basis_list_init(qnlist);

    // 外层 L 循环（S-D: 0 和 2）
    for (int L = L_min; L <= L_max; L += step) {
        // parity 检查（可选，Mathematica 无，但物理上必要）
        if (parity(L) != P) continue;

        // 为当前 L 构建子基并推入总列表
        build_sub_basis_for_L(qnlist, s1, s2, J, L, use_jj_mode, nmax, rmin, rmax);
    }

    // 日志输出（验证）
    basis_list_logs(*qnlist);
}