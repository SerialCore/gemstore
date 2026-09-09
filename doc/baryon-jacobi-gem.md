# 重子 Jacobi GEM：当前实现与未完成项

日期：2026-09-09  
范围：`SPECTRA` + `BARYON` + `GEM`  
对照代码：`src/model/cbaryon.c` 及下文文件表

本文记录第一版重子谱的交付边界、RR 路径怎么接到介子积分上、以及下次改三通道基 / 完整 GI 时不要踩的坑。不是用户手册。

公式、通道映射矩阵、配平方、固谐加法和矩阵元逐步推导见 **`doc/baryon-rr-gem-derivation.md`**（PDF：`doc/baryon-rr-gem-derivation.pdf`）。

```bash
pandoc doc/baryon-rr-gem-derivation.md -o doc/baryon-rr-gem-derivation.pdf \
  --pdf-engine=xelatex \
  -V CJKmainfont="Noto Serif CJK SC" \
  -V CJKmonofont="Noto Sans CJK SC" \
  --toc --number-sections
```

---

## 1. 现在能算什么

单 Jacobi 标架 \(c=1\) 上的中心势谱，含 P 波。

| 项目 | 状态 |
|---|---|
| 入口 | `gemstore --compute` → `SYSTEM_BARYON` → `compute_spectra_baryon` → `spectra_baryon_GEM` |
| 基 | 只建 \(c=1\)（\(\rho=r_1-r_2\)，\(\lambda\) 指向粒子 3） |
| \(L_{\max}\) | 0：S 波；1：P 波 \((l_\rho,l_\lambda)=(1,0)\) 与 \((0,1)\) |
| 动能 | \(T_\rho+T_\lambda\)：`GIVt`（对内两夸克）+ `GIVt_quark`（spectator），相对论 \(\sqrt{m^2+p^2}\) |
| 中心势 | \(V_{12}+V_{13}+V_{23}\) 的库仑 / 禁闭 / 接触径向都积 |
| 道内（`pair == c`） | 1D `GRnlr` / `GRnlp`，与介子相同；自旋–自旋、SO、张量走 `baryon_op_apply` |
| 道外（`pair != c`） | 中心势走 complete-the-square + 球谐加法定理（RR 路径 1）；\(\mathbf{s}_i\cdot\mathbf{s}_j\) 有自旋 recouple |
| 色因子 | \(C_{ij}=-2/3\)（介子是 \(-4/3\)） |
| 正交化 | 与介子相同：随机对称矩阵对 \(N\) 做广义本征，再变到正交基上拼 \(H\) |

**不是**完整 Godfrey–Isgur 重子，也**不是**三通道展开。

### 1.1 怎么跑

```bash
# 若 make 因时钟偏斜跳过链接：
rm -f gemstore obj/src/model/cbaryon.o obj/src/math/solidharm.o
make

./gemstore --compute test/Spectra-Baryon/opal.json    # S 波 1/2+
./gemstore --compute test/Spectra-Baryon/opal_P.json  # P 波 1/2-
```

运行一开始会做恒等检验（失败即 `exit(1)`）：

1. 道内 Coulomb：1D `GRnlr` vs RR `solidharm_central_me`（`pair == c`）
2. \(V=1\)：同一对高斯在 `pair=c` 与另一 pair 标架上的矩阵元（标架无关 = overlap）

### 1.2 已核对的数字

在 `GISTRING_MESON`、uds（`f1,f2,f3 = 1,2,3`）、\(J=1/2\)、`sym12=-1`、`nmax=4`、`rmin=0.2`、`rmax=2.0` 下：

| JSON | \(J^P\) | \(L_{\max}\) | 基态 | 状态数 |
|---|---|---|---|---|
| `opal.json` | \(1/2^+\) | 0 | **2.193 GeV** | 32，\(c^TNc=1\)，无零模、无负质量 |
| `opal_P.json` | \(1/2^-\) | 1 | **2.526 GeV** | 96，同上 |

P 波高于 S 波。若 JSON 里 `nmax/rmin/rmax` 已改，数字会变，以恒等检验 + 正质量为准。

独立（不跑全谱）还核对过 RR 空间元：

- \(V=1\)：pair 1/2/3 相对误差 \(\sim 10^{-16}\)，且与 \(M\) 无关
- \(V=1/r\)、\(\ell=1\)：与笛卡尔闭式一致（\(\alpha'^2\) 的 \(r\) 矩 + \(\beta^2\) 的 \(R\) 矩）

---

## 2. 两套 “\(c\)” 不要混

| 名字 | 代码 | 现在做什么 |
|---|---|---|
| 波函数的 Jacobi 标架 | `basis_qnum.c`，`baryon_basis_build` | **只建 \(c=1\)** |
| 势作用的夸克对 | 循环 `pair = 1,2,3` | **三对都积**。`pair==c` 走 1D；否则走 RR |

因此「P 波能算」= 在 \(c=1\) 的基上积了 \(V_{12}+V_{13}+V_{23}\) 的中心项，**不是**只算了 \(V_{12}\)。

通道约定（与 `jacobi.h` 一致）：

- \(c=1\)：对 \((1,2)\)，spectator 3
- \(c=2\)：对 \((3,1)\)，spectator 2
- \(c=3\)：对 \((2,3)\)，spectator 1

无质量权重坐标：

\[
\rho_c = r_i-r_j,\qquad
\lambda_c=\frac{m_i r_i+m_j r_j}{m_i+m_j}-r_k
\]

线性映射 \(\rho_{\mathrm{from}}=\alpha\rho_{\mathrm{to}}+\beta\lambda_{\mathrm{to}}\) 等，\(|\det(\alpha\beta\gamma\delta)|=1\)。

一套 Jacobi 坐标覆盖整个三体位形。\(L_{\max}\) 足够时，**只建 \(c=1\) 在数学上完备**。\(c=2,3\) 的高斯是另一套展开，与 \(c=1\) 线性相关，不是新的物理空间。

---

## 3. RR 路径（中心力，路径 1）

未采用 SCDK。道外中心力：

1. 把 bra/ket 的高斯从各自 `from` 标架映到势的 `pair` 标架
2. 配成平方：\(R'=R+\kappa r\)，\(\kappa=a_{rR}/(2a_{RR})\)，得到 \(e^{-b_{11}r^2-a_{RR}R'^2}\)
3. 移位后 \(\rho=(\alpha-\beta\kappa)r+\beta R'\)（`jacobi_gaussian_shift` 的 `al,be,ga,de`）
4. 固谐加法 \(\mathcal{Y}_{\ell m}(\alpha r+\beta R)\)，径向 \(\int r^{n}V(r)e^{-b_{11}r^2}dr\) 与 \(\int R^{n}e^{-a_{RR}R^2}dR\) 分离
5. 乘四个 `gem_pref`（`GRnlr` 去掉 \(r^\ell\) 的那一段，幂次改由 \(\mathcal{Y}\) 提供）

\(\ell=1\) 时 \(\mathcal{Y}_1(\alpha r+\beta R)=\alpha\mathcal{Y}_1(r)+\beta\mathcal{Y}_1(R)\) 精确成立。Bra 用 \(\mathcal{Y}_{\ell m}^*=(-1)^m\mathcal{Y}_{\ell,-m}\)，负磁量子数的相位用 `(m & 1)`，不要 `m % 2`（C 里负奇数 `% 2 == -1`）。

标量 \(\langle LM|V|LM\rangle\) 与 \(M\) 无关；实现里用拉伸态 \(M=L\)，已用 \(M=0,1\) 对过。

### 3.1 曾经的假态 / 零模（不要改回去）

| 现象 | 原因 | 现状 |
|---|---|---|
| 质量 \(\sim -16\)、\(-180\) GeV | 道外把 \(b_{11}\) 当 GEM \(\nu\) 送进 `integral_nlr_hamilton`，Hermite 权重与 `GRnlr` 的 \(\nu^{3/4}\) 双重计数 | 道外改 `integral_exp_rn` × `gem_pref` × spectator 高斯矩 |
| P 波 16 个精确零本征值，打开道外后变负质量 | `GIVt` 乘 `OCent`；同一通道若用介子 `operator_center_sl` recouple，1D 径向会把不同 \((l_\rho,l_\lambda,j_l)\) 搅在一起，正交化后出现 \(\dim=n_{\max}^2\) 的核；核上再叠道外吸引就塌缩 | 同一通道 `OCent = baryon_qn_match` |
| 势矩阵垃圾 | `matrix_init` 是 `malloc` 不置零，三对 `+=` | 每个 \((i,j)\) 在 `+=` 前先写成 0 |
| 道外 P 波曾被关掉 | \(\beta\neq 0\) 的 \(R\) 矩当时未验证 | 已验证，门已撤 |

`baryon_qn_match`：同一通道中心力 / 动能在 \((c,l_\rho,l_\lambda,L,s_{ij},j_l,J)\) 上 Kronecker。道外中心仍允许固定 \(L,s_{ij},j_l\) 下 \((1,0)\leftrightarrow(0,1)\) 混合，空间部分由 solidharm 出。

---

## 4. 改了哪些文件

| 文件 | 角色 |
|---|---|
| `include/gemstore/model/cbaryon.h`、`src/model/cbaryon.c` | 重子谱：基、\(N,T,V\)、RR 调用、恒等检验 |
| `include/gemstore/basis/jacobi.h`、`src/basis/jacobi.c` | 通道质量、\(r\) 映射、移位、`raynal_revai`（预备，谱里未用） |
| `include/gemstore/math/solidharm.h`、`src/math/solidharm.c` | 固谐加法、\(Y_{lm}\) 角积分、`solidharm_central_me` |
| `include/gemstore/math/integral.h`、`src/math/integral.c` | `integral_exp_rn` / `integral_exp_r2` |
| `src/model/gimodel.c`、`include/gemstore/model/gimodel.h` | `GIVt_quark`（单个 spectator） |
| `src/model/compute.c`、`src/entry.c` | `BARYON` 分发 |
| `test/Spectra-Baryon/opal.json`、`opal_P.json` | S / P 输入 |

`Makefile` 对 `src/math/*.c` 通配，`solidharm.c` 进树即可编。不要提交仓库根目录跑出来的 `opal.state.json`、`opal_P.state.json`。

`jacobi_angle`、`raynal_revai`、`jacobi_gaussian_overlap` 目前几乎无调用，是多通道 / RR 系数的预备，留着即可。

---

## 5. 关键代码位置

基只开 \(c=1\)：

```c
/* src/model/cbaryon.c  baryon_basis_build */
for (int c = 1; c <= 1; c++) {
```

全同对 Pauli 只作用在**当前通道的那一对**上（`pair_identical(f1,f2,f3,c)` + `sym12`）。uuc 把两个 u 放在粒子 1、2 时，\(c=1\) 的筛就是 uu 的 \((-1)^{s_{12}+l_\rho}\)。

矩阵元路由：

```c
/* me_pair_spatial */
if (qa->c == pair && qb->c == pair) { /* 1D GRnlr / GRnlp */ }
if (momentum) return 0.0;             /* 道外 β,δ 未做 */
return me_pair_reduced_r(...);        /* 道外中心 RR */
```

`baryon_set_operators`：道外 `OLSi=OLSj=OTens=0`，`OSdS` 仍算。

重叠：

```c
Nfi = I_ρ I_λ * (baryon_qn_match ? 1 : 0);
```

因此 \(c_i\neq c_j\) 时 \(N_{ij}=0\)。**在改 `c <= 3` 之前必须先改这里**，否则就是错误的直和。

---

## 6. 尚未完成

### 6.1 三通道基（过完备 GEM，不是直和）

**不要**把三个标架的谱 \(H_1\oplus H_2\oplus H_3\) 加起来。正确的是同一套广义本征问题里的分块：

\[
N=\begin{pmatrix}N_{11}&N_{12}&N_{13}\\N_{21}&N_{22}&N_{23}\\N_{31}&N_{32}&N_{33}\end{pmatrix},\quad
Hc=ENc,\quad N_{12}=\langle\phi^{(c=1)}|\phi^{(c=2)}\rangle\neq 0.
\]

只把 `for (c=1;c<=1)` 改成 `c<=3` **会算错**：

- `baryon_qn_match` 把 \(N_{12}\) 强制为 0
- \(T_{ij}\) 把不同标架的 \(\rho_i,\rho_j\) 当同一矢量做 1D 积分，再被 `OCent` 乘成 0
- 势的 RR 路径倒是按 `from_a,from_b → pair` 映射的，会与 \(N,T\) 不一致

若要做，至少：

1. **基** `baryon_basis_build`：打开 \(c=1,2,3\)（或只加少量重排道高斯）。每个 \(c\) 的 Pauli 用该对是否全同。uuc 的 \(c=2,3\) 在 \(\rho\) 上不是 uu，两个 u 的反对称要写成 \(c=2\) 与 \(c=3\) 的组合，或干脆只用 \(c=1\)。
2. **\(N\)**：同通道维持 \(I_\rho I_\lambda\)；异通道用已有 `jacobi_gaussian_shift` + `solidharm_central_me(V=1)`（或给 `jacobi_gaussian_overlap` 补 P 波）。
3. **\(T\)**：交叉项必须把 \(p_\rho,p_\lambda\) 变到同一标架，或对动能做与中心势相同的 RR。禁止对不同 \(c\) 的 `rho[i],rho[j]` 直接 `integral_nlp_hamilton`。
4. **势 / `OCent`**：道间不要再 `qn_match`。自旋 recouple 从「同一 \(c\)」推广到 \((c_a,c_b,\mathrm{pair})\)。
5. **数值**：过完备后 \(N\) 近奇异。现有 random-\(N\) 会大量 `singular vector`。需要丢掉小本征值，或 \(c=1\) 用全套径向、\(c=2,3\) 只用少数宽高斯。
6. **回归**：单通道 \(c=1\) 与「\(c=1\) + 少量 \(c=2\)」的低态质量在截断后应一致。差一截多半是 \(N_{12}\) 或 \(T_{12}\) 仍按直和写。

什么时候才需要：uud/uuu 的交换对称，或重排 / 散射。uds、以及把全同对放在粒子 1、2 上的 uuc，**单通道 \(c=1\) 即可**，不是功能缺失。若只想「主通道选 uu 还是 uc」，换 JSON 里谁当 `f1,f2,f3` 更便宜，仍是单通道。

### 6.2 完整 GI 重子

当前相对完整 GI 缺的是算符，不是三体运动学。

| 缺项 | 代码 | 后果 |
|---|---|---|
| 道外 \(\mathbf{L}_{ij}\cdot\mathbf{S}_{ij}\)、张量 | `baryon_set_operators` 里 `pair!=c` 时 `OLSi=OLSj=OTens=0` | P 波精细结构少两对；不能拿劈裂和实验 / recycle 逐态比 |
| 道外动量夹心 \(\beta(p),\delta(p)\) | `me_pair_spatial(..., momentum=1)` 道外 `return 0` | \(V_{13},V_{23}\) 的库仑 / 接触是**不 smear** 的 \(V(r)\)；\(V_{12}\) 仍与介子相同 |
| 同一通道 SO/张量的 recouple | `baryon_op_apply` + 介子 `operator_*_sl` | 中心项已不用这条；自旋力没有单独金标 |
| 三对先求和再 \(\beta V\beta\) | 沿用介子 `matrix_productT` | 重子会混进「这对的 \(\beta\) 夹那对的 \(V\)」 |

道外 SO/张量不能靠现在的**标量** solid-harmonic 直接推广，要在 \(\mathbf{r}_{ij}\) 上造向量 / 二阶张量。道外 \(\beta(p)\) 还要把动量搬到另一套 Jacobi，工作量与 SCDK 同级，但语言仍可留在 GEM+RR，**不必为此换 SCDK**。

建议顺序：需要精细结构时先做道外 \(\mathbf{L}\cdot\mathbf{S}\)/张量（接触 \(\mathbf{s}_i\cdot\mathbf{s}_j\) 三对已有）；要严格 GI 短程再做道外 \(\beta,\delta\)；要全同夸克或散射再做 6.1。

---

## 7. 下次动手前

- 改矩阵元后：`rm gemstore` 再 `make`，然后跑 `opal.json` 与 `opal_P.json`。恒等失败会直接退出。
- 不要在没做交叉 \(N,T\) 时把基循环改成 `c<=3`。
- 根目录 `*.state.json` 是运行产物，不要进版本库。
- 正交化用随机矩阵：本征值应稳定，矢量符号 / 顺序不必 bit-for-bit。
- 改 `cbaryon.c` 里与介子共用的拼 \(H\) 逻辑时，用现有 meson JSON 做一次 GEM 回归。
