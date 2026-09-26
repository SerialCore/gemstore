scdkme# 重子 Jacobi GEM（SCDK）

日期：2026-09-26  
范围：`SPECTRA` + `BARYON` + `GEM`  
对照代码：`src/model/cbaryon.c`、`src/math/sumckdk.c`、`src/math/scdkme.c`、`src/basis/basis.c`、`src/basis/threebody.c`

重子谱走 **SCDK**（`recycle/sumckdk.h` + `inteCenV.h` + `vtype.h`）。介子仍用 1D GEM，矩阵元路径不共用。

Jacobi 映射矩阵与配平方（SCDK 的 `getT*` 用同一套 \(\alpha\beta\gamma\delta\)、\(b_{11}\)）见 **`doc/baryon-rr-gem-derivation.md`**。该文后半的「固谐加法 / RR 路径 1」是已废弃的实现，生产代码不再调用 `solidharm_central_me`。

```bash
pandoc doc/baryon-rr-gem-derivation.md -o doc/baryon-rr-gem-derivation.pdf \
  --pdf-engine=xelatex \
  -V CJKmainfont="Noto Serif CJK SC" \
  -V CJKmonofont="Noto Sans CJK SC" \
  --toc --number-sections
```

---

## 1. 现在能算什么

三张 Jacobi 连接图 \(c=1,2,3\)（对 \(12/31/23\)，即 123 / 312 / 231）写在**同一套** \(Hc=ENc\) 里。求 \(V_{13}\) 时 SCDK 把 bra/ket 的高斯从各自标架映到 \(\rho_{31}\) 再积。过完备 \(N\) 丢掉近零模后再解 \(H\)。

| 项目 | 现状 |
|---|---|
| 入口 | `--compute` → `SYSTEM_BARYON` → `compute_spectra_baryon` → `spectra_baryon_GEM` |
| 基 | `basis`。可分辨：三标架**独立**。全同对：该对 Pauli + 重排图 \(\lvert c_a\rangle+\eta\lvert c_b\rangle\) |
| \(L_{\max}\) | \(l_\rho+l_\lambda\le L_{\max}\) 且宇称对；\(L_{\max}=0\) S 波，\(1\) 含 \((1,0)\) 与 \((0,1)\) |
| 角向多项式 | `sumckdk_scdk_calc`，对 `qnlist_spfy` 缓存，径向点共用 |
| 运动学 | `scdkme`：`getTir` / `getT1p` / `getTpi` 等，配平方后对 \(t_j\) 求值 |
| 势 | recycle 径向核 `inteVogeG` / `inteVstring` / `inteVcont` / 张量 / SO；参数来自 JSON `model` |
| GISCREEN | \(\mu>0\) 时禁闭与 Thomas 走 smeared+screened 的 `InteCenV` |
| 夹心 | **每对**先 \(\beta V\beta\)，再三对相加（recycle `transP`）；介子是先对求和再夹 |
| 正交化 | `threebody_overlap_basis`：对角 \(N\)，丢掉 \(\lambda\le 10^{-8}\lambda_{\max}\)，行按 \(1/\sqrt{\lambda}\) 归一 |
| 输出 | 质量、三对 RMS（fm）、\(c^TNc\)；可选 `*.pot.dat`、`*.basis.dat`、`*.wfn.N.dat` |

**没有** \(S_3\) Young 投影（recycle 也没有）。sss 的 S 波 \(3/2^+\) 在 \(S_2(1,2)\) 下已经全对称；其它 \(S_3\) 类下次再做。

### 1.1 怎么跑

```bash
make
./gemstore --compute test/Spectra-Baryon/opal.json    # uds，1/2+，Lmax=0
./gemstore --compute test/Spectra-Baryon/opal_P.json  # uds，1/2-，Lmax=1
```

JSON：`system.f1,f2,f3`（1=n，2=s，3=c），`J`，`P`，`sym12`，`Lmax`；`basis.nmax,rmin,rmax`（fm）。

### 1.2 基：可分辨 vs 全同

通道约定与 `jacobi.h` 一致：

- \(c=1\)：\(\rho=r_1-r_2\)，spectator 3（123）
- \(c=2\)：\(\rho=r_3-r_1\)，spectator 2（312；与 132 同一对，差 \(\rho\to-\rho\)）
- \(c=3\)：\(\rho=r_2-r_3\)，spectator 1（231）

无质量权重：

\[
\rho_c=r_i-r_j,\qquad
\lambda_c=\frac{m_i r_i+m_j r_j}{m_i+m_j}-r_k.
\]

线性映射 \(\rho_{\mathrm{from}}=\alpha\rho_{\mathrm{to}}+\beta\lambda_{\mathrm{to}}\)，\(|\alpha\delta-\beta\gamma|=1\)。SCDK 的 `getT*` 用这套矩阵把任意 \(c_a,c_b\) 映到势所在的 pair。

| JSON | 构造 |
|---|---|
| **csn**（三味都不同，如 `f1,f2,f3=1,2,3`） | 三标架各成一套独立基，**不用** `sym12` |
| **ssn**（两个相同放在 1、2） | \(c=1\)：该对 Pauli；\(c=2\) 与 \(c=3\) 收成 \(\lvert 2\rangle+\eta\lvert 3\rangle\) |
| **sss**（三个相同） | 与 ssn 相同：只显式做 \(S_2(1,2)\)，不是 \(S_3\) |

Pauli / 叠加系数（recycle）：

\[
\eta=f_{12}\,(-1)^{1+s_{ij}+l_\rho},\qquad
\text{留下 } \eta=+1.
\]

`sym12=-1` 且 \(l_\rho=0\) 时留下 \(s_{ij}=0\)；`sym12=+1` 留下 \(s_{ij}=1\)。色单态 \(\varepsilon_{abc}\) 要求空间–自旋对称时，全同对一般用 `sym12=-1`；sss 的 \(3/2^+\) S 波应对 `sym12=+1`（\(s_{ij}=1\)）。

\(2\leftrightarrow 3\)、\(3\leftrightarrow 1\) 全同时，把对应的两张重排图按同样 \(\eta\) 组合。实现：`threebody_pair_identical`、`threebody_exchange_eta`，组装在 `baryon_basis_build`。

不要把三个标架的谱直和相加。交叉块

\[
N_{12}=\langle\phi^{(c=1)}\lvert\phi^{(c=2)}\rangle\neq 0
\]

由 SCDK 的 `inteNfi`（\(V=1\)，映到某一对标架）给出。

### 1.3 SCDK 矩阵元

1. `qnlist_spfy`：角动量通道（无径向 \(n\)）。`mlsj` 把 \(\lvert(s_{ij}L)j_l,s_3;J,M=J\rangle\) 拆成投影。
2. 21 组多项式（7 类算符 × 3 对）：`sumckdk_scdk_vtype` → `vcent/vcont/vtens/vsoii/jj/ji/ij`。与 \(\nu\) 无关。
3. `qnlist_full`：每个 spfy 态 \(\times(n_\rho,n_\lambda=1..n_{\max})\)，\(\nu=\mathrm{getnu}\)（`GEMSTORE_FM`）。
4. `getmfi`：用 `map1/map2` 取多项式，`getT*` 填 \(t_j\)（含 \(b_{11}\)、spectator 矩、四个 GEM 归一），`inteVcenPartA` 收缩 × 径向核。

动能是三夸克相对论单粒子能量之和（`tpi_cent` + `inteTi`）。动量夹心 \(\beta,\delta\) 用 `t1p_cent`，**每对** `matrix_sandwich` 后再相加。

### 1.4 已核对的数字

`GISTRING_MESON`、uds（`1,2,3`）、\(J=1/2\)、三标架独立、`nmax=3`、`rmin=0.1`、`rmax=3.0` fm：

| JSON | \(J^P\) | \(L_{\max}\) | 角动量通道 | 全基 | 基态 |
|---|---|---|---|---|---|
| `opal.json` | \(1/2^+\) | 0 | 6（每标架 2 个 \(s_{ij}\)） | 54 | **2.573 GeV**，\(r_{12}\approx 0.52\) fm，\(c^TNc=1\) |
| `opal_P.json` | \(1/2^-\) | 1 | 18 | 162 | **2.918 GeV** |

P 波高于 S 波。单通道 \(c=1\) 时 S 波基态约 2.59 GeV；三图独立后略降，是过完备展开，不是直和。`nmax/rmin/rmax` 改了数字会变。

recycle 实际跑的是**只建 \(c=1\)**（\(c=2,3\) 写了又注释）。要对齐 recycle 的旧数字，需关三图；当前默认按 Hiyama 三连接图进同一套 \(H\)。

---

## 2. 文件

| 文件 | 角色 |
|---|---|
| `include/gemstore/basis/basis.h`、`src/basis/basis.c` | 三体基（recycle `basis.h`） |
| `include/gemstore/basis/threebody.h`、`src/basis/threebody.c` | 全同对、SCDK 表、\(N\) 截断（分子态可复用） |
| `include/gemstore/basis/jacobi.h`、`src/basis/jacobi.c` | 通道质量、\(r\) 映射；`raynal_revai` 留给 \(S_3\) |
| `include/gemstore/math/sumckdk.h`、`src/math/sumckdk.c` | \(Y_{\ell m}\to\sum c(\mathbf{d}\cdot\hat n)^\ell\) |
| `include/gemstore/math/scdkme.h`、`src/math/scdkme.c` | `vtype` + `getT*` + 径向积分 |
| `include/gemstore/model/cbaryon.h`、`src/model/cbaryon.c` | 夸克 GI 装配：基、21 算符、夹心、本征 |
| `src/model/compute.c`、`src/print.c` | 分发、打印、`*.state.json`、势/波函数 |
| `include/gemstore/math/matrix.h` | `matrix_init` **置零**；`copy` / `symmetrize` / `expect` / `sandwich` |
| `test/Spectra-Baryon/opal.json`、`opal_P.json` | S / P 输入 |

`solidharm.c` 仍在树里，重子谱不再调用。不要提交根目录的 `*.state.json`。

---

## 3. 装配要点

```
baryon_basis_build          三图；全同时 Pauli + η 组合
basis_list_push_full        nρ,nλ × getnu
baryon_mlsj_jl              M=J 投影
threebody_scdk_table_alloc  21 × (spfy)⁴
thread  calc_scdk_mt        角向多项式
thread  getmfi              N,T,V,p,⟨r²⟩
matrix_symmetrize
threebody_overlap_basis     丢掉 N 的核
每对 uMu 与 βVβ 再求和
eigen_standard(H)
RMS = sqrt(⟨r_ij²⟩) / GEMSTORE_FM
```

`thread_load` 必须把**任务数组基址**传给每个 worker（recycle `mt_load`）；传 `&arg[i]` 会越界。

---

## 4. 尚未做

- **\(S_3\)**：sss 的全对称空间–自旋。系数要用 `raynal_revai` + 自旋 6j，不能把三图里同一套 \((l_\rho,l_\lambda,s_{ij})\) 直接相加。recycle 也没做。
- **GISCREEN 的 \(\mu\)** 已进禁闭 / Thomas；其它短程核仍是 recycle GI-string 形状。
- 分子三体：复用 `basis` / `sumckdk` / `getT*` / `threebody_*`；不要复用 21 个夸克 GI 算符和 `mlsj`。

---

## 5. 下次动手前

- 改矩阵元：`make`，跑 `opal.json` 与 `opal_P.json`。低态质量应稳定；\(N\)-正交后 \(c^TNc\approx 1\)；矢量符号不必 bit-for-bit。
- 改 `cbaryon.c` 或 `matrix_init` 时，用现有 meson JSON 做一次 GEM 回归。
- 根目录 `*.state.json` 是运行产物，不要进版本库。
