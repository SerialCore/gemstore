---
title: 重子三体：Jacobi 坐标、GEM 与 RR 路径推导
author: GEMSTORE
date: 2026-09-09
documentclass: ctexart
classoption:
  - a4paper
  - UTF8
geometry: margin=2.4cm
colorlinks: true
toc: true
numbersections: true
header-includes:
  - \usepackage{amsmath,amssymb,bm}
  - \allowdisplaybreaks
---

对照实现：`src/basis/jacobi.c`、`src/math/solidharm.c`、`src/model/cbaryon.c`、`src/basis/orbit.c`。  
工程边界与未完成项见 `doc/baryon-jacobi-gem.md`。

本文写的是**代码里实际用的约定**，不是通用教科书的质量加权 Jacobi。核对时以函数名为准。

---

## 1. 三体位形与 Jacobi 通道

三个夸克坐标 $\mathbf{r}_1,\mathbf{r}_2,\mathbf{r}_3$，质量 $m_1,m_2,m_3$。去掉质心后剩两个相对矢量。代码用**无质量权重**的 Jacobi 坐标（`jacobi.h`）：

通道 $c$ 指定一对 $(i,j)$ 和 spectator $k$：

| $c$ | 对 $(i,j)$ | spectator $k$ |
|---|---|---|
| 1 | $(1,2)$ | 3 |
| 2 | $(3,1)$ | 2 |
| 3 | $(2,3)$ | 1 |

$$
\vec{\rho}_c = \mathbf{r}_i - \mathbf{r}_j,
\qquad
\vec{\lambda}_c
= \frac{m_i\mathbf{r}_i+m_j\mathbf{r}_j}{m_i+m_j} - \mathbf{r}_k.
$$

`jacobi_pair_mass` 给出 $(m_i,m_j,m_k)$。约化质量（预备，动能并不用非相对论 $p^2/2\mu$）

$$
\mu_\rho = \frac{m_i m_j}{m_i+m_j},
\qquad
\mu_\lambda = \frac{(m_i+m_j)m_k}{m_i+m_j+m_k}.
$$

体积元：$\lvert\det\rvert=1$ 的线性映射下 $\mathrm{d}^3\rho\,\mathrm{d}^3\lambda$ 在通道之间不变。

当前基只在 $c=1$ 上展开。一套坐标已经覆盖整个相对位形；$c=2,3$ 的高斯是另一套过完备展开，不是直和空间。

---

## 2. 通道之间的线性映射

把基通道 `from` 的 $(\vec{\rho},\vec{\lambda})$ 写成算符通道 `to` 的线性组合（`jacobi_r_map`）：

$$
\vec{\rho}_{\mathrm{from}}
= \alpha\,\vec{\rho}_{\mathrm{to}} + \beta\,\vec{\lambda}_{\mathrm{to}},
\qquad
\vec{\lambda}_{\mathrm{from}}
= \gamma\,\vec{\rho}_{\mathrm{to}} + \delta\,\vec{\lambda}_{\mathrm{to}}.
$$

下文把算符通道的相对坐标写成 $\mathbf{r}\equiv\vec{\rho}_{\mathrm{to}}$（势依赖的那一对）、spectator 写成 $\mathbf{R}\equiv\vec{\lambda}_{\mathrm{to}}$。

### 2.1 矩阵（行 = from，列 = to，下标 0 对应 $c=1$）

$\alpha$：

$$
\begin{pmatrix}
1 & -\dfrac{m_3}{m_1+m_3} & -\dfrac{m_3}{m_2+m_3} \\[8pt]
-\dfrac{m_2}{m_1+m_2} & 1 & -\dfrac{m_2}{m_2+m_3} \\[8pt]
-\dfrac{m_1}{m_1+m_2} & -\dfrac{m_1}{m_1+m_3} & 1
\end{pmatrix}
$$

$\beta$：

$$
\begin{pmatrix}
0 & 1 & -1 \\
-1 & 0 & 1 \\
1 & -1 & 0
\end{pmatrix}
$$

$\gamma$：

$$
\begin{pmatrix}
0 &
-1+\dfrac{m_2 m_3}{(m_1+m_2)(m_1+m_3)} &
1-\dfrac{m_1 m_3}{(m_1+m_2)(m_2+m_3)} \\[8pt]
1-\dfrac{m_2 m_3}{(m_1+m_2)(m_1+m_3)} &
0 &
-1+\dfrac{m_1 m_2}{(m_1+m_3)(m_2+m_3)} \\[8pt]
-1+\dfrac{m_1 m_3}{(m_1+m_2)(m_2+m_3)} &
1-\dfrac{m_1 m_2}{(m_1+m_3)(m_2+m_3)} &
0
\end{pmatrix}
$$

$\delta$：

$$
\begin{pmatrix}
1 & -\dfrac{m_2}{m_1+m_2} & -\dfrac{m_1}{m_1+m_2} \\[8pt]
-\dfrac{m_3}{m_1+m_3} & 1 & -\dfrac{m_1}{m_1+m_3} \\[8pt]
-\dfrac{m_3}{m_2+m_3} & -\dfrac{m_2}{m_2+m_3} & 1
\end{pmatrix}
$$

恒有 $\lvert\alpha\delta-\beta\gamma\rvert=1$。逆映射（$\det=+1$ 时）

$$
\mathbf{r} = \delta\,\vec{\rho}-\beta\,\vec{\lambda},
\qquad
\mathbf{R} = -\gamma\,\vec{\rho}+\alpha\,\vec{\lambda}.
$$

### 2.2 抽查：$c=1\to c=2$

$\vec{\rho}_1=\mathbf{r}_1-\mathbf{r}_2$，$\vec{\rho}_2=\mathbf{r}_3-\mathbf{r}_1$，
$\vec{\lambda}_2=\dfrac{m_3\mathbf{r}_3+m_1\mathbf{r}_1}{m_1+m_3}-\mathbf{r}_2$。

令 $\vec{\rho}_1=\alpha\vec{\rho}_2+\beta\vec{\lambda}_2$。$\mathbf{r}_2$ 的系数给出 $\beta=1$；$\mathbf{r}_1$ 的系数给出 $\alpha=-m_3/(m_1+m_3)$。与上表 $a_{12},b_{12}$ 一致。

等质量 $m_1=m_2=m_3$：

$$
c=1\to 2:\quad
\alpha=-\tfrac12,\; \beta=1,\; \gamma=-\tfrac34,\; \delta=-\tfrac12.
$$

### 2.3 质量加权角（未用于中心势路径）

`jacobi_angle` 定义

$$
\phi = \mathrm{atan2}\!\left(
\sqrt{\mu_\rho^{\mathrm{from}}}\,\beta/\sqrt{\mu_\lambda^{\mathrm{to}}},\;
\sqrt{\mu_\rho^{\mathrm{from}}}\,\alpha/\sqrt{\mu_\rho^{\mathrm{to}}}
\right),
$$

供 `raynal_revai` 用。当前中心力不走 RR 系数，而走第 7 节的固谐加法。

---

## 3. GEM 径向约定（与介子相同）

坐标空间高斯（`GRnlr`）**不写指数**，指数在半程 Hermite 权重里：

$$
R_{nl}(r;\nu)
= 2^{l/2+5/4}\,\nu^{l/2+3/4}\,
\frac{r^{l}}{\sqrt{\Gamma(l+3/2)}}
\qquad
\bigl(\times\; e^{-\nu r^2}\bigr).
$$

代码：`pow(2, l/2+1.25)`、`pow(nu, l/2+0.75)`。与 $n$ 无关；不同径向点只改 $\nu$。

完整轨道（含角向）

$$
\phi_{\nu\ell m}(\mathbf{x})
= N_{\ell}(\nu)\,
\mathcal{Y}_{\ell m}(\mathbf{x})\,
e^{-\nu x^2},
\qquad
N_{\ell}(\nu)
= 2^{l/2+5/4}\nu^{l/2+3/4}/\sqrt{\Gamma(l+3/2)}.
$$

$N_\ell$ 即 `gem_pref(l, nu)`：把 $r^l$ 从 `GRnlr` 里剥掉，改由固谐 $\mathcal{Y}_{\ell m}(\mathbf{x})=x^\ell Y_{\ell m}(\hat{\mathbf{x}})$ 提供。

动量空间（`GRnlp`）对应 $e^{-p^2/(4\nu)}$：

$$
\tilde R_{nl}(p;\nu)
= 2^{-l/2-1/4}\nu^{-l/2-3/4}
\frac{p^{l}}{\sqrt{\Gamma(l+3/2)}}
(-i)^{l}
\qquad
\bigl(\times\; e^{-p^2/(4\nu)}\bigr).
$$

Hermite 节点缩放：

- 坐标：$\mathrm{node\_factor}=1/\sqrt{\nu_a+\nu_b}$
- 动量：$\mathrm{node\_factor}=\sqrt{4\nu_a\nu_b/(\nu_a+\nu_b)}$

积分是 $\int_0^\infty r^2\,\mathrm{d}r$ 的半程 Hermite（`OHP=50`）。

径向宽度（`getnu`）：`rmin,rmax` 以 fm 输入，乘 $1/\hbar c = 5.06773\,\mathrm{GeV}^{-1}/\mathrm{fm}$，

$$
\nu_n
= \left(\frac{r_{\max}}{r_{\min}}\right)^{(2-2n)/(n_{\max}-1)}
\big/ (r_{\min}/\hbar c)^{2}.
$$

重子每个角动量通道是 $n_\rho,n_\lambda=1,\ldots,n_{\max}$ 的直积，维数 $n_{\max}^2$。

---

## 4. 角动量耦合与基

单通道基（`baryon_basis_build`）

$$
\bigl\lvert
(l_\rho l_\lambda)L,\; s_{ij};\; j_l,\; s_3=\tfrac12;\; J
\bigr\rangle
\otimes
\lvert n_\rho n_\lambda;\nu_\rho\nu_\lambda\bigr\rangle,
$$

空间部分

$$
\sum_{m_\rho m_\lambda}
\langle l_\rho m_\rho,\, l_\lambda m_\lambda \mid L M\rangle
\,\phi_{\nu_\rho l_\rho m_\rho}(\vec{\rho})
\,\phi_{\nu_\lambda l_\lambda m_\lambda}(\vec{\lambda}).
$$

筛选：

- 宇称 $P=(-1)^{l_\rho+l_\lambda}$，故 $L_{\max}=0$ 只有 S 波，$L_{\max}=1$ 且 $P=-1$ 只有 $(l_\rho,l_\lambda)=(1,0)$ 或 $(0,1)$。
- $L\in[\lvert l_\rho-l_\lambda\rvert,\,l_\rho+l_\lambda]$，$j_l\in[\lvert L-s_{ij}\rvert,\,L+s_{ij}]$，再与 $s_3$ 耦到 $J$。
- 若该通道的一对全同：$\mathrm{sym}_{12}\cdot(-1)^{s_{ij}+l_\rho}=+1$ 才保留（`pair_identical` + `f12`）。

`baryon_qn_match` 是同一标架上的 Kronecker：

$$
\delta_{c,c'}\delta_{l_\rho l_\rho'}\delta_{l_\lambda l_\lambda'}
\delta_{LL'}\delta_{s_{ij}s_{ij}'}\delta_{j_l j_l'}\delta_{JJ'}.
$$

---

## 5. 重叠 $N$ 与动能 $T$

同一 Jacobi 标架（当前所有基）：

$$
N_{ab}
= \langle R_{n_\rho^a l_\rho} R_{n_\rho^b l_\rho}\rangle_\rho
\,\langle R_{n_\lambda^a l_\lambda} R_{n_\lambda^b l_\lambda}\rangle_\lambda
\times
[\texttt{baryon\_qn\_match}].
$$

1D 重叠**没有**角向 $\delta_{ll'}$（`GRnlr` 对 $l\neq l'$ 径向积分一般非零），角向正交全靠 `baryon_qn_match`。异通道若将来打开，不能再用这条 1D 公式。

动能是 GI 相对论单夸克能量之和。在通道 $c$ 的坐标里：

$$
T
= \sqrt{m_i^2+p_\rho^2}+\sqrt{m_j^2+p_\rho^2}
+\sqrt{m_k^2+p_\lambda^2}.
$$

代码：`GIVt` 用 $m_i,m_j$ 积在 $\rho$ 上，乘 $\lambda$ 重叠；然后 `mi=mj=m_k`，`GIVt_quark` 积在 $\lambda$ 上，乘 $\rho$ 重叠。两者都乘 `OCent`。因此同一通道的 `OCent` 必须是第 4 节的全套 Kronecker；若只用介子 `operator_center_sl`（只认 $s_{ij},l_\rho$），1D 径向会把不同 $(l_\rho,l_\lambda,j_l)$ 搅在一起，$H$ 出现 $n_{\max}^2$ 维核。

色因子：重子对 $C_{ij}=-2/3$（介子 $-4/3$），经 `args_dynmc->Cij` 进入 GI 势。

---

## 6. 道内中心势（`pair == c`）

势只依赖 $\rho=r_{ij}$，$\lambda$ 上是重叠。与介子相同：

$$
\langle a\lvert V(\rho)\rvert b\rangle
=
\Bigl(\int_0^\infty r^2\mathrm{d}r\;
R_a(r)V(r)R_b(r)\Bigr)
\times
\langle\lambda_a\lvert\lambda_b\rangle.
$$

动量空间夹心 $\beta(p),\delta(p)$ 同样只在 $\rho$ 上走 `GRnlp`。自旋–轨道、张量、自旋–自旋在 $\rho$ 上用介子算符，经

$$
\bigl\lvert(l_\rho l_\lambda)L,s_{ij};j_l\bigr\rangle
\;\longrightarrow\;
\bigl\lvert(s_{ij} l_\rho)j_\rho,\,l_\lambda;j_l\bigr\rangle
$$

接到 `operator_*_sl`（`baryon_op_apply`）。实现里

$$
\langle\mathrm{rec}\rangle
= (-1)^{s_{ij}+l_\lambda+j_l+L}
\sqrt{(2L+1)(2j_\rho+1)}
\begin{Bmatrix}
l_\rho & l_\lambda & L \\
j_l & s_{ij} & j_\rho
\end{Bmatrix}.
$$

`sixJ_symbol(j1,j2,j12,j3,j,j23)` 即 $\{j_1 j_2 j_{12};\, j_3\, j\, j_{23}\}$。中心项已不再走这条 recouple。核对 SO/张量时应用标准 Racah 公式对照 6j 参数顺序。

---

## 7. 道外中心势：配平方 + 固谐加法（RR 路径 1）

算符 $V(\lvert\mathbf{r}\rvert)$，$\mathbf{r}=\vec{\rho}_{\mathrm{pair}}$。bra/ket 的高斯写在 `from_a`、`from_b`（当前都是 $c=1$）。

### 7.1 二次型

$$
\nu_\rho^a\rho_a^2+\nu_\lambda^a\lambda_a^2
+\nu_\rho^b\rho_b^2+\nu_\lambda^b\lambda_b^2
= a_{rr}\,r^2 + a_{rR}\,\mathbf{r}\cdot\mathbf{R} + a_{RR}\,R^2,
$$

$$
\begin{aligned}
a_{rr}
&= \nu_\rho^a\alpha_a^2+\nu_\lambda^a\gamma_a^2
+\nu_\rho^b\alpha_b^2+\nu_\lambda^b\gamma_b^2,\\
a_{rR}
&= 2\bigl(
\nu_\rho^a\alpha_a\beta_a+\nu_\lambda^a\gamma_a\delta_a
+\nu_\rho^b\alpha_b\beta_b+\nu_\lambda^b\gamma_b\delta_b
\bigr),\\
a_{RR}
&= \nu_\rho^a\beta_a^2+\nu_\lambda^a\delta_a^2
+\nu_\rho^b\beta_b^2+\nu_\lambda^b\delta_b^2.
\end{aligned}
$$

### 7.2 移位

$$
\kappa = \frac{a_{rR}}{2 a_{RR}},
\qquad
\mathbf{R}'=\mathbf{R}+\kappa\mathbf{r}.
$$

则

$$
a_{rr}r^2+a_{rR}\mathbf{r}\cdot\mathbf{R}+a_{RR}R^2
= b_{11} r^2 + a_{RR} (R')^{2},
\qquad
b_{11}=a_{rr}-\frac{a_{rR}^2}{4 a_{RR}}.
$$

$\mathrm{d}^3R=\mathrm{d}^3R'$。坐标本身

$$
\vec{\rho}
= (\alpha-\beta\kappa)\,\mathbf{r} + \beta\,\mathbf{R}'
\equiv \alpha'\mathbf{r}+\beta\mathbf{R}',
\qquad
\vec{\lambda}
= (\gamma-\delta\kappa)\,\mathbf{r} + \delta\,\mathbf{R}'.
$$

即 `jacobi_shift_t` 的 `al, be, ga, de`（`be` 仍是原 $\beta$）。要求 $b_{11}>0$、$a_{RR}>0$。

等质量、$\nu_\rho=0.5$、$\nu_\lambda=0.4$、$c=1\to\mathrm{pair}=2$ 时代码给出

$$
b_{11}=\tfrac23,\quad a_{RR}=1.2,\quad
\alpha'=-\tfrac13,\;\beta=1,\;
\gamma'=-\tfrac56,\;\delta=-\tfrac12.
$$

### 7.3 固谐加法

$$
\mathcal{Y}_{\ell m}(\alpha\mathbf{r}+\beta\mathbf{R})
= \sqrt{4\pi}
\sum_{k=0}^{\ell}
\alpha^{k}\beta^{\ell-k}
\sqrt{\frac{(2\ell+1)!}{(2k+1)!\,(2\ell-2k+1)!}}
\sum_{m_1}
\langle k m_1,\,\ell-k,\,m-m_1\mid \ell m\rangle
\,Y_{k m_1}(\hat{\mathbf{r}})
\,Y_{\ell-k,\,m-m_1}(\hat{\mathbf{R}}).
$$

代码 `fact_int(2*ell+1)` 是 $(2\ell+1)!$。$\ell=0$：系数 $\sqrt{4\pi}\,Y_{00}Y_{00}=Y_{00}=\mathcal{Y}_{00}$。$\ell=1$ 精确为

$$
\mathcal{Y}_{1m}(\alpha\mathbf{r}+\beta\mathbf{R})
= \alpha\,\mathcal{Y}_{1m}(\mathbf{r})+\beta\,\mathcal{Y}_{1m}(\mathbf{R}).
$$

### 7.4 Bra 共轭

$$
\mathcal{Y}_{\ell m}^*=(-1)^m\mathcal{Y}_{\ell,-m}.
$$

对 bra 的 $m_\rho,m_\lambda$：CG 仍用 $(m_\rho,m_\lambda;LM)$，加法展开 $\mathcal{Y}_{\ell,-m}$，再乘 $(-1)^{m_\rho+m_\lambda}$。相位用 `(m & 1)`，**不要** `m % 2`（C 里负奇数 `% 2 == -1`）。

### 7.5 6D 积分

移位后指数分离。四个固谐（bra 的 $\rho,\lambda$ 与 ket 的 $\rho,\lambda$）各拆成 $r$ 与 $R'$ 上的 $Y_{lm}$。令

$$
n_r=\sum_{i=1}^{4} \ell_r^{(i)},
\qquad
n_R=\sum_{i=1}^{4} \ell_R^{(i)}.
$$

体积元 $\mathrm{d}^3x=x^2\mathrm{d}x\,\mathrm{d}\Omega$，故

$$
\begin{aligned}
I_R
&= \int_0^\infty R^{n_R+2} e^{-a_{RR} R^2}\,\mathrm{d}R
= \tfrac12 a_{RR}^{-(n_R+3)/2}\Gamma\bigl(\tfrac{n_R+3}{2}\bigr),\\
I_r
&= \int_0^\infty r^{n_r+2} V(r)\,e^{-b_{11} r^2}\,\mathrm{d}r
\quad\text{(Hermite, integral\_exp\_rn)}.
\end{aligned}
$$

$n_R=0$ 时 $I_R=\sqrt{\pi}/(4 a_{RR}^{3/2})$，与介子 spectator 约定一致。

角向：四个 $Y_{lm}$ 的 $\int\mathrm{d}\Omega$。乘积

$$
Y_{l_1 m_1}Y_{l_2 m_2}
=\sum_{l}
\sqrt{\frac{(2l_1+1)(2l_2+1)}{4\pi(2l+1)}}
\langle l_1 0, l_2 0\mid l 0\rangle
\langle l_1 m_1, l_2 m_2\mid l m\rangle
Y_{lm}.
$$

$\langle l_1 0, l_2 0\mid l 0\rangle=0$ 除非 $l_1+l_2+l$ 为偶。从 $1=\sqrt{4\pi}\,Y_{00}$ 起逐个乘，取出 $Y_{00}$ 系数再乘 $\sqrt{4\pi}$。$n=0$ 个 $Y_{lm}$ 时积分为 $4\pi$。

总空间元（尚无 `OCent`）

$$
\langle a\lvert V\rvert b\rangle_{\mathrm{sp}}
=
N_{l_\rho^a}(\nu_\rho^a)\,N_{l_\lambda^a}(\nu_\lambda^a)
\,N_{l_\rho^b}(\nu_\rho^b)\,N_{l_\lambda^b}(\nu_\lambda^b)
\times
\texttt{solidharm\_central\_me}(\ldots;\,b_{11},a_{RR},V).
$$

标量 $\langle LM\lvert V\rvert LM\rangle$ 与 $M$ 无关；实现取拉伸态 $M=L$。

**禁止**把 $b_{11}$ 当作 GEM $\nu$ 送进 `integral_nlr_hamilton`：Hermite 权重已含 $e^{-b_{11}r^2}$，`GRnlr` 再带 $\nu^{l/2+3/4}$，会双重计数，出现 $\mathcal{O}(10^2)\,\mathrm{GeV}$ 假束缚。

### 7.6 $\ell=1$、$M=0$、$V=1/r$ 的笛卡尔闭式

$(l_\rho,l_\lambda)=(1,0)$，$\mathcal{Y}_{10}=\sqrt{3/4\pi}\,z$。交叉项 $z_r z_R$ 对中心 $V(\lvert\mathbf{r}\rvert)$ 积分为 0。剩

$$
\langle V\rangle
= N_1(\nu_\rho)^2 N_0(\nu_\lambda)^2
\Bigl[
{\alpha'}^2
\int_0^\infty r^{3}V_{\times}\,e^{-b_{11}r^2}\mathrm{d}r
\int_0^\infty R^{2}e^{-a_{RR}R^2}\mathrm{d}R
+
\beta^2
\int_0^\infty r\,V_{\times}\,e^{-b_{11}r^2}\mathrm{d}r
\int_0^\infty R^{4}e^{-a_{RR}R^2}\mathrm{d}R
\Bigr],
$$

其中 $V=1/r$ 时 $V_{\times}$ 已并进幂次（$\int r^4\cdot r^{-1}=\int r^3$ 等）。solidharm 与此式相对误差 $\sim 10^{-14}$。

$(0,1)$ 把 $\alpha',\beta$ 换成 $\gamma',\delta$。

### 7.7 $V=1$ 标架无关

$V=1$ 就是 overlap，与把哪一对当 $\mathbf{r}$ 无关。同一对高斯、`pair=1,2,3` 的 RR 元必须相等（相对误差 $\sim 10^{-16}$）。道内 Coulomb 的 1D `GRnlr` 必须等于 `pair==c` 的 RR。这两条在 `baryon_check_reduce_identity` 里，失败即 `exit(1)`。

---

## 8. 道外自旋与尚未接入的 GI 块

中心力的 `OCent`（道外）：$\delta_{JJ'}\delta_{LL'}\delta_{j_l j_l'}\delta_{s_{ij}s_{ij}'}$。**不**要求 $\delta_{l_\rho l_\rho'}$：固定 $L$ 时 $(1,0)$ 与 $(0,1)$ 可由 $V_{13}$ 混合，空间部分由 solidharm 出。

$\mathbf{s}_i\cdot\mathbf{s}_j$（道外）：先耦到总自旋 $S$，再耦到该 pair 的 pair-spin $s'$，

$$
\mathbf{s}_i\cdot\mathbf{s}_j
= \tfrac12\bigl(s'(s'+1)-\tfrac32\bigr).
$$

`recouple_jl_to_S`：

$$
(-1)^{s_{ij}+L+s_3+J}
\sqrt{(2j_l+1)(2S+1)}
\begin{Bmatrix}s_{ij}&L&j_l\\ s_3&J&S\end{Bmatrix}.
$$

`recouple_spin_pair`：同通道 $\delta_{s s'}$；异通道 $\sqrt{(2s+1)(2s'+1)}\,\{1/2\,1/2\,s;\,1/2\,S\,s'\}$，邻接通道 $(1\leftrightarrow 3,2\leftrightarrow 3)$ 再乘 $(-1)^{3/2+S}$。

未实现（显式 0）：

- 道外 $\mathbf{L}_{ij}\cdot\mathbf{S}_{ij}$、张量（标量加法定理不够）
- 道外 $\beta(p),\delta(p)$：`me_pair_spatial(..., momentum=1)` 在 `pair!=c` 时返回 0。故 $V_{13},V_{23}$ 的库仑/接触是不 smear 的 $V(r)$，$V_{12}$ 仍与介子相同。

---

## 9. 哈密顿装配（与介子同构）

在原始 GEM 基上构造 $N,T,V_{\mathrm{conf}},V_{\mathrm{coul}},\beta_{\mathrm{coul}},\ldots$。三对势对 `pair=1,2,3` **先求和**。

`matrix_init` 不置零，凡 `+=` 的矩阵必须在每个 $(i,j)$ 先写成 0。

然后：随机对称矩阵对 $N$ 做广义本征，把矢量按 $\sqrt{(v N v^T)_{kk}}$ 归一，得到 $N$-正交行 $v$。各矩阵变到该基：$tM = v\,M\,v^{T}$（`matrix_productT`）。

$$
H = T + V_{\mathrm{conf}}
+ \beta_{\mathrm{coul}} V_{\mathrm{coul}} \beta_{\mathrm{coul}}^{T}
+ \delta_{\mathrm{cont}} V_{\mathrm{cont}} \delta_{\mathrm{cont}}^{T}
+ \cdots
$$

再对 $H$ 做标准本征。三对先求和再夹心，会混进「这对的 $\beta$ 夹那对的 $V$」，是沿用介子写法的近似。

本征值应稳定；随机正交化不保证矢量符号 bit-for-bit。

---

## 10. 算法流程（对照代码）

```
baryon_basis_build          c=1 only; parity, Pauli, J
basis_list_push_full        n_rho, n_lam = 1..nmax, nu=getnu
baryon_check_reduce_identity
  same-channel Coulomb: 1D GRnlr vs solidharm (pair==c)
  V=1: pair=c vs another pair
for i,j:
  zero += matrices
  N = I_rho I_lam * qn_match
  set OCent = qn_match for pair==c; T_rho(GIVt)+T_lam(GIVt_quark)
  for pair=1..3:
    baryon_set_operators
    me_pair_spatial: pair==c -> 1D; else if r-space -> RR; else 0
random-N orthogonalize
PVP sandwich -> eigen_standard(H)
```

入口：`entry_compute` → `compute_spectra_baryon` → `spectra_baryon_GEM`。

---

## 11. 建议的手算 / 代码对照清单

1. 等质量 $1\to 2$：$\alpha,\beta,\gamma,\delta$ 是否为 $(-1/2,1,-3/4,-1/2)$。
2. 同一组 $\nu$ 的 $\kappa,\alpha',b_{11},a_{RR}$ 是否与 `jacobi_gaussian_shift` 一致。
3. $\ell=0$ 加法：四条 $\mathcal{Y}_{00}$ 与两个 $\int\mathrm{d}\Omega=1/(4\pi)$ 相乘得 1，径向只剩 $I_r[V]\,I_R$。
4. $\ell=1$ 加法：$\alpha\mathcal{Y}(r)+\beta\mathcal{Y}(R)$，无交叉项时笛卡尔 7.6 成立。
5. `gem_pref` 与 `GRnlr/r^l` 逐因子相同。
6. 道内 Coulomb：1D 与 RR 相对误差 $<10^{-8}$。
7. $V=1$：三个 pair 标架的 RR 元相同。
8. 对角 `OCent=1`、`T_{ii}>0`；P 波不应出现恰好 $n_{\max}^2$ 个零本征值。
9. S 波基态应低于 P 波（同模型、同径向格点）。

---

## 12. 明确不在本推导内的

- 三通道过完备基的交叉 $N_{12},T_{12}$（需要把 $p_\rho,p_\lambda$ 变到同一标架，或对 $V=1$ 与动能再用本节 RR）。
- 道外向量 / 张量算符（SCDK 或固谐的梯度形式）。
- 质量加权 Jacobi 与标准 RR 系数展开（`raynal_revai` 已实现但中心势未用）。
- CRG / SHO 重子。

## 13. 编译 PDF

显示公式用 `$$...$$`，行内用 `$...$`（pandoc 不会把 `\[...\]` 当数学）。矢量 Jacobi 坐标写成 `\vec{\rho}`，以便 `xelatex` 处理希腊字母。

```bash
pandoc doc/baryon-rr-gem-derivation.md -o doc/baryon-rr-gem-derivation.pdf \
  --pdf-engine=xelatex \
  -V CJKmainfont="Noto Serif CJK SC" \
  -V CJKmonofont="Noto Sans CJK SC" \
  --toc --number-sections
```
