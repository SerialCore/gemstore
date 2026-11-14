<div align="center">

![GEMSTORE Logo](doc/logo.png)

# GEMSTORE: Hadron Spectroscopy Simulation Tools

[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0-or-later)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Platform: Linux/Unix](https://img.shields.io/badge/Platform-Linux%2FUnix-lightgrey.svg)](https://en.wikipedia.org/wiki/Unix)

**A C program for hadron spectra, S-wave scattering lengths, meson parameter fits, and \(^3P_0\) decay widths.**

</div>

GEMSTORE solves quark-model meson and baryon spectra with the Gaussian expansion method (GEM), MIT bag masses from mesons through hexaquarks, those bags in symmetric nuclear matter, and two- and three-body hadronic molecules. A separate task fits built-in Godfrey–Isgur meson datasets with Minuit2. Two further tasks, outside the GEM eigenvalue problem, return an S-wave scattering length and a single-oscillator \(^3P_0\) width.

Strings in the JSON input are case-sensitive. Masses in spectrum files are in GeV. Lengths printed as fm are femtometres.

## Contents

1. [What this version computes](#what-this-version-computes)
2. [Requirements and build](#requirements-and-build)
3. [How to run](#how-to-run)
4. [Input file](#input-file)
5. [Parameter files](#parameter-files)
6. [Output](#output)
7. [Solvers](#solvers)
8. [Test suite](#test-suite)
9. [Source layout](#source-layout)
10. [Optional agent skill](#optional-agent-skill)
11. [License and citation](#license-and-citation)

## What this version computes

| Task | System | Model | Basis | Result |
|------|--------|-------|-------|--------|
| `SPECTRA` | `MESON` | `GISCREEN`, `GISTRING`, `NRSTRING`, `NRSCREEN` | `GEM` or `SHO` | Masses, RMS radii, eigenvectors |
| `SPECTRA` | `BARYON` | `GISCREEN`, `GISTRING` | `GEM` | Masses, three pair radii, eigenvectors |
| `SPECTRA` | `BAG` | `BAGMIT`, `QMCSTAD`, `QMCMOD` | `NONE` or `CMI` | Bag masses, radius, color–spin mixing; QMC also writes the nuclear-matter fields |
| `SPECTRA` | `MOLECUL2` | `OBESTAT`, `YUKAWA`, `EFTCONT`, `HALQCD` | `GEM` | Two-body mass, binding, RMS |
| `SPECTRA` | `MOLECUL3` | `HYBRID` of three pair models | `GEM` | Three-body mass, binding, pair RMS |
| `FITTING` | chosen by `fit.target` | pinned `GISCREEN` set | `GEM` or `SHO` | Minuit2 \(\chi^2\) and parameters |
| `SCATTER` | `MOLECUL2` | `OBESTAT`, `YUKAWA`, `EFTCONT`, `HALQCD` | none | S-wave scattering length and effective range |
| `3P0DECAY` | `MESON` | `QPC3P0` | none | Partial width of one parent into two daughters |

`NRSTRING` and `NRSCREEN` are non-relativistic Cornell Hamiltonians. They check the meson radial solver. Baryon spectra in this version are Godfrey–Isgur only. `SHO` is accepted for mesons and for meson fits. Bag and molecular spectra use the bases in the table.

The following problems are outside this version: non-relativistic baryon spectra, coupled channels, folded potentials, \(P\leftrightarrow P^*\) transitions, Coulomb scattering, a coupled \(\ell=0\) and \(\ell=2\) scattering equation, and a \(^3P_0\) width built from a GEM wave function. The decay task uses one harmonic-oscillator function for each meson.

## Requirements and build

- Linux or another Unix with GNU make
- GCC or Clang with C99 and a matching C++ compiler (`gcc` and `g++` in the Makefile)
- CMake, used only to build the bundled Minuit2
- POSIX threads and `libm`

Clone with the Minuit2 submodule and build:

```bash
git clone --recursive https://github.com/serialcore/gemstore.git
cd gemstore
make
./gemstore --help
./gemstore --version
```

If `lib/Minuit2` is empty, run `git submodule update --init --recursive` and then `make` again. `make install` copies the binary to `~/.local/bin`. `make uninstall` removes that copy. `make clean` deletes `obj/`, the binary, and the Minuit2 build directory.

The executable links `libm`, `libpthread`, and the static Minuit2 libraries. It does not link LAPACK or OpenBLAS. Python is not required to run the program. The helper `app/gemstore-assistant/scripts/generate_meson_inputs.py` needs Python 3.

## How to run

```text
gemstore [--compute FILE] [--debug UNIT]
gemstore --help
gemstore --version
```

`--compute` reads one JSON file and dispatches on `task`. Output is written in the process working directory, under the `project` name. A relative `model.file` is opened from that directory, so run the program from the directory that holds the parameter files:

```bash
./gemstore --compute test/Spectra-Quarkonium/charmonium_1S0.json

cd test/Spectra-Molecule2
../../gemstore --compute dstar_xi_I0.json
```

`--debug` runs an internal check. The unit names are `su3_product`, `soc_operator`, `casimir_operator`, `color_wfn`, `spin_wfn`, `isospin_wfn`, `orbit_wfn`, `eigen_system`, and `eigen_system_complex`.

## Input file

Every file has a string `project` and a string `task`. The other objects depend on the task.

| Task | Required objects | Omit |
|------|------------------|------|
| `SPECTRA` | `model`, `system`, `basis`, `print` | `fit` |
| `FITTING` | `basis`, `fit` | `model`, `system`, `print` |
| `SCATTER` | `model`, `system` | `basis`, `print`, `fit` |
| `3P0DECAY` | `model`, `system` | `basis`, `print`, `fit` |

`print.pot` and `print.wfn` are the strings `"true"` or `"false"`. Any other string turns that dump off. `SPECTRA` still requires the `print` object.

### Meson spectrum

Flavors are `1` = \(n\) (up or down), `2` = \(s\), `3` = \(c\), `4` = \(b\). The quark masses are model parameters, not these indices.

```json
{
  "project": "charmonium_1S0",
  "task": "SPECTRA",
  "model": { "type": "GISCREEN", "param": "GISCREEN_CCBAR" },
  "system": { "type": "MESON", "f1": 3, "f2": 3, "S": 0, "L": 0, "J": 0 },
  "basis": { "type": "GEM", "nmax": 20, "rmax": 30.0, "rmin": 0.1 },
  "print": { "pot": "false", "wfn": "false" }
}
```

`S`, `L`, and `J` are the spin, orbital angular momentum, and total angular momentum of the \(q\bar q\) state. Worked inputs: `test/Spectra-Quarkonium/`, `test/Spectra-SHO/`, `test/Spectra-NRmodel/`, `test/Spectra-CustomParameter/`.

### Baryon spectrum

```json
"system": {
  "type": "BARYON",
  "f1": 1, "f2": 1, "f3": 3,
  "J": 0.5, "P": 1, "sym12": -1, "Lmax": 0, "jl": 0
}
```

`P` is \(+1\) or \(-1\). `sym12` is the sign under exchange of quarks 1 and 2. `Lmax` keeps channels with \(l_\rho+l_\lambda\le L_{\max}\) on the Jacobi charts. The parser requires `jl`. The basis is every \(j_l\) that couples to \(J^P\); the stored `jl` does not select a chart. Use `GISTRING_BARYON` for the smeared linear set. `GISTRING_MESON` is a different fit. Inputs: `test/Spectra-Baryon/`.

### MIT bag and quark-meson coupling

`nn`, `ns`, `nc`, and `nb` are the flavor counts. For `CMI`, `f1`–`f6` give the quark order and must match those counts. `NONE` does not need `f1`–`f6`.

```json
"system": { "type": "BAG", "nn": 2, "ns": 0, "nc": 2, "nb": 0,
            "f1": 3, "f2": 3, "f3": 1, "f4": 1 },
"basis": {
  "type": "CMI",
  "color": ["tetra11", "tetra88"],
  "spin": [ {"s12": 1, "s34": 1, "st": 2, "st3": 2} ]
}
```

Color names: `meson`, `baryon`, `tetra11`, `tetra33`, `tetra66`, `tetra88`, `penta688`, `penta388`, `penta11`, `hexa1010`, `hexa6868`, `hexa6838`, `hexa3868`, `hexa3838`, `hexa11`. A basis has 1–8 colors and 1–20 spin channels. Unused couplings in a spin object stay 0. Channels with different total spin `st` are diagonalized separately. For five and six quarks the program keeps the Pauli-allowed color\(\otimes\)spin products. Inputs: `test/Spectra-Bag/`.

`QMCSTAD` and `QMCMOD` use the same `BAG` system and add `rho_b` (fm\(^{-3}\)) and `temp` (GeV). Both are required and must be \(\ge 0\). `BAGMIT` does not read them. Inputs: `test/Spectra-QMC/`.

### Two-body molecule

```json
"model": { "type": "OBESTAT", "param": "OBESTAT_CUSTOM", "file": "param_OBESTAT_dstarxi.json" },
"system": {
  "type": "MOLECUL2",
  "type1": "VECTOR", "type2": "DIRAC",
  "m1": 2.010, "m2": 1.318,
  "t1": 0.5, "t2": 0.5, "I": 0,
  "S": 0.5, "L": 0, "J": 0.5
},
"basis": { "type": "GEM", "nmax": 20, "rmax": 30.0, "rmin": 0.1 }
```

`type1` and `type2` are `DIRAC` (spin-1/2 baryon), `VECTOR` (\(1^-\) meson), or `PSEUDOSCALAR` (elastic \(0^-\) or \(0^+\)). `m1` and `m2` are the physical masses in GeV. `t1`, `t2`, and the total isospin `I` are required. Two identical \(I=0\) baryons set `t1` and `t2` to `0`.

`L` is the lowest partial wave. The solver adds \(L+2\) and every total spin that couples to `J`. The JSON value of `S` does not choose that list. For example, \(J=1/2\) is \(^2S+^4D\) and \(J=3/2\) is \(^4S+^2D+^4D\).

The same system block is used with `YUKAWA`, `EFTCONT`, and `HALQCD`. Those models do not read the OBE factors `cc` and `eta`. `YUKAWA_MT13` has no parameter file. For the Malfliet–Tjon samples, the nucleon mass is `0.938942` GeV so that \(\hbar^2/m_N=41.47\) MeV·fm\(^2\).

`SCATTER` uses this system with `L` equal to 0 and with no `basis` or `print`. `S` and `J` select the partial wave.

### Three-body molecule

`model.type` is `HYBRID`. `model1`, `model2`, and `model3` are the pair interactions on legs 12, 13, and 23. Each pair is `OBESTAT`, `YUKAWA`, `EFTCONT`, or `HALQCD`, with the same `param` keys as the two-body model.

```json
"system": {
  "type": "MOLECUL3",
  "type1": "PSEUDOSCALAR", "type2": "PSEUDOSCALAR", "type3": "VECTOR",
  "m1": 1.867, "m2": 0.494, "m3": 0.892,
  "t1": 0.5, "t2": 0.5, "t3": 0.5,
  "J": 1.0, "P": -1, "I": 0.5, "Lmax": 2, "sym12": 0
}
```

Intrinsic parity is \(+1\) for a Dirac leg and \(-1\) for a vector or pseudoscalar. The total parity is \(P_1 P_2 P_3(-1)^{l_\rho+l_\lambda}\). `I` is the total isospin. Pair isospin is recoupled to `I`. `sym12` equal to `0` keeps three Jacobi charts. `sym12` equal to \(-1\) on three identical Dirac legs builds each component as the rearrangement sum of charts 1, 2, and 3, with the fermion phase \((-1)^{l_\rho+s_{ij}+t_{ij}}=-1\).

`nmax` is the number of Gaussian widths on each Jacobi axis. The dimension is \(n_{\mathrm{ang}}\times n_{\max}^2\). It is not a spectroscopic radial quantum number.

### Basis

| `basis.type` | Fields | Where |
|--------------|--------|--------|
| `GEM` | `nmax`, `rmin`, `rmax` (fm) | Meson, baryon, molecule, meson fit |
| `SHO` | `nmax`, `beta` (GeV) | Meson spectrum and meson fit |
| `CMI` | `color`, `spin` | Bag |
| `NONE` | none | Bag energy without color–magnetic mixing |

GEM functions are \(r^\ell\exp(-\nu r^2)\), with \(\nu\) placed from `rmin` to `rmax`. A single-`beta` SHO basis does not grow with excitation the way that GEM ladder does. Near-threshold molecular states with RMS close to `rmax` are box states of that basis.

### Fit

```json
{
  "project": "fit_csbar",
  "task": "FITTING",
  "fit": { "target": "GISCREEN_CSBAR" },
  "basis": { "type": "GEM", "nmax": 20, "rmax": 30.0, "rmin": 0.1 }
}
```

`fit.target` is one of `GISCREEN_MESON`, `GISCREEN_BBBAR`, `GISCREEN_BCBAR`, `GISCREEN_BSBAR`, `GISCREEN_CCBAR`, `GISCREEN_CSBAR`. The target selects the model and the experimental list. Starting values, steps, bounds, and which parameters are fixed live in `src/param/f*.c`.

### \(^3P_0\) decay

```json
{
  "project": "rho_pipi",
  "task": "3P0DECAY",
  "model": { "type": "QPC3P0", "param": "QPC3P0_CUSTOM", "file": "param_qpc3p0_rho.json" },
  "system": {
    "type": "MESON",
    "F": 1, "fl1": 0.7071067811865476, "fl2": -0.7071067811865476,
    "initial": { "f1": 1, "f2": 1, "n": 1, "S": 1, "L": 0, "J": 1, "mass": 0.77, "beta": 0.397 },
    "final": [
      { "f1": 1, "f2": 1, "n": 1, "S": 0, "L": 0, "J": 0, "mass": 0.138, "beta": 0.397 },
      { "f1": 1, "f2": 1, "n": 1, "S": 0, "L": 0, "J": 0, "mass": 0.138, "beta": 0.397 }
    ]
  }
}
```

Each meson is one SHO. `n` starts at 1 and is the radial quantum number of that oscillator, not a GEM index. `beta` is the SHO scale in GeV. `mass` is the hadron mass in GeV. `fl1` and `fl2` are the signed flavor weights of the two quark-line diagrams. `F` multiplies the width and counts charge states that share the same \(|M|^2\). There is no built-in parameter set, so `model.file` is required. Run from the directory that contains that file. Inputs: `test/Decay-3P0/`.

## Parameter files

A preset named `*_CUSTOM` needs `model.file`. The file is a JSON object with one `param` object. Built-in presets that need no file are the `GISCREEN_*` and `GISTRING_*` names other than `*_CUSTOM`, `NRSTRING_MESON`, `NRSCREEN_MESON`, `BAGMIT_DEFAULT`, `QMCSTAD_DEFAULT`, and `YUKAWA_MT13`.

**Godfrey–Isgur** (`GISTRING_CUSTOM`, `GISCREEN_CUSTOM`), all required except `mu`:

| Field | Meaning |
|-------|---------|
| `mn`, `ms`, `mc`, `mb` | Quark masses (GeV) |
| `b`, `c` | String tension and constant |
| `mu` | Screening mass. Missing means 0. Used by the screened set |
| `sigma_0`, `s` | Gaussian smearing |
| `epsilon_cont`, `epsilon_sov`, `epsilon_sos`, `epsilon_tens` | Contact, spin–orbit, Thomas, and tensor strengths |

**Cornell** (`NRSTRING_CUSTOM`, `NRSCREEN_CUSTOM`): `mn`, `ms`, `mc`, `mb`, `b`, `c`, `alpha_s`, `sigma`. `NRSCREEN` also requires `mu` (GeV). `NRSTRING` ignores `mu` and keeps a linear confinement.

**MIT bag** (`BAGMIT_CUSTOM`): `mn`, `ms`, `mc`, `mb` (GeV; the built-in \(m_n\) is 0), `B` (GeV\(^4\)), `Z`, and the pair bindings `Bcs`, `Bbs`, `Bcc`, `Bbc`, `Bbb` (GeV).

**Standard QMC** (`QMCSTAD_CUSTOM`): the vacuum masses, `B`, `Z`, `m_sigma`, `m_omega`, `m_rho`, `g_sigma`, `g_omega`, `g_rho`. `lambda3` is optional and defaults to 0. A file that contains `g_sigma_B` or `delta` is rejected. `QMCSTAD_DEFAULT` copies the vacuum bag, \(m_\sigma=0.550\) GeV, the physical \(\omega\) and \(\rho\) masses, and zero medium couplings. It is not a saturation fit.

**Modified QMC** (`QMCMOD_CUSTOM`): the `QMCSTAD` fields plus `g_sigma_B` \(\neq 0\) and `delta` \(> 0\). There is no modified-QMC preset. `QMCMOD` cannot use `QMCSTAD_DEFAULT`.

**Static OBE** (`OBESTAT_CUSTOM`). Legs, masses, and `I` stay in the system file.

| Field | Meaning |
|-------|---------|
| `lambda` | Common monopole cutoff (GeV) |
| `tensor_mass1`, `tensor_mass2` | Mass \(M\) in \(f/(2M)\), not the kinetic pole mass |
| `g_pi1`, `g_pi2` | Pion couplings. The vertex conversion is internal |
| `g_sigma1`, `g_sigma2` | Scalar couplings |
| `g_rho1`, `g_rho2`, `f_rho1`, `f_rho2` | Vector \(g\) and magnetic \(f=\kappa g\) |
| `g_omega1`, `g_omega2`, `f_omega1`, `f_omega2` | The same for \(\omega\) |
| `cc1`, `cc2` | Meson charge conjugation. Write 0 when unused. A nonzero value on a meson leg flips \(\pi\) and \(\omega\). A Dirac `cc` is ignored |
| `eta` | Meson–meson C-parity. Write \(+1\) or \(-1\). A 0 is read as \(+1\) |

**Yukawa** (`YUKAWA_CUSTOM`): `param.terms`, an array of `{ "V", "mu", "wigner", "bartlett" }`. `V` is in MeV·fm and `mu` in fm\(^{-1}\). `wigner` defaults to 1 and `bartlett` to 0. `YUKAWA_MT13` is the Malfliet–Tjon I/III potential and takes no file.

**EFT contact** (`EFTCONT_CUSTOM`): `Ca` (MeV) and `Rc` (fm),

\[
V(r)=-C_a\exp\left[-(r/R_c)^2\right].
\]

The Gaussian that replaces \(\delta^{(3)}(\mathbf r)\) is already inside \(C_a\). \(C_a<0\) is repulsion. There is no spin operator. One file is one projected channel.

**HAL QCD** (`HALQCD_CUSTOM`): `param.V0`, `param.Vsig`, and `param.VT`. Each present key is an array of `{ "r", "V" }` with `r` in fm and `V` in MeV. `r` is non-negative and strictly increasing, at most 256 points. A missing array is zero, and at least one array is required. The kernel interpolates linearly. \(r\) at or below the first node keeps \(V(r_1)\). \(r\) at or above the last node is zero. A three-body file may instead key those objects by the pair isospin, `"0"`, `"1"`, `"0.5"`, or `"1.5"`. Each recoupled \(t_{ij}\) loads its own block. An unlabelled file is accepted only when that pair has one \(t_{ij}\).

**\(^3P_0\)** (`QPC3P0_CUSTOM`): `gamma`, `mq1`, `mq2`, `mq3`, all masses in GeV and positive. `mq1` is the parent quark, `mq2` the parent antiquark, and `mq3` the created quark.

## Output

All result files are written in the working directory.

| File | Task | Contents |
|------|------|----------|
| `<project>.state.json` | `SPECTRA` | `states`, and the extra root fields listed below |
| `<project>.fit.json` | `FITTING` | `chi2`, `valid`, `edm`, parameters, and per-state masses in MeV |
| `<project>.scatter.json` | `SCATTER` | `mu` (GeV), `a0` (fm), `r0` (fm) |
| `<project>.decay.json` | `3P0DECAY` | `open`, `p`, `EB`, `EC` (GeV), `width` (MeV), and `waves` with `L`, `S`, `M` |
| `<project>.pot.dat` | `SPECTRA` with `print.pot` | Radial potential. Bag and QMC runs do not write it |
| `<project>.wfn.N.dat` | `SPECTRA` with `print.wfn` | State \(N\), counting from 1. Bag and QMC runs do not write it |
| `<project>.basis.dat` | baryon or three-body `SPECTRA` with `print.wfn` | Quantum numbers of each basis function |

A spectrum file does not copy the input back. Meson and two-body `states` have `index`, `mass` (GeV), `rms_radius` (fm), and `eigenvector`. The meson vector has length `nmax`. The two-body vector has length \(n_{\mathrm{ch}}\times n_{\max}\). Baryon states have `mass`, `rms_r12`, `rms_r13`, `rms_r23`, and `eigenvector`. Three-body files add a root field `threshold` (GeV) and, on each state, `ebind` together with the three pair radii and the eigenvector. A fresh three-body file includes that eigenvector. Some copies under `test/Spectra-Molecule3/` omit it, because the full vector makes the file tens to hundreds of megabytes.

Bag files may carry a root string `configuration`. Each state has `st`, `mass`, variational `R` (GeV\(^{-1}\)), an array `x` of four frequencies, and an `eigenvector` of objects `{color, s12, s34, s123, s45, s456, st, st3, c}`. `NONE` stores `{c}` only. A QMC file adds, on the root, `sigma`, `omega`, `B`, `mn_star`, `m_free`, `m_nucleon`, `R_nucleon`, `rho_s`, `mu`, `eps`, `press`, and `ebind`. `ebind` is \(E/A=\varepsilon/\rho_B-M_N\) in GeV, and is 0 at \(\rho_B=0\). `m_nucleon` is the three-light-quark source. It equals a spectrum mass only for an `nnn` bag with basis `NONE`.

Potential columns, from 0.01 fm to 10 fm:

- Meson: `r`, `V`
- Static OBE, Yukawa, and EFT: `r`, `Vc`, `Vss`, `Vt`, `V`
- HAL QCD: `r`, `V0`, `Vsig`, `VT`, and \(V(^1S_0)=V_0-3V_\sigma\), in GeV

A meson `.wfn.N.dat` runs from 0.01 fm to \(\max(10,r_{\max})\) fm and contains `r`, \(R\), \(u=rR\), and \(r^2|R|^2\). An SHO spectrum has no `rmax`, so that upper limit stays 10 fm. A two-body molecule file uses the same grid and writes one radial column for each partial wave. A baryon or three-body `.wfn.N.dat` is the list of expansion coefficients, one line per basis function, and `<project>.basis.dat` names those functions. The baryon and three-body terminal tables print the first 20 states and leave the rest in `.state.json`.

The fit terminal prints `call N  chi2 = ...` on every evaluation, then `chi2`, `valid`, and `edm`. The scatter terminal repeats `mu`, `a0`, and `r0`. Below a decay threshold, `open` is 0 and `width` is 0.

## Solvers

### Quark-model mesons and baryons

`GISCREEN` and `GISTRING` are the Godfrey–Isgur Hamiltonian: relativistic kinematics, linear or screened confinement, a running Coulomb term, and smeared contact, spin–orbit, and tensor interactions. `NRSTRING` and `NRSCREEN` replace that with a non-relativistic Cornell potential and a constant \(\alpha_s\).

The radial problem is the generalized eigenproblem \(Hc=ENc\). Overlap and kinetic matrix elements are analytic. Potential matrix elements are integrated numerically. The symmetric eigenproblem is solved by the routines in `src/math/eigen.c`.

A baryon uses three Jacobi charts and the same GI operator on each pair. Three identical quarks are reduced to one chart. That reduction is not a full \(S_3\) projection, and the three pair radii need not be equal.

### MIT bag and nuclear matter

At a cavity radius \(R\) in GeV\(^{-1}\),

\[
E_0(R)=\sum_f n_f\omega_f(R)+\frac{4\pi}{3}BR^3-\frac{Z}{R},
\qquad
\omega_f=\sqrt{m_f^2+x_f^2/R^2}.
\]

\(x_f\) is the ground-state bag frequency. With `CMI`,

\[
H(R)=E_0(R)\,I-\mathrm{CMI}(R)+\mathrm{Bind}.
\]

The color and spin factor matrices are built once. \(R\) is varied in \((1,10)\), and \(H(R)\) is diagonalized inside one `st` block. `NONE` varies \(E_0(R)\) only.

`QMCSTAD` and `QMCMOD` call the same bag routine after a mean-field solve. The source is always a color-free bag of three light quarks, not the spectrum requested in the input. Only those quarks couple to \(\sigma\). With \(S_q=\partial\omega/\partial m\) and \(g_N=3g_\sigma S_q-(4\pi/3)R^3\,dB/d\sigma\),

\[
m_\sigma^2\sigma+\frac{\lambda_3}{2}\sigma^2=g_N\rho_s,
\qquad
m_n^*=m_n-g_\sigma\sigma,
\qquad
\omega=\frac{3g_\omega\rho_B}{m_\omega^2}.
\]

`QMCSTAD` keeps \(B=B_0\). `QMCMOD` uses the free source mass \(M_N\), the three-quark bag at \(\sigma=0\) and \(B=B_0\):

\[
\frac{B}{B_0}=\left(1-\frac{4g_\sigma^B\sigma}{\delta M_N}\right)^\delta.
\]

A root with \(B\le 0\) is rejected. \(\omega\) does not enter the bag mass or the boundary condition. Symmetric matter sets \(\rho_3=0\), so `g_rho` does not shift the spectrum. Temperature enters only the nucleon Fermi–Dirac sea, with degeneracy 4. It does not heat the cavity. Net density uses \(f-\bar f\). Scalar density, energy, and pressure use \(f+\bar f\). With \(U=\tfrac12 m_\sigma^2\sigma^2+(\lambda_3/6)\sigma^3\),

\[
\varepsilon=\varepsilon_{\mathrm{kin}}+U+\tfrac12 m_\omega^2\omega^2,
\qquad
P=P_{\mathrm{kin}}-U+\tfrac12 m_\omega^2\omega^2,
\qquad
E/A=\varepsilon/\rho_B-M_N.
\]

The spectrum is then the bag at \(m_n^*\), \(B(\sigma)\), and the vacuum \(m_s\), \(m_c\), \(m_b\), and \(Z\). Pair bindings are zero. That mass is not \(M_N+3g_\omega\omega\).

### Molecules

The two-body static one-boson-exchange Hamiltonian is

\[
H=m_1+m_2+\frac{p^2}{2\mu}+V_C+V_{ss}+V_T,
\]

with \(V_{\ell s}=0\). The tensor mixes \(L\) with \(L+2\). Elastic \(0^-\) pairs have \(V_{ss}=V_T=0\) and no pion. Binding reported by the program is \(E_{\mathrm{bind}}=M-(m_1+m_2)\). Treat a state as bound when \(E_{\mathrm{bind}}<0\) and the RMS radius is clearly smaller than `rmax`.

A Yukawa term is \((V/r)e^{-\mu r}(W+B P_\sigma)\) with \(P_\sigma=(1+\sigma_1\cdot\sigma_2)/2\). The EFT contact is the Gaussian written above, and it does not mix \(L\) with \(L+2\). The HAL QCD kernel is

\[
V(r)=V_0(r)+V_\sigma(r)\,(\sigma_1\cdot\sigma_2)+V_T(r)\,S_{12}.
\]

\(^1S_0\) therefore sees \(V_0-3V_\sigma\). There is no \(\mathbf L\cdot\mathbf S\) term and no coupled-channel potential matrix. Each table is already one isospin.

The three-body Hamiltonian is \(\sum_i m_i+\sum_i p_i^2/(2m_i)\) plus the three pair potentials. Binding is \(M-(m_1+m_2+m_3)\).

### Fit, scattering length, and \(^3P_0\)

A fit minimizes

\[
\chi^2=\sum_i\frac{(M_{\mathrm{calc}}^i-M_{\mathrm{exp}}^i)^2}{\sigma_i^2}
\]

with Minuit2 MIGRAD. The masses in the datasets and in `<project>.fit.json` are in MeV. Each \(\chi^2\) evaluation solves the meson Hamiltonian on the requested basis.

`SCATTER` integrates the \(\ell=0\) radial equation with Numerov from 0 to 80 fm, step 0.01 fm, and \(\hbar c=197.3\) MeV·fm. It does not use a GEM eigenvector. The effective-range fit is

\[
k\cot\delta_0=\frac{1}{a_0}+\frac12 r_0 k^2
\]

for ten momenta from 0.05 to 0.50 MeV. \(a_0>0\) is attraction. If the channel list also contains \(L=2\), only the diagonal \(\ell=0\) potential is integrated. A Coulomb tail stops the run.

The decay width is

\[
\Gamma=2\pi\frac{E_B E_C}{M_A}\,p\,F\sum_{LS}|M_{LS}|^2.
\]

Momenta and energies are in GeV. The factor 1000 converts \(\Gamma\) to MeV. Each meson is one SHO of the `beta` and `n` given in the input. A GEM eigenvector is not summed into \(M_{LS}\). The analytic radial sum can lose digits at high \(n\) or large \(\beta\), so a high bottomonium width from this task should be checked by quadrature before it is quoted.

## Test suite

Each directory is a complete input, and the matching result file is already in the tree. Run `--compute` from a directory whose relative `model.file` names resolve.

| Directory | What it checks |
|-----------|----------------|
| `test/Spectra-Quarkonium/` | Built-in `GISCREEN_CCBAR` and `GISCREEN_BBBAR` |
| `test/Spectra-CustomParameter/` | `GISCREEN_CUSTOM` and `GISTRING_CUSTOM` |
| `test/Spectra-SHO/` | Meson SHO basis |
| `test/Spectra-NRmodel/` | Cornell GEM and SHO |
| `test/Spectra-PrintPotWfn/` | Meson `.pot.dat` and `.wfn.N.dat` (`diamond.json`) |
| `test/Spectra-Baryon/` | GI baryons on GEM |
| `test/Spectra-Bag/` | MIT bag, `NONE` and `CMI`, through hexaquarks |
| `test/Spectra-QMC/` | `QMCSTAD` and `QMCMOD`, including a warm run and a CMI nucleon |
| `test/Spectra-Molecule2/` | Two-body static OBE |
| `test/Spectra-Molecule3/` | Three-body `HYBRID` |
| `test/Spectra-Yukawa/` | Malfliet–Tjon deuteron and triton |
| `test/Spectra-EFT/` | Yuan–Liu contact, two- and three-body |
| `test/Spectra-HALQCD/` | Table reader, including an isospin-keyed three-body file |
| `test/ScatterTask/` | S-wave effective range |
| `test/Decay-3P0/` | \(\rho\to\pi\pi\) and open-charm charmonium widths |
| `test/FittingTask/` | `GISCREEN_CSBAR` input. The `ScreenFitting-*` logs are earlier fit records |

Shipped checks, in program units:

- `test/Decay-3P0/rho_pipi.decay.json` records a width of 79.429 MeV.
- `test/ScatterTask/dd_scatter.scatter.json` records \(a_0=-0.250006\) fm. `lclc_scatter.scatter.json` records \(a_0=-0.210000\) fm.
- `test/Spectra-Yukawa/nn_yukawa.state.json` is the MT-I/III deuteron, binding about \(-2.231\) MeV.
- `test/Spectra-QMC/nucleon_qmc.json` is `QMCSTAD_DEFAULT` at \(0.16\,\mathrm{fm}^{-3}\) and \(T=0\). With every medium coupling equal to 0, \(\sigma=0\) and \(E/A\) is the Fermi kinetic energy of the free bag.

The HAL QCD sample `param_HAL_V0.json` stores the same Gaussian as the \(DD\) contact. It checks the table reader. It is not a lattice HAL QCD potential.

## Source layout

```text
gemstore/
├── Makefile
├── include/gemstore/     public headers
├── src/                  parse.c, write.c, print.c, entry.c
│   ├── basis/            GEM, SHO, color, spin, CMI
│   ├── math/             quadrature, matrices, eigenproblems
│   ├── model/            spectra, QMC, OBE, Yukawa, EFT, HAL QCD, scatter, 3P0
│   └── param/            built-in parameters, fit datasets, Minuit2 glue
├── lib/Minuit2/          CERN Minuit2, LGPL-2.1
├── app/gemstore-assistant/
└── test/
```

`src/parse.c` reads the run file. `src/param/argset.c` loads the coefficient set. `src/model/compute.c` dispatches the task. `src/write.c` writes the JSON and text dumps. `src/print.c` writes the terminal tables.

## Optional AI agent skill

`app/gemstore-assistant/` is an optional agent skill: `SKILL.md`, JSON templates, and `scripts/generate_meson_inputs.py`. It prepares a parser file, runs `gemstore --compute`, and reads the result JSON. The solver does not call it, and a calculation does not need an agent. The skill text is distributed under the MIT license. The program itself remains under GPL-3.0-or-later.

## License and citation

GEMSTORE is free software under the GNU General Public License, version 3 or any later version. See `LICENSE`. Copyright © 2026 Wen-Xuan Zhang. The bundled Minuit2 library is under the GNU Lesser General Public License, version 2.1. See `lib/Minuit2/LICENSE`.

The Computer Physics Communications citation will replace this software record once a DOI is assigned:

```bibtex
@software{gemstore2026,
  author = {Zhang, Wen-Xuan},
  title  = {GEMSTORE: Hadron Spectroscopy Simulation Tools},
  year   = {2026},
  url    = {https://github.com/serialcore/gemstore},
  license = {GPL-3.0-or-later}
}
```

Author: Wen-Xuan Zhang, [serialcore@outlook.com](mailto:serialcore@outlook.com).
