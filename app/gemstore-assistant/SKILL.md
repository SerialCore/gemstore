---
name: gemstore-assistant
description: Expert agent for running hadron spectroscopy simulations using the gemstore program. Handles JSON input generation for the current parser, meson and baryon GEM/SHO spectra, two-body molecular GEM with static OBE (MOLECUL2 + OBESTAT), three-body molecular GEM with HYBRID pair OBESTAT (MOLECUL3), MIT bag spectra (BAG + BAGMIT with NONE or CMI), in-medium QMC on that bag (QMCSTAD and QMCMOD at rho_b and temp), built-in meson parameter fits, basis selection, preset or custom parameter-file models, potential/wavefunction output control, and structured JSON + text outputs.
license: MIT
compatibility: opencode
metadata:
  audience: researchers, hadron physicists, computational particle physics
  domain: hadron spectroscopy, quark models, Gaussian expansion method
  tools: bash, file operations, subprocess execution
  keywords: gemstore, hadron spectroscopy, JSON input, GISCREEN, GISTRING, GISTRING_BARYON, NRSTRING, NRSCREEN, BAGMIT, BAG, QMCSTAD, QMCMOD, QMC_DEFAULT, QMC_CUSTOM, QMC, quark-meson coupling, rho_b, CMI, NONE, MIT bag, tetraquark, pentaquark, hexaquark, meson spectra, baryon spectra, MOLECUL2, OBESTAT, OBE, Dstar Xi, molecular states, FITTING, fit.target, GISCREEN_CCBAR, charmonium, bottomonium, GEM, SHO, custom parameter file, print section, potential output, wavefunction output
---

# Gemstore Hadron Spectra Skill

You are an expert agent specialized in **hadron spectroscopy** using the `gemstore` binary. Your role is to interpret natural language user requests, prepare correct JSON input files for the current parser in `src/parse.c`, execute the program safely, parse the JSON output, and deliver clean, structured results.

## Current Parser Contract

Use JSON input files with `--compute FILE`.

The parser expects:

- top-level `project`
- top-level `task`
- object `basis`
- for `SPECTRA`: object `model`, object `system`, object `print`
- for `FITTING`: object `fit` (`target` only). Do not write `model`, `system`, or `print`.

Use exact uppercase strings where shown below.

## Supported Input Values

### Tasks

- `SPECTRA`
- `FITTING`

### Fitting

Use `templates/meson_fitting_template.md`. A fit uses `gemstore --compute FILE`.

For `task` `FITTING`:

- `fit.target` is one of `GISCREEN_MESON`, `GISCREEN_BBBAR`, `GISCREEN_BCBAR`, `GISCREEN_BSBAR`, `GISCREEN_CCBAR`, `GISCREEN_CSBAR`
- The pinned target sets the model and the meson dataset. Do not write `model` or `system`
- `basis` is `GEM` or `SHO`. Use `nmax` 20, `rmax` 30, `rmin` 0.1 unless the user chooses another basis
- Do not write `print`

The state list and the parameter starts live in `src/param/f*.c`. Do not put experimental masses or starting values in the JSON. Each evaluation prints `call N  chi2 = ...`. Read `<project>.fit.json` for `chi2`, `valid`, `edm`, parameters, and per-state `mass_calc`, `mass_exp`, and `residual`. Masses there are in MeV.

### Systems

- `MESON`
- `BARYON`
- `BAG`
- `MOLECUL2`
- `MOLECUL3` — three-body molecule; `model.type` is `HYBRID` with nested `model1`/`model2`/`model3` (pairs 12/13/23). Basis `GEM`. Example: `test/Spectra-Molecule3/dkkstar.json`

Baryon `SPECTRA` accepts `GISTRING` and `GISCREEN` with basis `GEM`. `NRSTRING` and `NRSCREEN` baryon spectra are not implemented. `SHO` is meson-only. Bag `SPECTRA` uses `BAGMIT`, `QMCSTAD`, or `QMCMOD` with basis `NONE` or `CMI`. Use `templates/bag_spectra_template.md` for `BAGMIT`. QMC adds `rho_b` and `temp` on the same `BAG` system; see the QMC section below. Two-body molecules use `MOLECUL2` + `OBESTAT` + `GEM`. Three-body molecules use `MOLECUL3` + `HYBRID` + `GEM`. Use `templates/molecule_spectra_template.md`.

### Models

- `GISTRING`
- `GISCREEN`
- `NRSTRING` — non-relativistic Cornell model, linear confinement
- `NRSCREEN` — non-relativistic Cornell model, screened confinement
- `BAGMIT` — MIT bag, for `system.type` `BAG` only
- `QMCSTAD` — standard quark-meson coupling on `BAG`. \(B\) stays at \(B_0\). Preset `QMC_DEFAULT` or `QMC_CUSTOM`
- `QMCMOD` — modified QMC on `BAG`. \(B/B_0=(1-4 g_\sigma^B\sigma/(\delta M_N))^\delta\). Preset `QMC_CUSTOM` only
- `OBESTAT` — static elastic one-boson exchange, for `system.type` `MOLECUL2`, or as a nested pair model under `HYBRID`
- `YUKAWA` — local Yukawa (Malfliet–Tjon), for `MOLECUL2` or as a nested pair model under `HYBRID`. NN/NNN binding. Preset `YUKAWA_MT13` or `YUKAWA_CUSTOM`
- `EFTCONT` — leading-order EFT contact \(V=-C_a\exp[-(r/R_c)^2]\) (arXiv:2610.28021 Eq. 40), for `MOLECUL2` or a nested `HYBRID` pair. Preset `EFTCONT_CUSTOM`
- `HYBRID` — three pair models `model1`/`model2`/`model3` (pairs 12/13/23), for `system.type` `MOLECUL3`. Each pair is `OBESTAT`, `YUKAWA`, or `EFTCONT`

Use `model.param` for presets:

- `GISTRING_MESON`
- `GISTRING_BARYON`
- `GISTRING_CUSTOM`
- `GISCREEN_MESON`
- `GISCREEN_BBBAR`
- `GISCREEN_BCBAR`
- `GISCREEN_BSBAR`
- `GISCREEN_CCBAR`
- `GISCREEN_CSBAR`
- `GISCREEN_CUSTOM`
- `NRSTRING_MESON`
- `NRSTRING_CUSTOM`
- `NRSCREEN_MESON`
- `NRSCREEN_CUSTOM`
- `BAGMIT_DEFAULT`
- `BAGMIT_CUSTOM`
- `QMC_DEFAULT`
- `QMC_CUSTOM`
- `OBESTAT_CUSTOM`
- `YUKAWA_MT13`
- `YUKAWA_CUSTOM`
- `EFTCONT_CUSTOM`

For custom parameter sets, the model object must also include:

- `file`: path to the JSON parameter file

- `mn`, `ms`, `mc`, `mb` — quark masses (GeV)
- `b` — string tension (GeV²)
- `c` — constant potential (GeV)
- `alpha_s` — constant strong coupling
- `sigma` — contact smearing (GeV)
- `mu` — screening mass (GeV). Required for `NRSCREEN`. `NRSTRING` does not read `mu`; confinement stays linear.

`BAGMIT_CUSTOM` fields are listed in `templates/bag_spectra_template.md`. `QMC_CUSTOM` fields are `mn`, `ms`, `mc`, `mb`, `B`, `Z`, `m_sigma`, `m_omega`, `m_rho`, `g_sigma`, `g_omega`, `g_rho`, `g_sigma_B`, `delta`, inside `param`. `lambda3` is optional and defaults to 0. There is no pair-binding block. `OBESTAT_CUSTOM` fields are listed in `templates/molecule_spectra_template.md` (`lambda`, `tensor_mass1/2`, `g_*`, `f_*`). `YUKAWA_CUSTOM` is `param.terms[]` with `V` (MeV·fm), `mu` (fm⁻¹), `wigner` (default 1), `bartlett` (default 0). `YUKAWA_MT13` is built in and takes no file. `EFTCONT_CUSTOM` is `param.Ca` (MeV) and `param.Rc` (fm), with \(V=-C_a\exp[-(r/R_c)^2]\). The Gaussian norm is already inside \(C_a\). \(C_a<0\) is repulsion. There is no spin operator. Yuan–Liu Table I: `test/Spectra-EFTCONT/`. Legs and physical masses stay in the system JSON. Yukawa ignores `iso_pi` / `iso_rho`; the parser still requires them on `MOLECUL2`. For NN/NNN use nucleon mass `0.938942` GeV so \(\hbar^2/m = 41.47\) MeV·fm². Deuteron: `J=1`, `L=0`. Triton: `J=0.5`, `P=1`, `I=0.5`, `sym12=-1`, `Lmax=0` (s-wave MT-I/III benchmark). Examples: `test/Spectra-Yukawa/nn_yukawa.json`, `test/Spectra-Yukawa/nnn_yukawa.json`.

Examples:

```json
"model": {
  "type": "GISCREEN",
  "param": "GISCREEN_CUSTOM",
  "file": "param_GISCREEN.json"
}
```

```json
"model": {
  "type": "GISTRING",
  "param": "GISTRING_CUSTOM",
  "file": "param_GISTRING.json"
}
```

```json
"model": {
  "type": "OBESTAT",
  "param": "OBESTAT_CUSTOM",
  "file": "param_OBESTAT_dstarxi.json"
}
```

### Basis Types

- `GEM`
- `SHO`
- `CMI` — bag color–spin space
- `NONE` — bag energy without CMI

Basis-specific required parameters:

- `GEM`: `nmax`, `rmax`, `rmin`
- `SHO`: `nmax`, `beta`
- `CMI`: `color` array and `spin` array (see `templates/bag_spectra_template.md`)
- `NONE`: no extra fields

### Print Control

```json
"print": {
  "pot": "false",
  "wfn": "false"
}
```

- `"true"` → enable output (1)
- `"false"` → disable output (0)
- Controls generation of `.pot.dat` and `.wfn.N.dat` files

### Meson Quantum Numbers

Inside `system` provide:

- `f1`
- `f2`
- `S`
- `L`
- `J`

Flavor mapping:

- `1` = `n`
- `2` = `s`
- `3` = `c`
- `4` = `b`

### Baryon Quantum Numbers

Inside `system` provide `f1`, `f2`, `f3`, `J`, `P`, `sym12`, `Lmax`. Basis must be `GEM`. Use `templates/baryon_spectra_template.md`. Examples live in `test/Spectra-Baryon/`.

- `P`: parity, `+1` or `-1`
- `sym12`: sign under exchange of quarks 1 and 2
- `Lmax`: keep channels with \(l_\rho+l_\lambda\le L_{\max}\) on three Jacobi charts
- `jl`: required by the parser. The basis includes every \(j_l\) that couples to \(J^P\); this field does not select a chart.

Prefer `GISTRING_BARYON` for the smeared linear GI set `argsGIString_baryon` in `src/param/argset.c`. `GISTRING_MESON` is a different fit.

`<project>.state.json` states have `mass`, `rms_r12`, `rms_r13`, `rms_r23`, and `eigenvector`. With `"print":{"wfn":"true"}` the program writes `<project>.basis.dat` and coefficient files `<project>.wfn.N.dat`. `"print":{"pot":"true"}` writes `<project>.pot.dat`.

### Bag Quantum Numbers

Inside `system` provide `nn`, `ns`, `nc`, `nb`. For `CMI` also provide `f1`–`f6` as needed so the ordered flavors match those counts. Use `templates/bag_spectra_template.md`. Examples live in `test/Spectra-Bag/`.

`<project>.state.json` has `configuration` and `states` with `mass`, `R`, `x`, and color–spin `eigenvector` components. Bag runs do not write `.pot.dat` or `.wfn.N.dat`.

### QMC on the bag

`model.type` is `QMCSTAD` or `QMCMOD`. `system.type` stays `BAG`. Add required `rho_b` (fm\(^{-3}\)) and `temp` (GeV). Both must be \(\ge 0\). Basis is `NONE` or `CMI`, with the same flavor and color–spin rules as `BAGMIT`. Examples: `test/Spectra-QMC/`.

```json
"model": { "type": "QMCSTAD", "param": "QMC_DEFAULT" },
"system": { "type": "BAG", "nn": 3, "ns": 0, "nc": 0, "nb": 0, "rho_b": 0.16, "temp": 0.0 }
```

- `QMC_DEFAULT` copies the vacuum MIT masses, \(B=4.42051\times 10^{-4}\) GeV\(^4\), and \(Z=1.83\), with \(m_\sigma=0.550\) GeV, physical \(\omega\) and \(\rho\) masses, and every medium coupling 0. \(\sigma\) stays 0. It is not a saturation fit. Use it only with `QMCSTAD`. `QMCMOD` rejects it because \(\delta=0\) and \(g_\sigma^B=0\)
- `QMC_CUSTOM` needs `model.file`. `QMCSTAD` keeps \(B=B_0\) and ignores `g_sigma_B`. `QMCMOD` needs `g_sigma_B` \(\neq 0\) and `delta` \(> 0\). A finite \(\sigma\) also needs `mn` \(> 0\) when `g_sigma` \(\neq 0\)
- The mean-field source is always a three-light-quark `NONE` bag. The requested `nn,ns,nc,nb` and `CMI` space are evaluated once at the solved \(m_n^*=m_n-g_\sigma\sigma\) and \(B(\sigma)\). \(m_s,m_c,m_b,Z\) stay at their vacuum values. Pair bindings stay 0
- \(\omega=3 g_\omega\rho_B/m_\omega^2\) does not enter the bag mass. Symmetric matter sets \(\rho_3=0\), so `g_rho` does not shift the spectrum
- Temperature enters only the nucleon Fermi–Dirac sea. It does not heat the cavity
- Run from the directory that holds `model.file`. `nucleon_qmc.json` is the zero-coupling preset. `nucleon_qmc_custom.json` plus `param_qmc.json` is the sketch \(m_n=0.40\), \(g_\sigma=4.40\), \(g_\omega=3.95\). `nucleon_mqmc.json` is that sketch with \(g_\sigma^B=2\), \(\delta=4\). `nucleon_qmc_cmi.json` is the color–magnetic nucleon. Those couplings are illustrative

`<project>.state.json` has the bag `states` plus root fields `sigma`, `omega`, `B`, `mn_star`, `m_free`, `m_nucleon`, `R_nucleon`, `rho_s`, `mu`, `eps`, `press`, `ebind`. `R` and `R_nucleon` are in GeV\(^{-1}\). `ebind` is \(E/A=\varepsilon/\rho_B-M_N\) in GeV, and is 0 at \(\rho_B=0\). `m_nucleon` is the three-light-quark source, which equals the spectrum mass only for `nnn` with `NONE`. The CLI adds a `QMC MATTER SUMMARY` before the bag table. QMC writes no `.pot.dat` or `.wfn.N.dat`.

At \(0.16\,\mathrm{fm}^{-3}\) and \(T=0\), the shipped samples give \(E/A=+18.9\) MeV (`QMC_DEFAULT`), \(\approx-16.7\) MeV (`param_qmc.json`), and \(\approx-57.8\) MeV (`param_mqmc.json`, \(B/B_0\approx 0.840\)).

### Two-Body Molecule Quantum Numbers

Inside `system` provide `type1`, `type2`, `m1`, `m2`, `iso_pi`, `iso_rho`, `S`, `L`, `J`. Optional `cc1`, `cc2`, `eta` default to 0. Optional `t1`, `t2` default to \(1/2\); two identical \(I=0\) baryons set them to `0`. Basis must be `GEM`. Model is `OBESTAT` with `param` `OBESTAT_CUSTOM` and `file`, `YUKAWA` with `YUKAWA_MT13` or `YUKAWA_CUSTOM`, or `EFTCONT` with `EFTCONT_CUSTOM` and `file`. Use `templates/molecule_spectra_template.md`. Example: `test/Spectra-Molecule2/dstar_xi_I0.json`. Yukawa NN: `test/Spectra-Yukawa/nn_yukawa.json`. EFT contact: `test/Spectra-EFTCONT/dd_eft.json`.

- `type1`, `type2`: `DIRAC` (baryon 1/2), `VECTOR` (1− meson), `PSEUDOSCALAR` (elastic 0− or 0+)
- `m1`, `m2`: physical masses in GeV (threshold and \(p^2/2\mu\))
- `iso_pi`, `iso_rho`: \(\langle\tau_1\cdot\tau_2\rangle\) on \(\pi\) and \(\rho\). Two \(I=1/2\) legs: \(I=0\to-3\), \(I=1\to+1\). Never emit a unified `iso` key
- `cc1`, `cc2`: meson C; `0` leaves G-parity; nonzero on a meson leg flips \(\pi,\omega\)
- `eta`: MM C-parity; `0` is read as `+1`; unused for meson–baryon
- `L` is the lowest partial wave. The solver adds \(L+2\) and every total spin that couples to \(J\). JSON `S` does not pick the channel list (\(J=1/2\) → \(^2S+^4D\); \(J=3/2\) → \(^4S+^2D+^4D\))

The Hamiltonian is static LO OBE: \(H=T+V_C+V_{ss}+V_t\), \(V_{ls}=0\). CLI `Ebind = Mass − (m1+m2)`. Treat a state as bound only if `Ebind < 0` and RMS is clearly smaller than `rmax`. Near-threshold RMS \(\approx\) `rmax` is a box state. Do not report Liu scattering lengths; this path is a bound GEM spectrum.

`"print":{"pot":"true"}` writes `<project>.pot.dat` as `r  Vc  Vss  Vt  V`. `"print":{"wfn":"true"}` writes per-channel radial functions from 0.01 fm to \(\max(10, r_{\max})\) fm.

### Three-Body Molecule Quantum Numbers

Inside `system` provide `type1`–`type3`, `m1`–`m3`, `J`, `P`, `I`, `Lmax`, `sym12`, and `pair1`/`pair2`/`pair3` (`iso_pi`, `iso_rho`, `eta`). Optional `cc1`–`cc3` default to 0; optional `t1`,`t2`,`t3` default to \(1/2\). Basis must be `GEM`. Model must be `HYBRID`. Each nested pair is `OBESTAT` + `OBESTAT_CUSTOM`, `YUKAWA` + `YUKAWA_MT13` / `YUKAWA_CUSTOM`, or `EFTCONT` + `EFTCONT_CUSTOM`. Examples: `test/Spectra-Molecule3/dkkstar.json`, `test/Spectra-Yukawa/nnn_yukawa.json`, `test/Spectra-EFTCONT/ddd_eft.json`.

- Legs: `DIRAC`, `VECTOR`, `PSEUDOSCALAR`. Intrinsic parity is `+1` for Dirac and `−1` for vector/pseudoscalar; total \(P=P_1 P_2 P_3 (-1)^{l_\rho+l_\lambda}\)
- `nmax` is the number of GEM widths per Jacobi axis. The space is \(n_{\mathrm{ang}}\times n_{\max}^2\). Do not treat `nmax` as a spectroscopic \(n_r\)
- `I` is required total isospin. Pair \(t_{ij}\) is a basis quantum number recoupled to \(I\); pair `iso_pi`/`iso_rho` are parsed but ignored. DKK* uses \(I=1/2\)
- `sym12=0` keeps three independent Jacobi charts (distinguishable hadrons). `sym12=-1` on three identical Dirac legs builds each component as the rearrangement sum of charts c=1, 2 and 3, with fermion Pauli \((-1)^{l_\rho+s_{ij}+t_{ij}}=-1\)
- Hamiltonian is NR \(\sum m_i+\sum p_i^2/(2m_i)\) plus pair \(V_C+V_{ss}+V_t\). Elastic \(0^-\) pairs have \(V_{ss}=V_t=0\) and no pion. \(P\leftrightarrow P^*\) transitions are not in OBESTAT
- CLI `Ebind = Mass − (m1+m2+m3)`. Bound only if `Ebind < 0` and pair RMS is clearly smaller than `rmax`. CLI prints the first 20 states; the rest live in `<project>.state.json`
- Run from the directory that holds the three pair parameter files. Use distinct `project` names so `.state.json` files do not overwrite each other

## Available gemstore CLI

```bash
gemstore [--compute FILE] [--debug UNIT]
```

Use `--compute <file.json>` for spectroscopy and for `FITTING`.

## Input Generation Rules

- Always generate `.json` input files for normal runs.
- Use `templates/meson_spectra_template.md` for meson spectra, `templates/meson_fitting_template.md` for meson fits, `templates/baryon_spectra_template.md` for baryons, `templates/bag_spectra_template.md` for MIT bag spectra, and `templates/molecule_spectra_template.md` for two-body molecules.
- Use `scripts/generate_meson_inputs.py` only as a helper for meson JSON inputs.
- Use `scripts/parse_meson_output.py` to parse and summarize gemstore JSON output files when helpful.
- Do not generate the old `.inp` section-based format unless the user explicitly asks for historical compatibility.

## Output Expectations

The `gemstore` program writes JSON output to `<project>.state.json` with the following structure:

```json
{
  "states": [
    {
      "index": 1,
      "mass": 3.1043137618250256,
      "rms_radius": 0.3248640927951712,
      "eigenvector": [0.4314763615098784, 0.565769603960552, ...]
    },
    {
      "index": 2,
      "mass": 3.670680972467695,
      "rms_radius": 0.5430166015432547,
      "eigenvector": [-0.31343804644812767, -0.28099014148751866, ...]
    }
  ]
}
```

Each entry in `states` contains:

- `index` — state number (1-indexed)
- `mass` — computed meson mass (GeV)
- `rms_radius` — root-mean-square radius (fm)
- `eigenvector` — array of expansion coefficients (meson: length `nmax`; molecule: length `nch × nmax`)

When `"print":{"pot":"true"}` or `"print":{"wfn":"true"}` is set, additional text files are generated:
- Meson: `<project>.pot.dat` (r, V) on 0.01–10 fm; `<project>.wfn.N.dat` (r, R, u, r²|R|²) from 0.01 fm to \(\max(10, r_{\max})\) fm (SHO \(r_{\max}=0\) stays at 10 fm)
- Molecule: `<project>.pot.dat` (r, Vc, Vss, Vt, V) on 0.01–10 fm; `<project>.wfn.N.dat` has one \(R_a\) per partial wave on the same radial grid as the meson wavefunction

## Workflow

### Understand the request

- Identify meson, baryon, bag, QMC-in-medium bag, or two-body molecule, and the flavor or hadron content.
- If the user asks to fit, set `task` to `FITTING` and choose `fit.target` from the six `GISCREEN_*` names. Write `project`, `task`, `basis`, and `fit` only. Do not ask the user for a mass table.
- For a meson spectrum extract `S`, `L`, `J`. For a baryon extract `J`, `P`, `sym12`, `Lmax`, and `jl` if the user fixes it. For a bag extract `nn,ns,nc,nb`, `NONE` versus `CMI`, and the full color and spin lists when mixing. For QMC also extract `QMCSTAD` versus `QMCMOD`, `rho_b`, and `temp`. For a molecule extract the two legs, physical masses, `iso_pi` and `iso_rho`, `J`, and lowest `L`.
- Detect task type.
- Detect requested model. For a fit, take the model from `fit.target`.
- Detect requested basis and any basis-specific parameters.

### Prepare execution

- Build a JSON input file.
- For heavy quarkonia prefer `GISCREEN_CCBAR` or `GISCREEN_BBBAR` when appropriate.
- Use `GISTRING_CUSTOM`, `GISCREEN_CUSTOM`, `NRSTRING_CUSTOM`, `NRSCREEN_CUSTOM`, `BAGMIT_CUSTOM`, `QMC_CUSTOM`, `OBESTAT_CUSTOM`, or `YUKAWA_CUSTOM` only when the user explicitly has external parameter files. `QMCSTAD` may use `QMC_DEFAULT` with no file. `QMCMOD` uses `QMC_CUSTOM` plus a file. `MOLECUL2` uses `OBESTAT_CUSTOM` plus a file, or `YUKAWA` with `YUKAWA_MT13` / `YUKAWA_CUSTOM`.
- For spectra runs, use a dedicated run directory and keep user files untouched unless they explicitly ask you to edit them. For `MOLECUL2`, put the system JSON and the OBE parameter file in that directory; `model.file` is opened from the process cwd.

### Execute safely

- Run `gemstore --compute <file.json>` from the directory that contains a relative `model.file`.
- Capture stdout and stderr.
- Use reasonable timeouts.

### Post-process and present results

- For `SPECTRA`, read `<project>.state.json`. For `FITTING`, read `<project>.fit.json` and report `chi2`, `valid`, `edm`, and the floated parameters.
- Check for new text outputs: `<project>.pot.dat` (potential) and `<project>.wfn.N.dat` (wavefunctions per state) when `"print":{"pot":"true"}` or `"print":{"wfn":"true"}` is set.
- Summarize masses, RMS radii, and (if requested) potential/wavefunction file locations. For bag, report `st`, mass, \(R\), and the leading color–spin components. For QMC, also report \(\sigma\), \(\omega\), \(B\), \(m_n^*\), \(M_N\), and \(E/A\) (`ebind`). For molecules, report `Mass`, `Ebind = Mass − (m1+m2)`, RMS, and the printed channel list (`2S 4D`, …).
- Report eigenvectors when relevant.
- Mention parse or validation errors with the exact offending field if gemstore rejects the input.

## Best Practices

- Confirm parameters before large systematic runs.
- Use the `"print"` section to control output of potential (`.pot.dat`) and wavefunction (`.wfn.N.dat`) files.
- Prefer `GEM` unless the user explicitly asks for `SHO`.
- Use exact parser spellings: `MESON`, `BARYON`, `BAG`, `MOLECUL2`, `GISCREEN`, `GISTRING`, `GISTRING_BARYON`, `NRSTRING`, `NRSCREEN`, `BAGMIT`, `BAGMIT_DEFAULT`, `QMCSTAD`, `QMCMOD`, `QMC_DEFAULT`, `QMC_CUSTOM`, `OBESTAT`, `OBESTAT_CUSTOM`, `CMI`, `NONE`, `NRSTRING_MESON`, `NRSCREEN_MESON`, `GISCREEN_CCBAR`, `DIRAC`, `VECTOR`, `PSEUDOSCALAR`, etc.
- Baryon spectra stay on `GEM`. Do not emit `SHO` or an NR model for a baryon.
- Bag spectra use `system` `BAG` and `model` `BAGMIT`, `QMCSTAD`, or `QMCMOD`. List a complete CMI color and spin set for each `st` that should mix. QMC also requires non-negative `rho_b` and `temp`. Do not put `rho_b` or `temp` on `BAGMIT`. Do not pair `QMCMOD` with `QMC_DEFAULT`.
- Two-body molecules use `system` `MOLECUL2`, basis `GEM`, and `OBESTAT` + `OBESTAT_CUSTOM`, `YUKAWA` + `YUKAWA_MT13` / `YUKAWA_CUSTOM`, or `EFTCONT` + `EFTCONT_CUSTOM`. Three-body molecules use `MOLECUL3` + `HYBRID` + three nested pair models (`OBESTAT`, `YUKAWA`, or `EFTCONT`). Write `iso_pi` and `iso_rho` separately. `EFTCONT` ignores those isospin factors; the file is one projected channel.
- Remember that `model.param` is the parameter-set key for `SPECTRA`.
- For `*_CUSTOM`, always include `model.file`.
- For `FITTING`, `fit.target` is the dataset key. Do not write `model`, `system`, or `print`.

## Example Requests

- "Calculate charmonium 1P state with J=1 using GISCREEN and print both potential and wavefunction"
- "Run bottomonium S and P waves with GEM basis, enable potential output only"
- "Compute charmonium with SHO basis, nmax = 16, beta = 0.8, and do not print potential"
- "Calculate the nnc 1/2+ baryon with GISTRING_BARYON, GEM, nmax = 6"
- "Fit GISCREEN_CSBAR with GEM, nmax 20, rmax 30, rmin 0.1"
- "MIT bag mass for nnn with no CMI"
- "QMC nucleon at 0.16 fm^{-3} and T=0 with QMC_DEFAULT"
- "Modified QMC nnn with g_sigma_B and delta from param_mqmc.json, basis NONE"
- "Bag tetraquark ccnn with tetra11 and tetra88, complete spin set in st = 0, 1, 2"
- "D* Ξ I=0 1/2− molecular GEM with Liu Table II couplings and Λ=1.10 GeV"

When this skill is triggered, generate the exact JSON input, run `gemstore --compute <file>`, and report the resulting physics output cleanly, mentioning any generated `.pot.dat` or `.wfn.N.dat` files.
