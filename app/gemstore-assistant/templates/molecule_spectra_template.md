# Gemstore Molecule Spectra JSON Template

Matches `src/parse.c`. Two-body molecular `SPECTRA` uses `system.type` `MOLECUL2`, basis `GEM`, and either `model.type` `OBESTAT` (`OBESTAT_CUSTOM`) or `YUKAWA` (`YUKAWA_MT13` or `YUKAWA_CUSTOM`). Example: `test/Spectra-Molecule2/`. Three-body `MOLECUL3` uses `model.type` `HYBRID`. Each of `model1`/`model2`/`model3` is `OBESTAT` or `YUKAWA`. Example: `test/Spectra-Molecule3/`.

The parameter file holds the OBE couplings for a pair of species. The system JSON holds the state (legs, masses, isospin, C, J).

```json
{
  "project": "dstar_xi_I0",
  "task": "SPECTRA",
  "model": {
    "type": "OBESTAT",
    "param": "OBESTAT_CUSTOM",
    "file": "param_OBESTAT_dstarxi.json"
  },
  "system": {
    "type": "MOLECUL2",
    "type1": "VECTOR",
    "type2": "DIRAC",
    "m1": 2.010,
    "m2": 1.318,
    "cc1": 0,
    "cc2": 0,
    "eta": 0,
    "iso_pi": -3,
    "iso_rho": -3,
    "S": 0.5,
    "L": 0,
    "J": 0.5
  },
  "basis": {
    "type": "GEM",
    "nmax": 20,
    "rmax": 30.0,
    "rmin": 0.1
  },
  "print": {
    "pot": "true",
    "wfn": "false"
  }
}
```

## Allowed Values

- `task`: `SPECTRA`
- `system.type`: `MOLECUL2`
- `model.type`: `OBESTAT`
- `model.param`: `OBESTAT_CUSTOM` (requires `model.file`)
- `basis.type`: `GEM` (`nmax`, `rmax`, `rmin`)
- `print.pot`, `print.wfn`: `"true"` or `"false"`
- `type1`, `type2`: `DIRAC` (baryon 1/2), `VECTOR` (1− meson), `PSEUDOSCALAR` (elastic 0− or 0+)

Do not write a unified `iso` key.

## System JSON (this state)

| Field | Meaning |
|-------|---------|
| `type1`, `type2` | OBE legs |
| `m1`, `m2` | Physical masses (GeV); kinetic \(T=m_1+m_2+p^2/2\mu\) and threshold |
| `cc1`, `cc2` | Meson C. `0` leaves G-parity. Nonzero on a meson leg flips \(\pi,\omega\). Dirac `cc` is ignored. |
| `eta` | MM C-parity on pion and magnetic vector. Write `+1` or `-1`. `0` is read as `+1`. Unused for meson–baryon. |
| `iso_pi`, `iso_rho` | \(\langle\tau_1\cdot\tau_2\rangle\) on \(\pi\) and \(\rho\). Two \(I=1/2\) legs: \(I=0\to-3\), \(I=1\to+1\). Write both. \(\sigma,\omega\) are isoscalar. |
| `S`, `L`, `J` | `L` is the lowest partial wave. The driver adds \(L+2\) and every total spin that couples to \(J\). JSON `S` does not select the channel list. |

Examples: \(J=1/2\) → \(^2S+^4D\); \(J=3/2\) → \(^4S+^2D+^4D\).

## Parameter File (how this pair couples to \(\pi\sigma\rho\omega\))

```json
{
  "param": {
    "lambda": 1.10,
    "tensor_mass1": 1.867,
    "tensor_mass2": 1.318,
    "g_pi1": 0.60,
    "g_pi2": -0.25,
    "g_sigma1": 3.4,
    "g_sigma2": 3.4,
    "g_rho1": 2.6,
    "g_rho2": 2.9,
    "f_rho1": 11.7,
    "f_rho2": -4.35,
    "g_omega1": 2.6,
    "g_omega2": 2.9,
    "f_omega1": 11.7,
    "f_omega2": -4.35
  }
}
```

| Field | Meaning |
|-------|---------|
| `lambda` | Common monopole cutoff (GeV) |
| `tensor_mass1/2` | \(M\) in \(f/(2M)\), Liu Table II convention, not the pole mass |
| `g_*` | Liu table couplings; pion vertex conversion is internal |
| `f_*` | Magnetic \(f=\kappa g\). Omit or `0` to drop magnetic \(V_{ss},V_t\) |

Do not put `iso_pi`, `type1`, `m1`, `cc`, or `eta` in this file.

`model.file` is opened from the process cwd. Run next to the parameter file (as in `test/Spectra-Molecule2/`) or give an absolute path.

## Hamiltonian and output

Static LO OBE: \(H=T+V_C+V_{ss}+V_t\), \(V_{ls}=0\). Tensor mixes \(L\) with \(L+2\).

CLI prints `Mass`, `Ebind = Mass − (m1+m2)`, `RMS`. Negative `Ebind` with RMS clearly below `rmax` is a bound state. Near-threshold RMS \(\approx\) `rmax` is a box state.

- `"print":{"pot":"true"}` writes `<project>.pot.dat`: `r[fm]  Vc  Vss  Vt  V` on 0.01–10 fm with unit spin matrix elements.
- `"print":{"wfn":"true"}` writes `<project>.wfn.N.dat` through `write_mol2_wfn`: one \(R_a(r)\) column per partial wave, from 0.01 fm to \(\max(10, r_{\max})\) fm.
- `<project>.state.json` reuses the meson writer: `mass`, `rms_radius`, `eigenvector` of length `nch * nmax`.

## Three-body (`MOLECUL3`)

```json
{
  "project": "dkkstar",
  "task": "SPECTRA",
  "model": {
    "type": "HYBRID",
    "model1": { "type": "OBESTAT", "param": "OBESTAT_CUSTOM", "file": "param_OBESTAT_dk.json" },
    "model2": { "type": "OBESTAT", "param": "OBESTAT_CUSTOM", "file": "param_OBESTAT_dkstar.json" },
    "model3": { "type": "OBESTAT", "param": "OBESTAT_CUSTOM", "file": "param_OBESTAT_kkstar.json" }
  },
  "system": {
    "type": "MOLECUL3",
    "type1": "PSEUDOSCALAR",
    "type2": "PSEUDOSCALAR",
    "type3": "VECTOR",
    "m1": 1.867,
    "m2": 0.494,
    "m3": 0.892,
    "J": 1.0,
    "P": -1,
    "I": 0.5,
    "Lmax": 2,
    "sym12": 0,
    "pair1": { "iso_pi": -3, "iso_rho": -3, "eta": 0 },
    "pair2": { "iso_pi": -3, "iso_rho": -3, "eta": 0 },
    "pair3": { "iso_pi": 1, "iso_rho": 1, "eta": 0 }
  },
  "basis": { "type": "GEM", "nmax": 6, "rmax": 30.0, "rmin": 0.1 },
  "print": { "pot": "false", "wfn": "false" }
}
```

- Nested models are pairs 12, 13, 23. `nmax` is widths per Jacobi axis; \(n_{\mathrm{full}}=n_{\mathrm{ang}}\times n_{\max}^2\).
- `I` is required total isospin. Pair \(t_{ij}\) is a basis quantum number; \(\pi/\rho\) use \(\langle\tau_i\cdot\tau_j\rangle\) of the operator pair after recoupling. `pair*.iso_pi`/`iso_rho` are parsed but ignored. Optional `t1`,`t2`,`t3` default to \(1/2\). DKK* uses \(I=1/2\).
- `sym12=0` uses three independent charts (distinguishable hadrons).
- `sym12=-1` on three identical Dirac legs builds each component as the rearrangement sum of charts c=1, 2 and 3, with fermion Pauli \((-1)^{l_\rho+s_{ij}+t_{ij}}=-1\).
- CLI `Ebind = Mass − (m1+m2+m3)` and three pair RMS. Bound only if `Ebind < 0` and RMS is clearly below `rmax`.
- `"print":{"pot":"true"}` writes three pair blocks in `<project>.pot.dat`.
- Run from `test/Spectra-Molecule3/` so relative pair files resolve. Use a distinct `project` per scan.
- A nested pair may be `YUKAWA` instead of `OBESTAT`. `YUKAWA_MT13` needs no file. `YUKAWA_CUSTOM` uses the same `terms` file as two-body Yukawa. Example: `test/Spectra-Yukawa/nnn_yukawa.json`.

## Yukawa (`YUKAWA`)

Local potential, not the molecular OBE. Two-body driver is `mol2_yukawa_GEM`. Three-body uses the same Jacobi GEM and picks the pair kernel from the nested model.

```json
"model": { "type": "YUKAWA", "param": "YUKAWA_MT13" }
```

`YUKAWA_CUSTOM` file:

```json
{
  "param": {
    "terms": [
      {"V": 1438.72, "mu": 3.11, "wigner": 1.0, "bartlett": 0.0},
      {"V": -570.4265, "mu": 1.55, "wigner": 1.0, "bartlett": 0.0},
      {"V": -56.4585, "mu": 1.55, "wigner": 0.0, "bartlett": 1.0}
    ]
  }
}
```

`V` is MeV·fm, `mu` is fm⁻¹. `wigner` defaults to 1 and `bartlett` to 0. \(P_\sigma=(1+\sigma_1\cdot\sigma_2)/2\). `YUKAWA_MT13` is that Malfliet–Tjon I/III set. With \(m_N=0.938942\) GeV (\(\hbar^2/m=41.47\) MeV·fm²), the deuteron (`J=1`, `L=0`, `nmax` 20) is −2.231 MeV. The triton (`Lmax=0`, `I=0.5`, `sym12=-1`) uses the local potential in every partial wave. Inputs are in `test/Spectra-Yukawa/`. With `rmin` 0.1 fm and `rmax` 10 fm, `nmax` 7 (`nnn_yukawa.json`) gives −8.67 MeV and the scan levels off at −8.88 MeV by `nmax` 10. The s-wave-projected Faddeev benchmark is −8.54 MeV. `iso_pi` and `iso_rho` are ignored. There is no tensor, so \(L\) and \(L+2\) do not mix.
