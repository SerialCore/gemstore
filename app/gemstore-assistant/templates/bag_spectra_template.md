# Gemstore Bag Spectra JSON Template

Matches `src/parse.c`. Bag `SPECTRA` uses `system.type` `BAG` and `model.type` `BAGMIT`. Basis is `NONE` (bag energy only) or `CMI` (color–spin mixing). Examples: `test/Spectra-Bag/`.

## NONE (no color-magnetic term)

Quark number is `nn+ns+nc+nb` and is not capped.

```json
{
  "project": "baryon_nnn_none",
  "task": "SPECTRA",
  "model": {
    "type": "BAGMIT",
    "param": "BAGMIT_DEFAULT"
  },
  "system": {
    "type": "BAG",
    "nn": 3,
    "ns": 0,
    "nc": 0,
    "nb": 0
  },
  "basis": {
    "type": "NONE"
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## CMI (color–spin mixing)

`nn,ns,nc,nb` are flavor counts. `f1`–`f6` are the quark-order slots used to build CMI (flavor index 1=n, 2=s, 3=c, 4=b). They must match the counts. CMI basis exists for 2–6 quarks. List every color and spin channel the space needs; the program mixes only within the same `st` and applies Pauli selection for five and six quarks.

```json
{
  "project": "tetra_ccnn_18",
  "task": "SPECTRA",
  "model": {
    "type": "BAGMIT",
    "param": "BAGMIT_DEFAULT"
  },
  "system": {
    "type": "BAG",
    "nn": 2,
    "ns": 0,
    "nc": 2,
    "nb": 0,
    "f1": 3,
    "f2": 3,
    "f3": 1,
    "f4": 1
  },
  "basis": {
    "type": "CMI",
    "color": ["tetra11", "tetra88"],
    "spin": [
      {"s12": 1, "s34": 1, "st": 2, "st3": 2},
      {"s12": 1, "s34": 1, "st": 1, "st3": 1},
      {"s12": 1, "s34": 1, "st": 0, "st3": 0},
      {"s12": 1, "s34": 0, "st": 1, "st3": 1},
      {"s12": 0, "s34": 1, "st": 1, "st3": 1},
      {"s12": 0, "s34": 0, "st": 0, "st3": 0}
    ]
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## Fields

- `model.type`: `BAGMIT`
- `model.param`: `BAGMIT_DEFAULT` or `BAGMIT_CUSTOM` (then `model.file` is required)
- `system.type`: `BAG`
- `nn`, `ns`, `nc`, `nb`: integer flavor counts (required)
- `f1`–`f6`: quark order for CMI; omit for `NONE`
- `basis.type`: `NONE` or `CMI`
- `color`: 1–8 names from `meson`, `baryon`, `tetra11`, `tetra33`, `tetra66`, `tetra88`, `penta688`, `penta388`, `penta11`, `hexa1010`, `hexa6868`, `hexa6838`, `hexa3868`, `hexa3838`, `hexa11`
- `spin`: 1–20 channels. Unused couplings stay 0. `st3` does not change the mass; set it equal to `st`
  - meson: `st`, `st3`
  - baryon: `s12`, `st`, `st3`
  - tetra: `s12`, `s34`, `st`, `st3`
  - penta: `s12`, `s123`, `s45`, `st`, `st3`
  - hexa: `s12`, `s123`, `s45`, `s456`, `st`, `st3`
- `print.pot` / `print.wfn`: bag spectra ignore these

## Mass

Units are GeV; bag radius \(R\) is in GeV\(^{-1}\). Ground-state bag energy is \(\sum_f n_f\omega_f+(4\pi/3)BR^3-Z/R\), with \(\omega=\sqrt{m^2+x^2/R^2}\). With `CMI`, \(H(R)=E_0(R)I-\mathrm{CMI}+\mathrm{Bind}\). \(R\) is minimized in \((1,10)\). Color and spin factor matrices are built once; mixing stays inside one `st`.

`BAGMIT_DEFAULT` (`src/param/argset.c`): \(m_n=0\), \(m_s=0.279\), \(m_c=1.641\), \(m_b=5.093\), \(B=4.42051\times10^{-4}\,\mathrm{GeV}^4\), \(Z=1.83\), and pair bindings \(B_{cs},B_{bs},B_{cc},B_{bc},B_{bb}\).

Custom parameter file:

```json
{
  "param": {
    "mn": 0.0,
    "ms": 0.279,
    "mc": 1.641,
    "mb": 5.093,
    "Bcs": -0.0246,
    "Bbs": -0.0319,
    "Bcc": -0.0772,
    "Bbc": -0.1013,
    "Bbb": -0.1279,
    "B": 0.000442051,
    "Z": 1.83
  }
}
```

## Output

`<project>.state.json` has `configuration` (`q`/`Q` string; omitted for `NONE`) and `states` with `index`, `st`, `mass`, `R`, `x` \([x_n,x_s,x_c,x_b]\), and `eigenvector` entries `{color,s12,s34,s123,s45,s456,st,st3,c}`. For `NONE`, each eigenvector entry is `{c}`.
