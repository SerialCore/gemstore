# Gemstore Meson Spectra JSON Template

This template matches the current parser in `src/parse.c`.

```json
{
  "project": "amethyst",
  "task": "SPECTRA",
  "system": {
    "type": "MESON",
    "f1": 3,
    "f2": 3,
    "S": 1,
    "L": 0,
    "J": 1
  },
  "model": {
    "type": "GISCREEN",
    "param": "GISCREEN_CCBAR"
  },
  "basis": {
    "type": "GEM",
    "nmax": 16,
    "rmax": 30.0,
    "rmin": 0.1
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## Allowed Values

- `task`: `SPECTRA`, `DECAY3P0`, `COUPLCHN`, `SCATTER`. A parameter fit uses `FITTING` and `templates/meson_fitting_template.md`.
- `system.type`: `MESON`, `BARYON`
- `model.type`: `GISTRING`, `GISCREEN`, `NRSTRING`, `NRSCREEN`
- `model.param`:
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
- `basis.type`: `GEM`, `SHO`
- `print.pot`, `print.wfn`: `"true"` or `"false"` (controls output of `.pot.dat` and `.wfn.N.dat` files)

## Custom Parameter Files

If `model.param` is `GISTRING_CUSTOM`, `GISCREEN_CUSTOM`, `NRSTRING_CUSTOM`, or `NRSCREEN_CUSTOM`, also provide parameter file.

- `mn`, `ms`, `mc`, `mb` — quark masses (GeV)
- `b` — string tension (GeV²)
- `c` — constant potential (GeV)
- `alpha_s` — constant strong coupling
- `sigma` — contact smearing (GeV)
- `mu` — screening mass (GeV), read only for `NRSCREEN`. `NRSTRING` leaves confinement linear and does not read `mu`.

```json
"file": "param_GISTRING.json"
```

or

```json
"file": "param_GISCREEN.json"
```

Example:

```json
{
  "project": "amethyst",
  "task": "SPECTRA",
  "system": {
    "type": "MESON",
    "f1": 3,
    "f2": 3,
    "S": 1,
    "L": 0,
    "J": 1
  },
  "model": {
    "type": "GISCREEN",
    "param": "GISCREEN_CUSTOM",
    "file": "param_GISCREEN.json"
  },
  "basis": {
    "type": "GEM",
    "nmax": 16,
    "rmax": 30.0,
    "rmin": 0.1
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## Basis Parameters

- `GEM`: `nmax`, `rmax`, `rmin`
- `SHO`: `nmax`, `beta`

## SHO Example

```json
{
  "project": "topaz",
  "task": "SPECTRA",
  "system": {
    "type": "MESON",
    "f1": 3,
    "f2": 3,
    "S": 1,
    "L": 0,
    "J": 1
  },
  "model": {
    "type": "GISCREEN",
    "param": "GISCREEN_CCBAR"
  },
  "basis": {
    "type": "SHO",
    "nmax": 16,
    "beta": 0.8
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## Meson Quantum Numbers

- `f1`, `f2`: quark flavors, where `1=n`, `2=s`, `3=c`, `4=b`
- `S`: total spin
- `L`: orbital angular momentum
- `J`: total angular momentum
