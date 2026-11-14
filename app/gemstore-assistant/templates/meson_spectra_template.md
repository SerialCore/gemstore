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

- `task`: `SPECTRA`. A parameter fit uses `FITTING` and `templates/meson_fitting_template.md`.
- `system.type`: `MESON`, `BARYON`. MIT bag uses `BAG` and `templates/bag_spectra_template.md`. Two-body molecules use `MOLECUL2` and `templates/molecule_spectra_template.md`.
- `model.type`: `GISTRING`, `GISCREEN`, `NRSTRING`, `NRSCREEN`. Bag uses `BAGMIT`. Molecules use `OBESTAT`.
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
- `basis.type`: `GEM`, `SHO`. Bag uses `NONE` or `CMI`.
- `print.pot`, `print.wfn`: `"true"` or `"false"` (controls output of `.pot.dat` and `.wfn.N.dat` files)

## Custom Parameter Files

If `model.param` is `GISTRING_CUSTOM`, `GISCREEN_CUSTOM`, `NRSTRING_CUSTOM`, or `NRSCREEN_CUSTOM`, also provide `model.file`.

`GISTRING_CUSTOM` and `GISCREEN_CUSTOM` require `mn`, `ms`, `mc`, `mb`, `b`, `c`, `sigma_0`, `s`, `epsilon_cont`, `epsilon_sov`, `epsilon_sos`, `epsilon_tens`. `mu` is optional and defaults to 0. Do not write `alpha_s` or `sigma`.

`NRSTRING_CUSTOM` and `NRSCREEN_CUSTOM` require `mn`, `ms`, `mc`, `mb`, `b`, `c`, `alpha_s`, `sigma`. `mu` is optional and defaults to 0. `NRSTRING` does not use `mu`. `NRSCREEN` uses `mu` as the screening mass in GeV. Do not write `sigma_0` or `epsilon_*`.

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
