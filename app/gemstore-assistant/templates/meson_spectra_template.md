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

- `task`: `SPECTRA`, `DECAY3P0`, `COUPLCHN`, `SCATTER`
- `system.type`: `MESON`, `BAYRON`
- `model.type`: `GISTRING`, `GISCREEN`
- `model.param`:
  - `GISTRING_MESON`
  - `GISTRING_CUSTOM`
  - `GISCREEN_MESON`
  - `GISCREEN_BBBAR`
  - `GISCREEN_BCBAR`
  - `GISCREEN_BSBAR`
  - `GISCREEN_CCBAR`
  - `GISCREEN_CSBAR`
  - `GISCREEN_CUSTOM`
- `basis.type`: `GEM`, `CRG`, `SHO`
- `print.pot`, `print.wfn`: `"true"` or `"false"` (controls output of `.pot.dat` and `.wfn.N.dat` files)

## Custom Parameter Files

If `model.param` is `GISTRING_CUSTOM` or `GISCREEN_CUSTOM`, also provide:

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
- `CRG`: `nmax`, `rmax`, `rmin`, `omega`
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
