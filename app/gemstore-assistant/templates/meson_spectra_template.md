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
  }
}
```

## Allowed Values

- `task`: `SPECTRA`, `DECAY3P0`, `COUPLCHN`, `SCATTER`
- `system.type`: `MESON`
- `model.type`: `GISTRING`, `GISCREEN`
- `model.param`:
  - `GISTRING_MESON`
  - `GISTRING_CUSTOM`
  - `GISCREEN_MESON`
  - `GISCREEN_CCBAR`
  - `GISCREEN_BBBAR`
  - `GISCREEN_CUSTOM`
- `basis.type`: `GEM`, `CRG`, `CSM`, `SHO`

## Custom Parameter Files

If `model.param` is `GISTRING_CUSTOM` or `GISCREEN_CUSTOM`, also provide:

```json
"file": "app/param_GISTRING.json"
```

or

```json
"file": "app/param_GISCREEN.json"
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
    "file": "app/param_GISCREEN.json"
  },
  "basis": {
    "type": "GEM",
    "nmax": 16,
    "rmax": 30.0,
    "rmin": 0.1
  }
}
```

## Basis Parameters

- `GEM`: `nmax`, `rmax`, `rmin`
- `CRG`: `nmax`, `rmax`, `rmin`, `omega`
- `CSM`: `nmax`, `rmax`, `rmin`, `theta`
- `SHO`: `beta`

## Meson Quantum Numbers

- `f1`, `f2`: quark flavors, where `1=n`, `2=s`, `3=c`, `4=b`
- `S`: total spin
- `L`: orbital angular momentum
- `J`: total angular momentum
