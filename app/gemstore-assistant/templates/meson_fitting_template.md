# Gemstore Meson Fitting JSON Template

A fit is `gemstore --compute FILE` with `task` set to `FITTING`.

```json
{
  "project": "fit_ccbar",
  "task": "FITTING",
  "model": { "type": "GISCREEN" },
  "system": { "type": "MESON" },
  "basis": { "type": "GEM", "nmax": 20, "rmax": 30.0, "rmin": 0.1 },
  "print": { "pot": "false", "wfn": "false" },
  "fit": { "target": "GISCREEN_CCBAR" }
}
```

## Fields

- `model` has `type` only. It must be the model named by `fit.target`.
- `system` has `type` only. Current targets are meson fits, so the type is `MESON`.
- `basis` is `GEM` (`nmax`, `rmax`, `rmin`) or `SHO` (`nmax`, `beta`). The GEM values above are the ones used with these datasets.
- `print` is required. Keep both flags `"false"` for a fit.
- `fit.target` selects the built-in state list and parameter table in `src/param/f*.c`.

## fit.target

- `GISCREEN_MESON`
- `GISCREEN_BBBAR`
- `GISCREEN_BCBAR`
- `GISCREEN_BSBAR`
- `GISCREEN_CCBAR`
- `GISCREEN_CSBAR`

## Run

```bash
gemstore --compute fit_ccbar.json
```

Each χ² evaluation prints `call N  chi2 = ...`. The program then prints the Minuit parameters, `chi2`, `valid`, and `edm`, and writes `<project>.fit.json`. Masses in that file are in MeV. A fixed parameter has error 0.

An example input is `test/FittingTask/giscreen_csbar.json`.
