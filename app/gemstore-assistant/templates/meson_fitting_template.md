# Gemstore Meson Fitting JSON Template

A fit is `gemstore --compute FILE` with `task` set to `FITTING`.

```json
{
  "project": "fit_ccbar",
  "task": "FITTING",
  "fit": { "target": "GISCREEN_CCBAR" },
  "basis": { "type": "GEM", "nmax": 20, "rmax": 30.0, "rmin": 0.1 }
}
```

## Fields

- Do not write `model`, `system`, or `print`. The pinned `fit.target` sets the model and the meson dataset. `print.pot` and `print.wfn` stay off.
- `basis` is `GEM` (`nmax`, `rmax`, `rmin`) or `SHO` (`nmax`, `beta`). The GEM values above are the ones used with these datasets.
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
