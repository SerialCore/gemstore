# Gemstore Baryon Spectra JSON Template

Matches `src/parse.c`. Baryon `SPECTRA` uses `GEM` with `GISTRING` or `GISCREEN`.

```json
{
  "project": "nnc_sm_1h2_Pp_L0_jl0",
  "task": "SPECTRA",
  "model": {
    "type": "GISTRING",
    "param": "GISTRING_BARYON"
  },
  "system": {
    "type": "BARYON",
    "f1": 1,
    "f2": 1,
    "f3": 3,
    "J": 0.5,
    "P": 1,
    "sym12": -1,
    "Lmax": 0,
    "jl": 0
  },
  "basis": {
    "type": "GEM",
    "nmax": 6,
    "rmin": 0.2,
    "rmax": 2.0
  },
  "print": {
    "pot": "false",
    "wfn": "false"
  }
}
```

## Fields

- Flavors `f1`, `f2`, `f3`: `1` = n, `2` = s, `3` = c, `4` = b
- `J`, `P` (`+1` or `-1`), `sym12` (`+1` or `-1`), `Lmax`
- Omit `jl` for three Jacobi charts with \(l_\rho+l_\lambda\le L_{\max}\)
- Set `jl` for the `recycle` chart only: \(c=1\), \(l_\rho=0\), \(l_\lambda=L_{\max}\)
- `GISTRING_BARYON` is the `recycle/debug.h` parameter set. `GISTRING_MESON` is a meson fit
- `NRSTRING` and `NRSCREEN` are rejected for baryon spectra
- Output `<project>.state.json`: `mass`, `rms_r12`, `rms_r13`, `rms_r23`, `eigenvector`
- Examples: `test/Spectra-Baryon/`
