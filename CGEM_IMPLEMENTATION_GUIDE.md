# Complex-Range Gaussian (CGEM) Basis Implementation Guide

## Current State

The codebase already has a placeholder `ORBIT_CGEM` enum in `include/gemstore/types.h:12`, but it is **entirely unimplemented**. The infrastructure for complex arithmetic is already present (via `<complex.h>`, `orbit_wfn_complex_t`, and `integral_matrix_element_complex`), but it's currently only used for the p-space Fourier-transform phase factor `(-i)^l`, not for genuinely complex exponents.

---

## What Complex-Range Gaussians Are

Instead of real exponent `ν_n ∈ ℝ` in `exp(−ν_n r²)`, you use **complex** exponents:
> ν_n = |ν_n| · e^{iθ_n} ∈ ℂ

This gives basis functions that are oscillatory-Gaussian, useful for resonances (Gamow states) via the Complex Scaling Method (CSM).

---

## Files to Modify and What to Change

### 1. `include/gemstore/param/argset.h` — Add complex range parameters

Add new fields to `argsInput_t` (alongside existing `nmax`, `rmax`, `rmin`):
```c
double theta;       // complex rotation angle (radians), e.g. for ν_n → ν_n * e^{2iθ}
orbit_type_t orbit; // GEM / CGEM / SHO selector
```

And modify `argsOrbit_t` to carry a complex scale:
```c
double complex scale; // was: double scale
```

### 2. `include/gemstore/basis/orbit.h` — Add complex-ν `getnu` variant

Add a new `getnu_complex()` function that returns `double complex`:
```c
static inline double complex getnu_complex(int n, int nmax, double rmax, double rmin, double theta)
{
    double nu_real = getnu(n, nmax, rmax, rmin);
    return nu_real * cexp(2.0 * I * theta);   // complex rotation
}
```

### 3. `src/basis/orbit.c` — Add complex-exponent basis functions

Add new variants `CGRnlr` and `CGRnlp` where the `scale` parameter is `double complex`:
```c
double complex CGRnlr(double r, int n, int l, double complex nu);
double complex CGRnlp(double p, int n, int l, double complex nu);
```

These return `double complex` because the pre-factor involves `ν^{l/2+3/4}` which is now complex.

Also update the typedef in `orbit.h`:
```c
typedef double complex (*orbit_wfn_cgem_t)(double x, int n, int l, double complex scale);
```

### 4. `src/math/integral.c` — Add CGEM matrix element integrator

Add a new function `integral_matrix_element_cgem()` that uses `double complex` throughout — `conj()` on the bra, complex accumulation, and returns `double complex` (the Hamiltonian becomes complex-symmetric under CSM):
```c
double complex integral_matrix_element_cgem(
    orbit_wfn_cgem_t wfn, potential_t pot,
    double complex node_factor,
    const argsOrbit_t *args_bra, const argsOrbit_t *args_ket, ...);
```

### 5. `src/model/spectra.c` — Dispatch on orbit type and use complex matrices

- When `input->orbit == ORBIT_CGEM`, call `getnu_complex()` instead of `getnu()` in the basis array setup (lines 44–49)
- Use `integral_matrix_element_cgem()` for all matrix elements
- The Hamiltonian matrices `mT`, `mVcoul`, etc. become `double complex` arrays
- The eigenvalue solver (`eigen.c`) must accept complex Hermitian (or complex-symmetric) matrices — you may need to swap `dsygv_` (real symmetric LAPACK) for `zhegv_` (complex Hermitian) or `zsygv_` (complex symmetric)

### 6. `src/entry.c` — Parse the new `theta` and `orbit` input keys

In the `SECTION_GAUSS` parser (lines 116–121), add:
```c
else if (strcmp(key, "theta") == 0) input->theta = atof(val);
else if (strcmp(key, "orbit") == 0) {
    if (strcmp(val, "GEM")  == 0) input->orbit = ORBIT_GEM;
    if (strcmp(val, "CGEM") == 0) input->orbit = ORBIT_CGEM;
    if (strcmp(val, "SHO")  == 0) input->orbit = ORBIT_SHO;
}
```

### 7. `src/math/eigen.c` (likely) — Add complex eigenvalue solver

If currently using LAPACK's `dsygv_` for real symmetric generalized eigenvalue problem, you'll need to add a path for `zhegv_` (complex Hermitian) or `zsygv_` (complex symmetric) depending on the CSM formulation you use.

---

## Summary Table

| File | Change |
|---|---|
| `include/gemstore/param/argset.h` | Add `theta`, `orbit` to `argsInput_t`; change `scale` to `double complex` in `argsOrbit_t` |
| `include/gemstore/basis/orbit.h` | Add `getnu_complex()`, add `orbit_wfn_cgem_t` typedef |
| `src/basis/orbit.c` | Add `CGRnlr()`, `CGRnlp()` with `double complex nu` |
| `src/math/integral.c` | Add `integral_matrix_element_cgem()` returning `double complex` |
| `src/model/spectra.c` | Dispatch on `ORBIT_CGEM`, use complex matrices, call complex eigensolver |
| `src/entry.c` | Parse `theta` and `orbit` keys in `&GAUSS` section |
| `src/math/eigen.c` | Add `zhegv_`/`zsygv_` LAPACK path for complex eigenvalue problem |

---

## Key Design Decision to Clarify

The CSM can be implemented two ways:

1. **Rotate the basis exponents**: `ν_n → ν_n · e^{2iθ}` (what's described above — only changes the Gaussian parameters, coordinate stays real)
   - Advantage: Minimal changes to potential functions in `gimodel.c`
   - Disadvantage: Only rotates the basis, not the potential

2. **Rotate the coordinate**: `r → r · e^{iθ}` (changes the potential evaluation too — more general but more invasive)
   - Advantage: Full CSM — both basis and potential are rotated
   - Disadvantage: All potential functions need to accept complex `r`

**Which formulation do you want to use?** This affects how much of `gimodel.c` needs to change (option 2 requires all potential functions to accept complex `r`; option 1 does not).

---

## Implementation Order (Recommended)

1. Start with `include/gemstore/param/argset.h` — add the new struct fields
2. Add `getnu_complex()` to `include/gemstore/basis/orbit.h`
3. Implement `CGRnlr()` and `CGRnlp()` in `src/basis/orbit.c`
4. Add `integral_matrix_element_cgem()` to `src/math/integral.c`
5. Modify `src/entry.c` to parse `theta` and `orbit` keys
6. Update `src/model/spectra.c` to dispatch on `ORBIT_CGEM` and use complex matrices
7. Add complex eigenvalue solver to `src/math/eigen.c`
8. Test with a simple input file using `orbit = CGEM` and a small `theta` value
