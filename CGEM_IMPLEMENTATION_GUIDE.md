# Complex-Range Gaussian (CRG) Basis Implementation Guide

## Current State

The codebase already has a placeholder `ORBIT_CRG` enum in `include/gemstore/types.h:12`, but it is **entirely unimplemented**. The infrastructure for complex arithmetic is already present (via `<complex.h>`, `orbit_wfn_complex_t`, and `integral_matrix_element_complex`), but it's currently only used for the p-space Fourier-transform phase factor `(-i)^l`, not for genuinely complex exponents.

---

## What Complex-Range Gaussians Are

Instead of real exponent `ν_n ∈ ℝ` in `exp(−ν_n r²)`, you use **complex** exponents:
> ν_n = |ν_n| · e^{iθ_n} ∈ ℂ

This gives basis functions that are oscillatory-Gaussian, useful for resonances (Gamow states) via the Complex Scaling Method (CSM).

---

## Files to Modify and What to Change

### 1. `include/gemstore/basis/orbit.h` — Add complex-ν `getnu` variant

Add a new `getnu_complex()` function that returns `double complex`:
```c
static inline double complex getnu_complex(int n, int nmax, double rmax, double rmin, double theta)
{
    double nu_real = getnu(n, nmax, rmax, rmin);
    return nu_real * cexp(2.0 * I * theta);   // complex rotation
}
```

### 2. `src/basis/orbit.c` — Add complex-exponent basis functions

Add new variants `CGRnlr` and `CGRnlp` where the `scale` parameter is `double complex`:
```c
double complex CGRnlr(double r, int n, int l, double complex nu);
double complex CGRnlp(double p, int n, int l, double complex nu);
```

These return `double complex` because the pre-factor involves `ν^{l/2+3/4}` which is now complex.

Also update the typedef in `orbit.h`:
```c
typedef double complex (*orbit_wfn_CRG_t)(double x, int n, int l, double complex scale);
```

### 3. `src/math/integral.c` — Add CRG matrix element integrator

Add a new function `integral_matrix_element_CRG()` that uses `double complex` throughout — `conj()` on the bra, complex accumulation, and returns `double complex` (the Hamiltonian becomes complex-symmetric under CSM):
```c
double complex integral_matrix_element_CRG(
    orbit_wfn_CRG_t wfn, potential_t pot,
    double complex node_factor,
    const argsOrbit_t *args_bra, const argsOrbit_t *args_ket, ...);
```

### 4. `src/model/spectra.c` — Dispatch on orbit type and use complex matrices

- When `input->orbit == ORBIT_CRG`, call `getnu_complex()` instead of `getnu()` in the basis array setup (lines 44–49)
- Use `integral_matrix_element_CRG()` for all matrix elements
- The Hamiltonian matrices `mT`, `mVcoul`, etc. become `double complex` arrays
- The eigenvalue solver must accept complex Hermitian (or complex-symmetric) matrices via `zhegv_` (complex Hermitian) or `zsygv_` (complex symmetric)

---

## Summary Table

### **Option 1: Complex Gaussian Exponent**

| File | Change |
|---|---|
| `include/gemstore/basis/orbit.h` | Add `getnu_complex()`, add `orbit_wfn_CRG_t` typedef |
| `src/basis/orbit.c` | Add `CGRnlr()`, `CGRnlp()` with `double complex nu` |
| `src/math/integral.c` | Add `integral_matrix_element_CRG()` returning `double complex` |
| `src/model/spectra.c` | Dispatch on `ORBIT_CRG`, use complex matrices, call complex eigensolver |

### **Option 2: Hiyama's Basis Decomposition**

| File | Change |
|---|---|
| `include/gemstore/basis/orbit.h` | Add `getnu_hiyama()`, add `HiyamaCos_nlr/p()`, `HiyamaSin_nlr/p()` typedefs |
| `src/basis/orbit.c` | Add `HiyamaCos_nlr()`, `HiyamaSin_nlr()`, and momentum variants (all real) |
| `src/math/integral.c` | Add `integral_matrix_element_hiyama()` for paired cos/sin integration |
| `src/model/spectra.c` | Dispatch on `ORBIT_HIYAMA`, build 2×2 block matrices, use real eigensolver |
| `include/gemstore/types.h` | Add `ORBIT_HIYAMA` enum value |
| `src/entry.c` | Parse `omega` parameter in `&GAUSS` section for Hiyama method |

---

## Design Decision: Non-Rotating Potentials

The project has chosen **not to rotate potentials**. This means:
- All potentials in `gimodel.c` remain **real-valued functions with real coordinates**
- Only the Gaussian basis exponents become complex
- Quadrature integration still uses real coordinate nodes (no complex path deformation)

This choice is justified because:

1. **Potentials are typed for real coordinates only** (gimodel.h:40)
   - All 18 potential functions: `double (*potential_t)(double x, ...)`
   - All use real math: `erf()`, `exp()`, `pow()` on real arguments
   - No complex overloads (`cerf`, `cexp`) exist

2. **Integration always uses real nodes** (integral.c:30-40)
   - Gauss-Legendre quadrature generates 50 real values in [0.0037, 10.8]
   - Integration loop calls `pot(node_factor * nodes[i], ...)` with real `nodes[i]`

3. **Minimal code footprint**: Only `scale` parameter becomes `double complex`

---

## Three Approaches Considered (Option 1 Selected)

### **Option 1: Complex Gaussian Exponent** ✅ **SELECTED**

Form: `exp(-ν · e^{iθ} · r²)` or equivalently `exp(-ν(1 + i·tan(θ)) · r²)`

**Implementation:**
- `ν_real = getnu(n, nmax, rmax, rmin)` (stays real)
- `scale_complex = ν_real · e^{iθ}` for each basis function
- Basis functions `CGRnlr()`, `CGRnlp()` return `double complex`
- Prefactor `ν^{l/2+3/4}` is now complex; use `cpow(scale, exp)` to compute it
- All potentials remain real, called only on real coordinates

**Advantages:**
- Standard Complex Scaling Method for resonances (Gamow states)
- Minimal changes to existing code—only `scale` becomes complex
- Potentials untouched—no need to modify `gimodel.c`
- Already has infrastructure: complex function pointers exist, complex eigenvalue solver ready
- Resonances appear as poles in the complex Hamiltonian

**Why it works:**
The complex rotation of the Gaussian exponent creates basis functions that oscillate and grow at large r, asymptotically matching Gamow states. The real potentials still accurately capture short-range physics where the basis functions have support.

---

### **Option 2: Hiyama's Basis Decomposition** ⭐ **ALTERNATIVE**

Form: Two separate real basis functions:
- `r^l · exp(-ν·r²) · cos(ω·ν·r²)`
- `r^l · exp(-ν·r²) · sin(ω·ν·r²)`

where `ω` is a parameter controlling oscillation frequency.

**Implementation:**
- Add two new function variants: `HiyamaCos_nlr()` and `HiyamaSin_nlr()` (and momentum-space versions)
- Both take real parameters: `n`, `l`, `nu` (real), `omega` (real)
- Both return `double` (real-valued functions)
- No complex arithmetic needed
- Basis functions form a pair that spans the complex plane implicitly

**Advantages:**
- Purely real basis functions—no complex arithmetic in basis evaluation
- Simpler numerical integration (real × real × real)
- All existing infrastructure works: potentials, quadrature nodes, integration unchanged
- Hiyama's method is proven effective in nuclear/hadronic spectroscopy
- Easy debugging with real numbers

**Disadvantages:**
- Requires pairing: cos and sin bases must be treated as coupled degrees of freedom
- Matrix construction is more complex (2×2 block structure for each pair)
- Eigenvalue problem is larger (2N basis functions instead of N)
- Less direct connection to complex scaling theory

**Physical interpretation:**
- `cos` and `sin` components together approximate a complex exponential: `e^{i·ω·ν·r²}`
- The decaying envelope `exp(-ν·r²)` provides convergence
- Oscillations controlled by parameter `ω`

**Connection to Option 1:**
For small `ω·ν·r²`, the decomposition approximates:
```
exp(-ν·e^{iθ}·r²) ≈ exp(-ν·r²) · [cos(ω·ν·r²) + i·sin(ω·ν·r²)]
```
where `θ ≈ ω` in the complex scaling formulation.

---

## Codebase Evidence for Option 1

**Basis functions** (orbit.c:11-26):
- Real Gaussian: `GRnlr(r, n, l, nu)` uses real `nu`, returns `double`
- Complex momentum basis: `GRnlp(p, n, l, nu)` uses real `nu`, returns `complex`
- Both have exponential decay terms **currently commented out** (orbit.c:16, 25)

**Potential interface** (gimodel.h:40, gimodel.c:115-167):
- All potentials: `double (*potential_t)(double x, const argsGIModel_t *args_model, ...)`
- Coulomb: `erf(x / mu)` on real `x`
- Confining: `exp(-x*mu)` on real `x`
- All arithmetic: only real operations

**Integration** (integral.c:30-40, 97-138):
- Quadrature nodes are real: `static const double nodes[QUAD_ORDER] = { ..., 0.0037, ..., 10.8 }`
- Potentials called as: `pot(node_factor * nodes[i], args_model, args_dynmc)` with real argument
- No complex node generation or complex-path integration

**Matrix structure** (eigen.c):
- Already uses LAPACK's `zhegv_` or `zsygv_` for complex Hermitian eigenvalue problems
- Complex Hamiltonian matrices are properly constructed and diagonalized

---

## Implementation Path (Option 1: Complex Gaussian Exponent)

With Option 1, the implementation becomes straightforward:

1. **Scale stays as stored parameter** (already completed in `argset.h`)
   - `argsOrbit_t.scale` is `double complex` (fixed in your version)

2. **getnu_complex()** generates complex scales
   - Real part from `getnu(n, nmax, rmax, rmin)`
   - Rotate by phase: `scale = nu_real · cexp(I · theta)`

3. **Complex basis functions** `CGRnlr()`, `CGRnlp()`
   - Use `cpow(scale, l/2 + 3/4)` instead of `pow(nu, ...)`
   - All other prefactors follow
   - Return `double complex`

4. **Matrix elements** `integral_matrix_element_CRG()`
   - Quadrature still uses real nodes
   - Potentials still called with real arguments
   - Basis functions are complex, result is complex
   - Conjugate the bra for Hermiticity: `conj(wfn_bra) * pot(...) * wfn_ket`

5. **spectra.c dispatch**
   - When `orbit == ORBIT_CRG`, use complex basis and `integral_matrix_element_CRG()`
   - Hamiltonian matrices become `double complex` arrays
   - Eigenvalue solver already handles this via `zhegv_`

**No changes to potentials, quadrature nodes, or integration loops.**

---

## Implementation Order (Remaining Tasks)

### **If pursuing Option 1 (Complex Gaussian Exponent):**

1. Add `getnu_complex()` to `include/gemstore/basis/orbit.h`
2. Implement `CGRnlr()` and `CGRnlp()` in `src/basis/orbit.c`
3. Add `integral_matrix_element_CRG()` to `src/math/integral.c`
4. Modify `src/model/spectra.c` to dispatch on `ORBIT_CRG` and use complex matrices
5. Test with a simple input file using `orbit = CRG` and a small `theta` value

### **If pursuing Option 2 (Hiyama's Basis Decomposition):**

1. Add `ORBIT_HIYAMA` enum to `include/gemstore/types.h`
2. Add `getnu_hiyama()` and function prototypes to `include/gemstore/basis/orbit.h`
3. Implement `HiyamaCos_nlr()`, `HiyamaSin_nlr()`, `HiyamaCos_nlp()`, `HiyamaSin_nlp()` in `src/basis/orbit.c`
4. Add `omega` parameter parsing to `src/entry.c` (in `SECTION_GAUSS` parser)
5. Add `integral_matrix_element_hiyama()` to `src/math/integral.c` (handles paired cos/sin integration)
6. Modify `src/model/spectra.c` to build 2×2 block matrices for Hiyama basis pairs
7. Test with a simple input file using `orbit = HIYAMA` and an `omega` value

### **Comparison for Decision:**

| Criterion | Option 1 (Complex) | Option 2 (Hiyama) |
|-----------|-------------------|-------------------|
| Code complexity | Moderate | Higher (block matrices) |
| Numerical difficulty | Uses complex arithmetic | Pure real arithmetic |
| Basis size | N | 2N (paired cos/sin) |
| Integration cost | O(N²) + complex ops | O((2N)²) pure real |
| Debugging | Harder (complex numbers) | Easier (real numbers) |
| Theory connection | Standard CSM | Proven hadronic method |
| Convergence | Direct complex eigenvalues | Real eigenvalues + pairing |

---

# Next Steps

This guide provides the high-level overview and design decisions. For detailed implementation instructions, see:

- **Option 1 (Complex-Range Gaussian)**: See `OPTION1_CRG_IMPLEMENTATION_GUIDE.md`
- **Option 2 (Hiyama's Basis)**: See `OPTION2_HIYAMA_IMPLEMENTATION_GUIDE.md`
