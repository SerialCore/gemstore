# Option 1: Complex-Range Gaussian (CRG) Exponent Implementation Guide

## Overview

Option 1 uses complex-valued Gaussian exponents: `exp(-ν·e^{iθ}·r²)` to generate complex basis functions. The entire basis and Hamiltonian matrices become complex, but potentials remain real and are only evaluated at real coordinates.

This is the **standard Complex Scaling Method (CSM)** approach for resonance spectroscopy.

---

## Step 1: Add `getnu_complex()` to `include/gemstore/basis/orbit.h`

**Location:** After the existing `getnu()` function (around line 25)

```c
/**
 * Complex-valued scale parameter for Complex-Range Gaussians (CRG).
 * 
 * @param n        Basis function index (1 to nmax)
 * @param nmax     Number of basis functions
 * @param rmax     Maximum range (fm or natural units)
 * @param rmin     Minimum range
 * @param theta    Complex scaling angle (radians), typically [0, π/4]
 * @return         ν_n · e^{i·θ}, where ν_n is the real scale from getnu()
 */
static inline double complex getnu_complex(int n, int nmax, double rmax, double rmin, double theta)
{
    double nu_real = getnu(n, nmax, rmax, rmin);
    return nu_real * cexp(I * theta);  // Pure rotation by angle theta
}
```

**Notes:**
- `cexp(I * theta)` computes `cos(θ) + i·sin(θ)` with unit magnitude
- The real part of `getnu_complex()` decreases from tight to loose basis (same spectrum as `getnu()`)
- All basis functions share the same rotation angle `θ` across all n

---

## Step 2: Implement `CGRnlr()` and `CGRnlp()` in `src/basis/orbit.c`

**Location:** Add after existing `GRnlp()` function (around line 26)

### Coordinate-space basis (CGRnlr):

```c
/**
 * Complex-Range Gaussian basis function in coordinate space.
 * Form: [prefactor] · r^l · exp(-ν·r²)
 * where ν is now complex.
 *
 * @param r        Coordinate (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Complex scale parameter (from getnu_complex)
 * @return         R_nl(r) as double complex
 */
double complex CGRnlr(double r, int n, int l, double complex nu)
{
    // Prevent unused parameter warnings
    (void)n;
    
    // Compute prefactor: 2^(l/2+5/4) · ν^(l/2+3/4) / √Γ(l+3/2)
    double factor_real = pow(2.0, (double)l / 2.0 + 1.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    
    double complex nu_power = cpow(nu, (double)l / 2.0 + 0.75);
    double complex prefactor = (factor_real / gamma_term) * nu_power;
    
    // Radial part: r^l · exp(-ν·r²)
    double r_power = pow(r, (double)l);
    double complex exponential = cexp(-nu * r * r);
    
    return prefactor * r_power * exponential;
}
```

### Momentum-space basis (CGRnlp):

```c
/**
 * Complex-Range Gaussian basis function in momentum space.
 * Form: (-i)^l · [prefactor] · p^l · exp(-p²/(4ν))
 * where ν is complex.
 *
 * @param p        Momentum (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Complex scale parameter (from getnu_complex)
 * @return         R̃_nl(p) as double complex
 */
double complex CGRnlp(double p, int n, int l, double complex nu)
{
    // Prevent unused parameter warnings
    (void)n;
    
    // Phase factor: (-i)^l = e^{-i·π·l/2}
    double complex phase = cexp(-I * M_PI * (double)l / 2.0);
    
    // Prefactor: 2^(-l/2-1/4) · ν^(-l/2-3/4) / √Γ(l+3/2)
    double factor_real = pow(2.0, -(double)l / 2.0 - 0.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    
    double complex nu_power = cpow(nu, -(double)l / 2.0 - 0.75);
    double complex prefactor = phase * (factor_real / gamma_term) * nu_power;
    
    // Momentum part: p^l · exp(-p²/(4ν))
    double p_power = pow(p, (double)l);
    double complex exponential = cexp(-p * p / (4.0 * nu));
    
    return prefactor * p_power * exponential;
}
```

**Key points:**
- Use `cpow()` for complex power (handles branch cuts correctly)
- Use `cexp()` for complex exponential
- The prefactor involves complex powers because `ν` is complex
- Both functions return `double complex`
- The `(-i)^l` phase in momentum space is naturally expressed as `cexp(-I * M_PI * l / 2.0)`

---

## Step 3: Add `integral_matrix_element_CRG()` to `src/math/integral.c`

**Location:** Add after `integral_matrix_element()` function (around line 138)

```c
/**
 * Compute matrix element for Complex-Range Gaussian basis.
 * <ψ_bra | V(r) | ψ_ket> where both ψ are complex-valued.
 *
 * @param wfn            Function pointer to complex basis function (CGRnlr or CGRnlp)
 * @param pot            Potential function (real-valued, takes real argument)
 * @param args_bra       Basis parameters for bra state
 * @param args_ket       Basis parameters for ket state
 * @param args_model     Potential parameters
 * @param args_dynmc     Dynamic potential parameters
 * @return               Matrix element as double complex
 */
double complex integral_matrix_element_CRG(
    double complex (*wfn)(double x, int n, int l, double complex nu),
    double (*pot)(double x, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc),
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    const argsGIModel_t *args_model,
    const argsGIModelDy_t *args_dynmc)
{
    double complex result = 0.0 + 0.0*I;
    
    // Gauss-Legendre quadrature: integral_0^inf = sum_i w_i * f(node_i)
    for (int i = 0; i < QUAD_ORDER; i++) {
        double node = nodes[i];
        double weight = weights[i];
        
        // Compute basis functions at this quadrature node
        // Bra: conjugate for Hermiticity
        double complex wfn_bra_val = conj(wfn(node, args_bra->n, args_bra->l, args_bra->scale));
        
        // Ket
        double complex wfn_ket_val = wfn(node, args_ket->n, args_ket->l, args_ket->scale);
        
        // Potential is always real-valued at real coordinate
        double pot_val = pot(node, args_model, args_dynmc);
        
        // Integrate: node^2 * wfn_bra * pot * wfn_ket (r² from r dr in spherical coords)
        result += weight * node * node * wfn_bra_val * pot_val * wfn_ket_val;
    }
    
    return result;
}
```

**Key points:**
- All quadrature nodes are real (unchanged from standard integration)
- Potential is only called on real `node` values
- Conjugate the bra function for Hermiticity: `conj(wfn_bra)`
- Result is `double complex`
- The `r²` Jacobian factor from spherical coordinate integration (already in quadrature setup)

---

## Step 4: Modify `src/model/spectra.c` to dispatch on `ORBIT_CRG`

**Location:** Find the basis loop (around lines 44–49) and potential loop (around lines 106–137)

### In the basis setup loop:

**Before:** (for real GEM basis)
```c
for (int i = 0; i < nbasis; i++) {
    basis[i].n = i + 1;
    basis[i].l = l;
    basis[i].scale = getnu(i + 1, nbasis, input->rmax, input->rmin);
}
```

**After:** (add dispatch for CRG)
```c
for (int i = 0; i < nbasis; i++) {
    basis[i].n = i + 1;
    basis[i].l = l;
    
    if (input->orbit == ORBIT_CRG) {
        // Complex-range Gaussian: scale becomes complex
        basis[i].scale = getnu_complex(i + 1, nbasis, input->rmax, input->rmin, input->theta);
    } else {
        // Standard GEM or SHO
        double nu_real = getnu(i + 1, nbasis, input->rmax, input->rmin);
        basis[i].scale = (double complex)nu_real;  // Cast to complex for uniform storage
    }
}
```

### In the matrix element computation loop:

**Before:** (for real GEM basis)
```c
for (int ibra = 0; ibra < nbasis; ibra++) {
    for (int iket = 0; iket < nbasis; iket++) {
        mT[ibra * nbasis + iket] = integral_matrix_element(
            GRnlr, NULL, &basis[ibra], &basis[iket], args_model, args_dynmc);
        
        mVcoul[ibra * nbasis + iket] = integral_matrix_element(
            GRnlr, GIVcoul, &basis[ibra], &basis[iket], args_model, args_dynmc);
        // ... more potentials
    }
}
```

**After:** (add dispatch for CRG)
```c
for (int ibra = 0; ibra < nbasis; ibra++) {
    for (int iket = 0; iket < nbasis; iket++) {
        if (input->orbit == ORBIT_CRG) {
            // Use complex basis functions and integral
            mT[ibra * nbasis + iket] = integral_matrix_element_CRG(
                CGRnlr, NULL, &basis[ibra], &basis[iket], args_model, args_dynmc);
            
            mVcoul[ibra * nbasis + iket] = integral_matrix_element_CRG(
                CGRnlr, GIVcoul, &basis[ibra], &basis[iket], args_model, args_dynmc);
            // ... more potentials
        } else {
            // Use standard real basis functions
            mT[ibra * nbasis + iket] = integral_matrix_element(
                GRnlr, NULL, &basis[ibra], &basis[iket], args_model, args_dynmc);
            
            mVcoul[ibra * nbasis + iket] = integral_matrix_element(
                GRnlr, GIVcoul, &basis[ibra], &basis[iket], args_model, args_dynmc);
            // ... more potentials
        }
    }
}
```

**Key points:**
- When `ORBIT_CRG`, call `getnu_complex()` to generate complex scales
- Use `integral_matrix_element_CRG()` instead of `integral_matrix_element()`
- All matrix element results are `double complex` for CRG
- Hamiltonian matrices (`mT`, `mVcoul`, etc.) automatically become `double complex` arrays
- The eigenvalue solver (`eigen.c`) already handles complex matrices via `zhegv_`

---

## Testing Option 1

Create a test input file `test_crg.inp`:
```
&GAUSS
  nmax = 5
  rmin = 0.1
  rmax = 20.0
  theta = 0.1      # Small rotation angle
  orbit = CRG
&END

&GI_MODEL
  ... (your model parameters)
&END
```

Run and check:
1. Complex eigenvalues appear (should have small imaginary parts for small θ)
2. Real parts approximate the real GEM spectrum
3. No crashes in `cexp()`, `cpow()`, or `conj()` operations

---

## Summary: Option 1 Implementation Tasks

1. Add `getnu_complex()` to `include/gemstore/basis/orbit.h`
2. Implement `CGRnlr()` and `CGRnlp()` in `src/basis/orbit.c`
3. Add `integral_matrix_element_CRG()` to `src/math/integral.c`
4. Modify `src/model/spectra.c` to dispatch on `ORBIT_CRG` and use complex matrices
5. Test with a simple input file using `orbit = CRG` and a small `theta` value
