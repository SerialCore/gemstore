# Option 2: Hiyama's Basis Decomposition Implementation Guide

## Overview

Option 2 uses two separate real basis functions per degree of freedom:
- `f_{cos}(r) = r^l · exp(-ν·r²) · cos(ω·ν·r²)`
- `f_{sin}(r) = r^l · exp(-ν·r²) · sin(ω·ν·r²)`

These form pairs that together span the oscillatory-Gaussian space. All functions are real-valued, leading to larger but purely real matrices and standard eigensolvers.

This approach is proven effective in nuclear and hadronic spectroscopy.

---

## Step 1: Add `ORBIT_HIYAMA` enum to `include/gemstore/types.h`

**Location:** In the `orbit_type_t` enum (around line 12)

```c
typedef enum {
    ORBIT_GEM,       // Gaussian Expansion Method (real exponentials)
    ORBIT_CRG,       // Complex-Range Gaussian (complex exponentials)
    ORBIT_SHO,       // Simple Harmonic Oscillator
    ORBIT_HIYAMA     // Hiyama's basis (cos/sin decomposition)
} orbit_type_t;
```

---

## Step 2: Add basis function declarations to `include/gemstore/basis/orbit.h`

**Location:** After the `CGRnlr` typedef section (if Option 1 is also present)

```c
/**
 * Hiyama's basis function with cosine oscillation.
 * Form: r^l · exp(-ν·r²) · cos(ω·ν·r²)
 */
typedef double (*orbit_wfn_hiyama_cos_t)(double x, int n, int l, double nu, double omega);

/**
 * Hiyama's basis function with sine oscillation.
 * Form: r^l · exp(-ν·r²) · sin(ω·ν·r²)
 */
typedef double (*orbit_wfn_hiyama_sin_t)(double x, int n, int l, double nu, double omega);

/**
 * Get scale parameter for Hiyama basis.
 * Same spectrum as standard getnu() since real part drives the size.
 */
#define getnu_hiyama(n, nmax, rmax, rmin) getnu(n, nmax, rmax, rmin)
```

---

## Step 3: Implement basis functions in `src/basis/orbit.c`

**Location:** Add after existing basis functions

### Coordinate-space Hiyama cosine:

```c
/**
 * Hiyama basis function (cosine component) in coordinate space.
 * Form: [prefactor] · r^l · exp(-ν·r²) · cos(ω·ν·r²)
 *
 * @param r        Coordinate (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Real scale parameter (from getnu)
 * @param omega    Oscillation parameter (dimensionless)
 * @return         Real-valued basis function
 */
double HiyamaCos_nlr(double r, int n, int l, double nu, double omega)
{
    (void)n;  // unused
    
    // Prefactor: 2^(l/2+5/4) · ν^(l/2+3/4) / √Γ(l+3/2)
    double factor_real = pow(2.0, (double)l / 2.0 + 1.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    double prefactor = (factor_real / gamma_term) * pow(nu, (double)l / 2.0 + 0.75);
    
    // Radial part: r^l
    double r_power = pow(r, (double)l);
    
    // Gaussian envelope with oscillation: exp(-ν·r²) · cos(ω·ν·r²)
    double r_sq = r * r;
    double gaussian = exp(-nu * r_sq);
    double oscillation = cos(omega * nu * r_sq);
    
    return prefactor * r_power * gaussian * oscillation;
}
```

### Coordinate-space Hiyama sine:

```c
/**
 * Hiyama basis function (sine component) in coordinate space.
 * Form: [prefactor] · r^l · exp(-ν·r²) · sin(ω·ν·r²)
 *
 * @param r        Coordinate (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Real scale parameter (from getnu)
 * @param omega    Oscillation parameter (dimensionless)
 * @return         Real-valued basis function
 */
double HiyamaSin_nlr(double r, int n, int l, double nu, double omega)
{
    (void)n;  // unused
    
    // Prefactor: 2^(l/2+5/4) · ν^(l/2+3/4) / √Γ(l+3/2)
    double factor_real = pow(2.0, (double)l / 2.0 + 1.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    double prefactor = (factor_real / gamma_term) * pow(nu, (double)l / 2.0 + 0.75);
    
    // Radial part: r^l
    double r_power = pow(r, (double)l);
    
    // Gaussian envelope with oscillation: exp(-ν·r²) · sin(ω·ν·r²)
    double r_sq = r * r;
    double gaussian = exp(-nu * r_sq);
    double oscillation = sin(omega * nu * r_sq);
    
    return prefactor * r_power * gaussian * oscillation;
}
```

### Momentum-space Hiyama cosine:

```c
/**
 * Hiyama basis function (cosine component) in momentum space.
 * Form: (-i)^l · [prefactor] · p^l · [Fourier transform of exp(-ν·r²)·cos(ω·ν·r²)]
 *
 * Note: Full Fourier transform of cos term is complex; using simplified form.
 * For practical use, stick to coordinate-space integration.
 *
 * @param p        Momentum (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Real scale parameter (from getnu)
 * @param omega    Oscillation parameter (dimensionless)
 * @return         Real-valued basis function (simplified)
 */
double HiyamaCos_nlp(double p, int n, int l, double nu, double omega)
{
    (void)n;      // unused
    (void)omega;  // Note: Full Fourier transform is complex; simplified version
    
    // Use real Gaussian approximation (omega effect ignored for now)
    double phase_factor = (l % 2 == 0) ? 1.0 : -1.0;  // (-i)^l = ±1 for even/odd l
    
    double factor_real = pow(2.0, -(double)l / 2.0 - 0.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    double prefactor = phase_factor * (factor_real / gamma_term) * pow(nu, -(double)l / 2.0 - 0.75);
    
    double p_power = pow(p, (double)l);
    double exponential = exp(-p * p / (4.0 * nu));
    
    return prefactor * p_power * exponential;
}
```

### Momentum-space Hiyama sine:

```c
/**
 * Hiyama basis function (sine component) in momentum space.
 * Form: Similar to cosine version; see notes above.
 *
 * @param p        Momentum (real, positive)
 * @param n        Basis function index
 * @param l        Orbital angular momentum
 * @param nu       Real scale parameter (from getnu)
 * @param omega    Oscillation parameter (dimensionless)
 * @return         Real-valued basis function (simplified)
 */
double HiyamaSin_nlp(double p, int n, int l, double nu, double omega)
{
    (void)n;      // unused
    (void)omega;  // Note: Full Fourier transform is complex; simplified version
    
    double phase_factor = (l % 2 == 0) ? 1.0 : -1.0;
    
    double factor_real = pow(2.0, -(double)l / 2.0 - 0.25);
    double gamma_term = sqrt(tgamma((double)l + 1.5));
    double prefactor = phase_factor * (factor_real / gamma_term) * pow(nu, -(double)l / 2.0 - 0.75);
    
    double p_power = pow(p, (double)l);
    double exponential = exp(-p * p / (4.0 * nu));
    
    return prefactor * p_power * exponential;
}
```

**Notes:**
- All functions are real-valued with only real arithmetic
- Coordinate-space versions are exact; momentum-space simplified (use coordinate space for integration)
- Both `cos` and `sin` share the same prefactor and envelope
- The oscillation frequency is controlled by `ω·ν·r²`

---

## Step 4: Parse `omega` parameter in `src/entry.c`

**Location:** In the `SECTION_GAUSS` parser (around lines 116–121)

**Before:**
```c
else if (strcmp(key, "theta") == 0) input->theta = atof(val);
else if (strcmp(key, "orbit") == 0) {
    if (strcmp(val, "GEM")  == 0) input->orbit = ORBIT_GEM;
    if (strcmp(val, "CRG") == 0) input->orbit = ORBIT_CRG;
    if (strcmp(val, "SHO")  == 0) input->orbit = ORBIT_SHO;
}
```

**After:**
```c
else if (strcmp(key, "theta") == 0) input->theta = atof(val);
else if (strcmp(key, "omega") == 0) input->omega = atof(val);
else if (strcmp(key, "orbit") == 0) {
    if (strcmp(val, "GEM")    == 0) input->orbit = ORBIT_GEM;
    if (strcmp(val, "CRG")    == 0) input->orbit = ORBIT_CRG;
    if (strcmp(val, "SHO")    == 0) input->orbit = ORBIT_SHO;
    if (strcmp(val, "HIYAMA") == 0) input->orbit = ORBIT_HIYAMA;
}
```

**Also add to `argsInput_t` in `include/gemstore/param/argset.h`:**
```c
typedef struct argsInput {
    // ... existing fields ...
    double theta;   // CRG rotation angle
    double omega;   // Hiyama oscillation parameter
    orbit_type_t orbit;
} argsInput_t;
```

---

## Step 5: Add `integral_matrix_element_hiyama()` to `src/math/integral.c`

**Location:** Add after `integral_matrix_element_CRG()` function

```c
/**
 * Compute matrix element for Hiyama basis (cos/sin decomposition).
 * All functions and potentials are real-valued.
 *
 * Matrices are built with 2×2 block structure:
 * [ <cos|V|cos>  <cos|V|sin> ]
 * [ <sin|V|cos>  <sin|V|sin> ]
 *
 * @param wfn_bra      Basis function for bra (HiyamaCos_nlr or HiyamaSin_nlr)
 * @param wfn_ket      Basis function for ket (HiyamaCos_nlr or HiyamaSin_nlr)
 * @param pot          Potential function (real-valued)
 * @param args_bra     Basis parameters for bra state (n, l, scale)
 * @param args_ket     Basis parameters for ket state (n, l, scale)
 * @param omega        Oscillation parameter
 * @param args_model   Potential parameters
 * @param args_dynmc   Dynamic potential parameters
 * @return             Matrix element as double (real)
 */
double integral_matrix_element_hiyama(
    double (*wfn_bra)(double x, int n, int l, double nu, double omega),
    double (*wfn_ket)(double x, int n, int l, double nu, double omega),
    double (*pot)(double x, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc),
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    double omega,
    const argsGIModel_t *args_model,
    const argsGIModelDy_t *args_dynmc)
{
    double result = 0.0;
    
    // Gauss-Legendre quadrature
    for (int i = 0; i < QUAD_ORDER; i++) {
        double node = nodes[i];
        double weight = weights[i];
        
        // Basis functions at quadrature node
        double wfn_bra_val = wfn_bra(node, args_bra->n, args_bra->l, 
                                      creal(args_bra->scale), omega);
        double wfn_ket_val = wfn_ket(node, args_ket->n, args_ket->l, 
                                      creal(args_ket->scale), omega);
        
        // Potential at real coordinate
        double pot_val = pot(node, args_model, args_dynmc);
        
        // Integrate: node^2 * wfn_bra * pot * wfn_ket
        result += weight * node * node * wfn_bra_val * pot_val * wfn_ket_val;
    }
    
    return result;
}
```

**Key points:**
- All operations are real (no `complex.h` needed)
- Extract real part of `args_bra->scale` and `args_ket->scale` using `creal()`
- Returns `double` (real matrix element)

---

## Step 6: Modify `src/model/spectra.c` for 2×2 block matrices

**Location:** In the basis setup and matrix element loops

### Basis array setup:

**Before:**
```c
for (int i = 0; i < nbasis; i++) {
    basis[i].n = i + 1;
    basis[i].l = l;
    basis[i].scale = getnu(i + 1, nbasis, input->rmax, input->rmin);
}
```

**After:**
```c
// For Hiyama, we need 2N basis functions (cos and sin pairs)
int basis_size = (input->orbit == ORBIT_HIYAMA) ? 2 * nbasis : nbasis;

for (int i = 0; i < basis_size; i++) {
    int basis_index = i / 2 + 1;  // For Hiyama: pair index
    
    basis[i].n = basis_index;
    basis[i].l = l;
    
    if (input->orbit == ORBIT_HIYAMA) {
        double nu_real = getnu_hiyama(basis_index, nbasis, input->rmax, input->rmin);
        basis[i].scale = (double complex)nu_real;  // Store real as complex
    } else if (input->orbit == ORBIT_CRG) {
        basis[i].scale = getnu_complex(basis_index, nbasis, input->rmax, input->rmin, input->theta);
    } else {
        double nu_real = getnu(basis_index, nbasis, input->rmax, input->rmin);
        basis[i].scale = (double complex)nu_real;
    }
}
```

### Matrix element computation:

**Before:**
```c
for (int ibra = 0; ibra < nbasis; ibra++) {
    for (int iket = 0; iket < nbasis; iket++) {
        mT[ibra * nbasis + iket] = integral_matrix_element( ... );
        // ...
    }
}
```

**After:**
```c
int mat_size = (input->orbit == ORBIT_HIYAMA) ? 2 * nbasis : nbasis;

for (int ibra = 0; ibra < mat_size; ibra++) {
    for (int iket = 0; iket < mat_size; iket++) {
        if (input->orbit == ORBIT_HIYAMA) {
            // Determine which component (cos=0, sin=1)
            int bra_pair_idx = ibra / 2;
            int ket_pair_idx = iket / 2;
            int bra_comp = ibra % 2;  // 0=cos, 1=sin
            int ket_comp = iket % 2;
            
            // Select appropriate basis functions
            double (*wfn_bra)(double, int, int, double, double);
            double (*wfn_ket)(double, int, int, double, double);
            
            wfn_bra = (bra_comp == 0) ? HiyamaCos_nlr : HiyamaSin_nlr;
            wfn_ket = (ket_comp == 0) ? HiyamaCos_nlr : HiyamaSin_nlr;
            
            // Create temporary basis args for this pair
            argsOrbit_t args_bra_tmp = basis[bra_pair_idx * 2];
            argsOrbit_t args_ket_tmp = basis[ket_pair_idx * 2];
            
            // Compute matrix elements (kinetic and potential)
            mT[ibra * mat_size + iket] = (double complex)integral_matrix_element_hiyama(
                wfn_bra, wfn_ket, NULL, &args_bra_tmp, &args_ket_tmp, 
                input->omega, args_model, args_dynmc);
            
            mVcoul[ibra * mat_size + iket] = (double complex)integral_matrix_element_hiyama(
                wfn_bra, wfn_ket, GIVcoul, &args_bra_tmp, &args_ket_tmp, 
                input->omega, args_model, args_dynmc);
            // ... more potentials
        } else {
            // Original logic for GEM, CRG, SHO
            // ...
        }
    }
}
```

**Key points:**
- Matrix size becomes `2N × 2N` for Hiyama vs `N × N` for other bases
- Use modulo (`%`) and integer division (`/`) to map between pair index and cos/sin component
- Each 2×2 block corresponds to one `(n, l)` pair with both cos and sin components
- Matrices remain `double complex` but contain real values for Hiyama

---

## Step 7: Testing Option 2

Create a test input file `test_hiyama.inp`:
```
&GAUSS
  nmax = 5
  rmin = 0.1
  rmax = 20.0
  omega = 0.5      # Oscillation parameter
  orbit = HIYAMA
&END

&GI_MODEL
  ... (your model parameters)
&END
```

Run and check:
1. Matrix size is `2N × 2N` (e.g., 10×10 for nmax=5)
2. All eigenvalues are real
3. Real part of eigenvalues approximates resonance energies
4. No complex arithmetic errors (all real)
5. Convergence behavior depends on `omega` and oscillation details

---

## Summary: Option 2 Implementation Tasks

1. Add `ORBIT_HIYAMA` enum to `include/gemstore/types.h`
2. Add `getnu_hiyama()` and function prototypes to `include/gemstore/basis/orbit.h`
3. Implement `HiyamaCos_nlr()`, `HiyamaSin_nlr()`, `HiyamaCos_nlp()`, `HiyamaSin_nlp()` in `src/basis/orbit.c`
4. Add `omega` parameter parsing to `src/entry.c` (in `SECTION_GAUSS` parser)
5. Add `integral_matrix_element_hiyama()` to `src/math/integral.c` (handles paired cos/sin integration)
6. Modify `src/model/spectra.c` to build 2×2 block matrices for Hiyama basis pairs
7. Test with a simple input file using `orbit = HIYAMA` and an `omega` value
