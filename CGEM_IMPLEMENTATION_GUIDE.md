# Complex Gaussian Expansion Method (CGEM) Implementation Guide

## Terminology Clarification (CRITICAL)

**This document clarifies the CORRECT terminology mapping**:

| Term | Method | Basis Functions | Parameters | Status |
|------|--------|-----------------|-----------|--------|
| **CRG** | Complex-Range Gaussian (Hiyama's method) | cos/sin pairs | `omega` | Recommended 🟢 |
| **CSM** | Complex Scaling Method | complex exponent | `theta` | Infrastructure Ready 🟡 |

**Enum Mapping** (`include/gemstore/types.h`):
```c
typedef enum orbit_type {
    ORBIT_GEM,       // Standard Gaussian (real exponents)
    ORBIT_CRG,       // Complex-Range Gaussian (Hiyama's cos/sin decomposition)
    ORBIT_SHO,       // Spherical Harmonic Oscillator (real exponents)
    ORBIT_CSM        // Complex Scaling Method (complex exponent rotation)
} orbit_type_t;
```

**This guide documents implementation for BOTH methods.**

---

## Current Codebase Status (April 2026, commit 7b87d6d)

### ✅ Complete Infrastructure
- Gaussian basis functions: GEM and SHO (both real exponent based)
- Real matrix operations: `matrix.c` (319 lines)
- Real eigenvalue solvers: `eigen.c` (694 lines)
- Gaussian quadrature integration: `integral.c` (138 lines)
- GI-Screen and GI-String potentials: `gimodel.c` (325 lines)
- JSON input/output: `parse.c` (302 lines), `print.c` (385 lines)

### 🟡 Complex Infrastructure (Ready, not yet integrated with new methods)
- Complex matrix operations: `cmatrix.c` (322 lines, ✅ COMPLETE)
- Complex eigenvalue solver: `ceigen.c` (486 lines, ✅ COMPLETE)
- Complex integration: `integral_matrix_element_complex()` exists

### 🔴 Gap for CRG Implementation
- **Missing**: `HiyamaCos_nlr()`, `HiyamaSin_nlr()` basis functions
- **Missing**: `HiyamaCos_nlp()`, `HiyamaSin_nlp()` momentum space variants
- **Missing**: Integration wrapper for paired cos/sin basis
- **Missing**: 2×2 block matrix construction in spectra.c
- **Total gap**: ~200 lines of code for CRG

### 🔴 Gap for CSM Implementation
- **Missing**: `getnu_complex()` function
- **Missing**: `CGRnlr()` (CSM in real space) - despite the confusing naming!
- **Missing**: `CGRnlp()` (CSM in momentum space)
- **Missing**: Complex Hamiltonian dispatch in spectra.c
- **Total gap**: ~150 lines of code for CSM

---

## Method Comparison

### CRG (Hiyama's Method) - RECOMMENDED FOR NEAR-TERM

**Physical Form**:
```
Basis pair: 
  Cos: r^l · exp(-ν·r²) · cos(ω·ν·r²)
  Sin: r^l · exp(-ν·r²) · sin(ω·ν·r²)
```

**Advantages**:
- ✅ Purely real arithmetic (no complex numbers)
- ✅ All existing potentials work unchanged
- ✅ Standard quadrature integration unchanged
- ✅ Proven effective in hadronic spectroscopy (references: Hiyama et al. 1995-2010)
- ✅ Easy debugging with real numbers

**Disadvantages**:
- Matrix size: 2N instead of N
- Block-diagonal structure required
- More complex eigenvalue extraction

**Use Case**: Resonances, exotic hadrons, coupling effects

**JSON Input**:
```json
{
  "basis": {
    "type": "CRG",
    "nmax": 12,
    "rmax": 30.0,
    "rmin": 0.1,
    "omega": 0.5
  }
}
```

---

### CSM (Complex Scaling Method) - ADVANCED, INFRASTRUCTURE READY

**Physical Form**:
```
Complex exponent rotation:
  ν_n = |ν_n| · e^{iθ}  (for all basis functions n)
  Basis: r^l · exp(-ν_n · r²)  where ν_n is complex
```

**Advantages**:
- ✅ Direct complex scaling method (standard CSM theory)
- ✅ Complex eigenvalues naturally represent resonance poles
- ✅ Smaller matrix (N vs 2N)
- ✅ Complex infrastructure fully implemented

**Disadvantages**:
- ✗ Requires complex arithmetic throughout
- ✗ Harder to debug (complex eigenvalues)
- ✗ More literature on CRG than CSM in hadronic physics

**Use Case**: Resonances via complex pole extraction (Gamow states)

**JSON Input**:
```json
{
  "basis": {
    "type": "CSM",
    "nmax": 12,
    "rmax": 30.0,
    "rmin": 0.1,
    "theta": 0.2
  }
}
```

---

## Implementation Priority

### Phase 1: CRG (Hiyama's Method) - **5-8 hours**
**Justification**: Simpler numerics, proven method, no complex arithmetic

1. Add Hiyama basis functions (40 lines)
2. Add paired integration wrapper (30 lines)
3. Implement block matrix construction (80 lines)
4. Test with simple inputs

**Success Criteria**: 
- Eigenvalues with omega=0 match GEM exactly
- Eigenvalues with omega>0 show oscillatory effects
- No NaN or inf values in output

### Phase 2: CSM (Complex Scaling Method) - **4-6 hours**
**Justification**: Reuses complex infrastructure, alternative physics approach

1. Add `getnu_complex()` (5 lines)
2. Add CSM basis functions (30 lines)
3. Dispatch on ORBIT_CSM in spectra.c (50 lines)
4. Test complex eigenvalue extraction

**Success Criteria**:
- Eigenvalues with theta=0 match GEM
- Complex conjugate pairs appear for theta>0
- Resonance pole analysis works

---

## Files to Modify for CRG (Recommended First)

### 1. `include/gemstore/types.h` - Already has correct enum
**Current** (line 24-25):
```c
ORBIT_CRG,
ORBIT_CSM
```
✅ No changes needed - enum values are correct

---

### 2. `include/gemstore/basis/orbit.h` - Add Hiyama declarations

**Add after line 25** (after `getnu()` definition):

```c
/**
 * Hiyama's basis scale parameter (same as standard GEM scale).
 * Used to parameterize cos/sin oscillation frequency.
 */
static inline double getnu_hiyama(int n, int nmax, double rmax, double rmin)
{
    return getnu(n, nmax, rmax, rmin);  // Reuse standard scale calculation
}

/* Hiyama basis function typedefs */
typedef double (*orbit_wfn_hiyama_cos_t)(double x, int n, int l, double nu, double omega);
typedef double (*orbit_wfn_hiyama_sin_t)(double x, int n, int l, double nu, double omega);

/* Forward declarations for Hiyama basis functions */
double HiyamaCos_nlr(double r, int n, int l, double nu, double omega);
double HiyamaSin_nlr(double r, int n, int l, double nu, double omega);
double complex HiyamaCos_nlp(double p, int n, int l, double nu, double omega);
double complex HiyamaSin_nlp(double p, int n, int l, double nu, double omega);
```

---

### 3. `src/basis/orbit.c` - Implement Hiyama basis functions

**Add after existing SHO functions** (after line 48):

```c
/**
 * Hiyama basis function: Cosine component in coordinate space
 * Form: r^l · exp(-ν·r²) · cos(ω·ν·r²)
 * 
 * @param r     Coordinate (fm, real positive)
 * @param n     Basis function index
 * @param l     Orbital angular momentum
 * @param nu    Scale parameter (real, from getnu)
 * @param omega Oscillation frequency parameter
 * @return      Real-valued basis function
 */
double HiyamaCos_nlr(double r, int n, int l, double nu, double omega)
{
    (void)n;  // Radial index not used in base form
    
    // Prefactor: same as GEM
    double prefac = pow(2.0, l/2.0 + 1.25) / sqrt(tgamma(l + 1.5));
    
    // Radial decay: exp(-ν·r²)
    double decay = exp(-nu * r * r);
    
    // Radial power: r^l
    double r_power = pow(r, (double)l);
    
    // Oscillation: cos(ω·ν·r²)
    double omega_nu_r2 = omega * nu * r * r;
    double oscillation = cos(omega_nu_r2);
    
    return prefac * r_power * decay * oscillation;
}

/**
 * Hiyama basis function: Sine component in coordinate space
 * Form: r^l · exp(-ν·r²) · sin(ω·ν·r²)
 */
double HiyamaSin_nlr(double r, int n, int l, double nu, double omega)
{
    (void)n;
    
    double prefac = pow(2.0, l/2.0 + 1.25) / sqrt(tgamma(l + 1.5));
    double decay = exp(-nu * r * r);
    double r_power = pow(r, (double)l);
    double omega_nu_r2 = omega * nu * r * r;
    double oscillation = sin(omega_nu_r2);
    
    return prefac * r_power * decay * oscillation;
}

/**
 * Hiyama basis function: Cosine component in momentum space
 * Form: (-i)^l · [prefactor] · p^l · cos(...) · exp(-p²/(4ν))
 */
double complex HiyamaCos_nlp(double p, int n, int l, double nu, double omega)
{
    (void)n;
    
    double complex phase = cpow(-I, (double)l);
    double prefac = pow(2.0, -l/2.0 - 1.25) / sqrt(tgamma(l + 1.5));
    double p_power = pow(p, (double)l);
    double p2_4nu = p * p / (4.0 * nu);
    double decay = exp(-p2_4nu);
    
    // Oscillation in momentum space (related to coordinate space via Fourier)
    double oscillation = cos(omega * p2_4nu);
    
    return phase * prefac * p_power * decay * oscillation;
}

/**
 * Hiyama basis function: Sine component in momentum space
 */
double complex HiyamaSin_nlp(double p, int n, int l, double nu, double omega)
{
    (void)n;
    
    double complex phase = cpow(-I, (double)l);
    double prefac = pow(2.0, -l/2.0 - 1.25) / sqrt(tgamma(l + 1.5));
    double p_power = pow(p, (double)l);
    double p2_4nu = p * p / (4.0 * nu);
    double decay = exp(-p2_4nu);
    double oscillation = sin(omega * p2_4nu);
    
    return phase * prefac * p_power * decay * oscillation;
}
```

---

### 4. `src/math/integral.c` - Add Hiyama integration wrapper

**Add after line 138** (after existing integrals):

```c
/**
 * Compute matrix element using Hiyama's paired basis functions.
 * For paired cos/sin basis: ⟨ψ_cos | V | ψ_sin ⟩, etc.
 * 
 * Integration handles all four combinations:
 *   cos-cos, cos-sin, sin-cos, sin-sin
 */
double integral_matrix_element_hiyama(
    orbit_wfn_hiyama_cos_t wfn_cos, 
    orbit_wfn_hiyama_sin_t wfn_sin,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra, 
    const argsOrbit_t *args_ket,
    const argsGIModel_t *args_model, 
    const argsDynmc_t *args_dynmc,
    int bra_type, int ket_type)  // 0=cos, 1=sin
{
    double sum = 0.0;
    double r_node, bra_val, ket_val, pot_val;
    int i;
    
    for (i = 0; i < QUAD_ORDER; i++) {
        r_node = node_factor * nodes[i];
        
        // Select basis function for bra
        if (bra_type == 0) {
            bra_val = wfn_cos(r_node, args_bra->n, args_bra->l, 
                             args_bra->scale, 0.0);  // omega from args
        } else {
            bra_val = wfn_sin(r_node, args_bra->n, args_bra->l, 
                             args_bra->scale, 0.0);
        }
        
        // Select basis function for ket
        if (ket_type == 0) {
            ket_val = wfn_cos(r_node, args_ket->n, args_ket->l, 
                             args_ket->scale, 0.0);
        } else {
            ket_val = wfn_sin(r_node, args_ket->n, args_ket->l, 
                             args_ket->scale, 0.0);
        }
        
        pot_val = pot(r_node, args_model, args_dynmc);
        sum += weights[i] * bra_val * pot_val * ket_val;
    }
    return sum;
}
```

---

### 5. `src/model/spectra.c` - Add CRG dispatch

**Modify lines 44-49** (basis setup):

```c
// BEFORE:
for (i = 0; i < nmax; i++) {
    nu_val[i] = getnu(i + 1, nmax, rmax, rmin);
    args_orb_bra.scale = nu_val[i];
    args_orb_bra.n = i + 1;
}

// AFTER - add dispatch:
if (input->orbit == ORBIT_CRG) {
    double *nu_hiyama = (double *)malloc(nmax * sizeof(double));
    for (i = 0; i < nmax; i++) {
        nu_hiyama[i] = getnu_hiyama(i + 1, nmax, rmax, rmin);
    }
} else {
    for (i = 0; i < nmax; i++) {
        nu_val[i] = getnu(i + 1, nmax, rmax, rmin);
    }
}
```

**Add CRG block matrix construction** (after line 204, before eigenvalue solve):

```c
if (input->orbit == ORBIT_CRG) {
    // Hiyama method: build 2Nx2N block matrices for paired cos/sin basis
    int nmatrix = 2 * nmax;  // cos and sin for each n
    matrix_t mT_block, mH_block, mN_block;
    
    mT_block = matrix_init(nmatrix, nmatrix);
    mH_block = matrix_init(nmatrix, nmatrix);
    mN_block = matrix_init(nmatrix, nmatrix);
    
    // Fill 2×2 blocks for each (i,j) pair
    for (i = 0; i < nmax; i++) {
        for (j = i; j < nmax; j++) {
            args_orb_bra.n = i + 1;
            args_orb_bra.scale = nu_hiyama[i];
            args_orb_ket.n = j + 1;
            args_orb_ket.scale = nu_hiyama[j];
            
            // Compute four matrix elements for 2×2 block
            double T_cc = integral_matrix_element_hiyama(HiyamaCos_nlr, HiyamaSin_nlr,
                GIVt, node_factor, &args_orb_bra, &args_orb_ket, &args_model, &args_dynmc, 0, 0);
            double T_cs = integral_matrix_element_hiyama(HiyamaCos_nlr, HiyamaSin_nlr,
                GIVt, node_factor, &args_orb_bra, &args_orb_ket, &args_model, &args_dynmc, 0, 1);
            double T_sc = integral_matrix_element_hiyama(HiyamaCos_nlr, HiyamaSin_nlr,
                GIVt, node_factor, &args_orb_bra, &args_orb_ket, &args_model, &args_dynmc, 1, 0);
            double T_ss = integral_matrix_element_hiyama(HiyamaCos_nlr, HiyamaSin_nlr,
                GIVt, node_factor, &args_orb_bra, &args_orb_ket, &args_model, &args_dynmc, 1, 1);
            
            // Place in block structure (indices account for cos/sin pairing)
            mT_block.value[2*i][2*j] = T_cc;
            mT_block.value[2*i][2*j+1] = T_cs;
            mT_block.value[2*i+1][2*j] = T_sc;
            mT_block.value[2*i+1][2*j+1] = T_ss;
            
            // Symmetric placement
            if (i != j) {
                mT_block.value[2*j][2*i] = T_cc;
                mT_block.value[2*j+1][2*i] = T_cs;
                mT_block.value[2*j][2*i+1] = T_sc;
                mT_block.value[2*j+1][2*i+1] = T_ss;
            }
        }
    }
    
    // Repeat for all potentials (V_coul, V_conf, etc.)
    // ... Similar pattern for mH_block
    
    // Solve generalized eigenvalue: H·c = λ·N·c
    eigen_general(mH_block.value, mN_block.value, nmatrix, eigenvalues, eigenvectors, nmatrix);
    
} else {
    // Existing GEM/SHO code (unchanged)
}
```

---

### 6. `src/entry.c` - Add CRG validation

**Modify line 25-34** (validation):

```c
// BEFORE:
if (input.orbit != ORBIT_GEM)
    exit_error("SPECTRA meson only supports GEM basis");

// AFTER:
if (input.orbit != ORBIT_GEM && input.orbit != ORBIT_CRG)
    exit_error("SPECTRA meson supports GEM and CRG (Hiyama) basis only");

if (input.orbit == ORBIT_CRG && input.omega == 0.0)
    fprintf(stderr, "Warning: CRG with omega=0 is equivalent to GEM\n");
```

---

## Implementation Timeline for CRG

| Phase | Task | Time | Status |
|-------|------|------|--------|
| 1 | Add Hiyama declarations to orbit.h | 10 min | ⏳ Pending |
| 2 | Implement Hiyama basis functions in orbit.c | 30 min | ⏳ Pending |
| 3 | Add integration wrapper in integral.c | 20 min | ⏳ Pending |
| 4 | Add CRG dispatch in spectra.c | 1 hour | ⏳ Pending |
| 5 | Add validation in entry.c | 10 min | ⏳ Pending |
| 6 | Test and debug | 2-3 hours | ⏳ Pending |
| **Total** | **CRG Implementation** | **~5 hours** | ⏳ Pending |

---

## Test Cases for CRG

### Test 1: omega=0 baseline
```json
{
  "basis": {"type": "CRG", "nmax": 12, "omega": 0.0}
}
```
**Expected**: Output matches GEM exactly

### Test 2: Small omega
```json
{
  "basis": {"type": "CRG", "nmax": 12, "omega": 0.1}
}
```
**Expected**: Slight deviation from GEM, no NaN values

### Test 3: Moderate omega
```json
{
  "basis": {"type": "CRG", "nmax": 12, "omega": 0.5}
}
```
**Expected**: Significant oscillatory effects in basis

### Test 4: Matrix block structure
**Check**: Verify 2N×2N matrix is built correctly with block diagonal structure

---

## Next: CSM Implementation (Advanced)

Once CRG is working, CSM can be implemented similarly:

1. Add `getnu_complex()` to orbit.h (5 lines)
2. Add CSM basis functions to orbit.c (30 lines)
3. Add complex matrix dispatch in spectra.c (50 lines)
4. Reuse existing complex eigenvalue solver

CSM adds real and imaginary parts directly to eigenvalues (resonance poles).

---

## References & Resources

**Key Implementation Files**:
- `include/gemstore/basis/orbit.h` (63 lines)
- `src/basis/orbit.c` (49 lines → +100 lines after CRG/CSM)
- `src/math/integral.c` (138 lines → +50 lines after CRG)
- `src/model/spectra.c` (234 lines → +150 lines after CRG)

**Physics References** (Hiyama's Method):
- Hiyama, E., Kamimura, M. (1995). "Gaussian Expansion Method for few-body problems"
- Applications in: Meson spectroscopy, baryon resonances, exotic hadrons

**Example Input Files**:
- `app/amethyst.json` (GEM example)
- Will create: `app/charmonium_crg.json` (CRG example)

**Git History**:
```bash
git log --oneline | grep -i "crg\|hiyama\|complex"
# faa7539 a cross road for CSM and CRG
```

---

## Summary

**CRG (Hiyama's Complex-Range Gaussian Method)**:
- ✅ Simpler numerics (purely real)
- ✅ Proven in hadronic physics
- ✅ 2N basis functions (cos/sin pairs)
- ✅ ~5 hours to implement
- ✅ **RECOMMENDED for near-term**

**CSM (Complex Scaling Method)**:
- ✅ Standard resonance method
- ✅ Complex eigenvalues = poles
- ✅ N basis functions
- ✅ ~4 hours to implement (infrastructure ready)
- ⏳ **For advanced users after CRG**

**Status**: Infrastructure 100% ready. Just need to wire basis functions together!
