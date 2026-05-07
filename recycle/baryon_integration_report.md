# Baryon GEM Integration Report

## Scope

This report compares the active meson implementation in `src/` with the legacy baryon implementation in `recycle/`, and proposes a concrete path for adding baryon computation into the current codebase without importing the old design wholesale.

## Code Paths Reviewed

### Active codebase

- Entry and dispatch: `src/main.c`, `src/entry.c`
- Input parsing: `src/parse.c`
- Meson compute dispatcher: `src/model/compute.c`
- Meson GEM and CRG solvers: `src/model/cmeson.c`
- Orbital basis helpers: `src/basis/orbit.c`, `include/gemstore/basis/orbit.h`
- GI potential interfaces: `src/model/gimodel.c`, `include/gemstore/model/gimodel.h`
- Shared argument types: `include/gemstore/param/argset.h`
- Output: `src/print.c`

### Legacy baryon code in `recycle/`

- Top-level baryon task: `recycle/task.c`
- Matrix-element assembly and basis generation: `recycle/mfi.c`, `recycle/mfi.h`
- Pairwise operator kernels: `recycle/vtype.c`, `recycle/vtype.h`
- Auxiliary integral/eigen helpers: `recycle/inteCenV.*`, `recycle/sumckdk.*`, `recycle/eigen.*`

## Main Differences Between Implementations

### 1. System architecture

The active codebase is organized around a small runtime pipeline:

1. Parse JSON input into `argsInput_t`
2. Dispatch by task and system in `src/entry.c`
3. Dispatch by basis in `src/model/compute.c`
4. Solve spectra in system-specific solver code
5. Write structured output in `src/print.c`

The legacy baryon implementation is structured as a self-contained procedural task:

1. Hardcode model parameters and basis settings in `task_baryon()`
2. Generate baryon quantum-number lists
3. Precompute many pairwise operator caches
4. Build `Nfi` and `Hfi`
5. Solve the generalized eigenproblem
6. Print a short numerical result

The active meson path is integrated into the application. The baryon path is currently a standalone prototype.

### 2. Input model

The active code accepts external JSON input via `parse_input_file()` in `src/parse.c`.

- It supports only `system.type = "MESON"`
- It parses basis configuration and model parameter set selection
- It stores the result in `argsInput_t`

The legacy baryon code does not read application input.

- Quark masses, potential parameters, `nmax`, `rmin`, and `rmax` are hardcoded in `recycle/task.c`
- Baryon quantum numbers are passed as positional arguments to `task_baryon()`
- There is no JSON schema, no parameter-set dispatch, and no output serialization

### 3. Basis representation

The active meson implementation uses a simple radial basis representation:

- `argsOrbit_t` with `n`, `l`, `scale`, `param`
- One orbital channel per run
- Basis size is directly `nmax`

The legacy baryon implementation uses a much richer coupled basis:

- Jacobi-coordinate quantum numbers for `rho` and `lambda`
- Pair index `c = 1, 2, 3`
- Coupled spin/orbital channels with permutation structure
- Two layered basis lists: `qnlist_spfy` and `qnlist_full`
- Radial expansion added after the spin-orbital basis is enumerated

This is the biggest structural gap. Baryons are not a minor extension of the meson radial basis. They require a different basis generator and different matrix-element assembly logic.

### 4. Hamiltonian construction

The active meson solver in `src/model/cmeson.c` builds Hamiltonian pieces directly as dense `nmax x nmax` matrices:

- kinetic term
- Coulomb and confinement terms
- contact, spin-orbit, Thomas, tensor terms
- overlap matrix

Each element is computed from generic integral functions with the current basis and GI model.

The baryon code uses a two-stage assembly:

1. Precompute symbolic-angular and Gaussian-combinatoric objects in `sumckdk_scdk`
2. Reuse those cached objects in `getmfi()` to accumulate matrix elements for each basis pair

This cache-heavy structure exists because baryon matrix elements are much more coupled and expensive than the meson radial case.

### 5. Reuse of GI model code

The current GI potential layer is partly prepared for baryons.

- `argsGIModelDy_t` already has `system`
- `src/model/gimodel.c` already contains `SYSTEM_BARYON` sign branches for some spin-orbit terms

However, the active runtime does not expose baryon execution.

- `src/entry.c` rejects all non-meson spectra runs
- `src/parse.c` rejects all non-meson systems
- `src/model/compute.c` only calls meson solvers

So low-level model support is incomplete but not absent.

### 6. Linear algebra strategy

The meson path orthogonalizes the overlap matrix with Cholesky and then solves a generalized eigenproblem in `src/model/cmeson.c`.

The baryon path directly forms `Nfi` and `Hfi` and solves using `eigv2Mul()` in `recycle/task.c`.

Functionally these are compatible at a high level, but the active code already has a cleaner matrix abstraction and should remain the target architecture.

### 7. Output behavior

The active meson path produces:

- console summaries
- JSON output via `write_meson_spectra()`

The legacy baryon path only prints a small amount of console output. It has no structured output model and no integration with the current project/reporting pipeline.

### 8. Code quality and maintainability

The active code is not perfect, but it is much easier to extend safely because it already follows the application structure.

The `recycle/` baryon implementation has several traits that should not be copied directly:

- hardcoded physics parameters inside compute code
- custom ad hoc memory graphs with four-level pointer trees
- weak separation between basis generation, physics setup, and runtime orchestration
- no parser or output integration
- legacy includes and naming that do not match current `gemstore/...` headers

The right approach is to port the baryon physics and basis logic, not the old application structure.

## What Already Exists in the Current Codebase

The following pieces reduce the amount of new work needed:

- `SYSTEM_BARYON` enum already exists in `include/gemstore/types.h`
- `system_type_str[]` already includes `"BARYON"` in `src/types.c`
- `argsInput_t` already contains `f1`, `f2`, `f3`, `S`, `L`, `jl`, `J`, `nmax`, `rmin`, `rmax`
- `argsGIModelDy_t` already carries `system`, pair masses, color factor, and operator values
- `src/model/gimodel.c` already distinguishes meson and baryon sign conventions in part of the GI spin-dependent potential code

This means the codebase already expects baryon support conceptually, but the main runtime and solver layers are still meson-only.

## Missing Pieces for Baryon Support

### 1. Input parsing and validation

`src/parse.c` must be extended to accept baryon systems.

At minimum, baryon parsing needs:

- `system.type = "BARYON"`
- `f1`, `f2`, `f3`
- baryon spin-coupling quantum numbers, likely including `S`, `L`, `J`
- one baryon-exchange symmetry selector corresponding to legacy `f12`
- optionally `Lmax` if the baryon basis is generated from a truncation shell instead of a single fixed `L`

The current `argsInput_t` does not have fields for `f12`, parity `P`, or `Lmax`, so those will need to be added if the recycled basis-generation logic is kept in recognizable form.

### 2. Runtime dispatch

`src/entry.c` and `src/model/compute.c` must stop assuming meson-only spectra.

Needed changes:

- allow `TASK_SPECTRA` with `SYSTEM_BARYON`
- add `compute_spectra_baryon()`
- keep basis validation system-specific

### 3. Baryon basis layer

The current active basis layer is radial and meson-oriented. Baryons need a dedicated basis module, likely under `src/basis/` and `include/gemstore/basis/`.

Recommended new modules:

- `basis/baryon.h`
- `basis/baryon.c`

Responsibilities:

- enumerate baryon spin-orbital channels
- generate Jacobi-coordinate Gaussian basis states
- represent channel metadata cleanly
- separate symbolic channel generation from radial expansion

The old `basis_list` machinery in `recycle/` is a useful reference, but it should be redesigned with current naming and ownership rules.

### 4. Baryon matrix-element layer

The old `sumckdk`, `vtype`, and `inteCenV` code contains the baryon-specific physics machinery. That machinery needs to be migrated into dedicated source files under `src/model/` or `src/math/`, not called directly from `recycle/`.

Recommended split:

- angular and spin recoupling helpers
- cached pairwise kernel generation
- radial/integral evaluation
- full `Nfi` and `Hfi` assembly

This port should hide old implementation details behind smaller interfaces instead of exposing `sumckdk_scdk ****` throughout the new solver.

### 5. Baryon solver

Create a new solver module parallel to `src/model/cmeson.c`.

Recommended files:

- `include/gemstore/model/cbaryon.h`
- `src/model/cbaryon.c`

Responsibilities:

- build baryon basis from parsed input
- evaluate overlap and Hamiltonian matrices
- solve the generalized eigenproblem with the current matrix/eigen layer
- optionally compute baryon observables such as RMS radii later

This solver should follow the current application style rather than the `task_baryon()` style.

### 6. Output writer

The project needs a baryon output path similar to `write_meson_spectra()`.

Recommended additions:

- `print_baryon_spectra()`
- `write_baryon_spectra()`

The JSON should include:

- baryon flavor content `f1`, `f2`, `f3`
- basis truncation settings
- quantum numbers used for the run
- eigenvalues and eigenvectors
- optional channel decomposition metadata if available

## Recommended Integration Strategy

### Guiding principle

Do not integrate baryons by exposing `recycle/task_baryon()` as a special-case runtime path.

That would create a second application architecture inside the repository. Instead, port the baryon physics into the existing `parse -> dispatch -> compute -> print` pipeline.

### Phase 1: Make the runtime baryon-aware

1. Extend `argsInput_t` with baryon-specific fields that are truly required.
2. Update `src/parse.c` to accept `system.type = "BARYON"`.
3. Update `src/entry.c` to dispatch baryon spectra.
4. Add a stub `compute_spectra_baryon()` that validates input and reports unsupported basis combinations clearly.

This phase should introduce no baryon physics yet. It just establishes the public integration points.

### Phase 2: Port baryon basis generation

1. Extract the logical content of `getlsj_sl()` and `basis_mlsj_sl()` from `recycle/mfi.c`.
2. Redesign the data structures for active code conventions.
3. Add clean ownership and cleanup helpers.
4. Add debug printing for generated baryon channels.

Goal: produce a clean in-memory baryon basis object that the solver can consume.

### Phase 3: Port matrix-element assembly

1. Move the useful `vtype`, `sumckdk`, and integral logic into namespaced `gemstore` modules.
2. Replace legacy headers and inconsistent types.
3. Keep the cache strategy where it materially reduces cost.
4. Wrap the cache implementation behind a smaller interface used by `cbaryon.c`.

Goal: obtain `Nfi` and `Hfi` for baryon channels using the current matrix layer.

### Phase 4: Add baryon solver and output

1. Implement `spectra_baryon_GEM()` or a similarly named baryon solve function.
2. Add console and JSON output.
3. Add at least one example baryon input file.
4. Validate the spectrum against trusted reference numbers from the old code.

### Phase 5: Cleanup and convergence

1. Remove duplicated or dead baryon support code that becomes obsolete.
2. Keep `recycle/` only as archival reference if still useful.
3. Document the baryon input schema in `README.md`.

## Recommended File-Level Additions

A minimal clean integration would likely add or modify the following files.

### Modify

- `include/gemstore/param/argset.h`
- `src/parse.c`
- `src/entry.c`
- `src/model/compute.c`
- `src/print.c`
- `include/gemstore/print.h`

### Add

- `include/gemstore/model/cbaryon.h`
- `src/model/cbaryon.c`
- `include/gemstore/basis/baryon.h`
- `src/basis/baryon.c`
- one or more baryon-specific helper modules for matrix-element assembly

## Proposed Input Extensions

The baryon parser should avoid overloading meson semantics where they do not fit. A practical JSON shape would look like this:

```json
{
  "project": "omega_sss_ground",
  "task": "SPECTRA",
  "model": {
    "type": "GISCREEN",
    "param": "GISCREEN_CUSTOM",
    "file": "params/baryon_giscreen.json"
  },
  "system": {
    "type": "BARYON",
    "f1": 2,
    "f2": 2,
    "f3": 2,
    "sym12": 1,
    "J": 1.5,
    "P": 1,
    "Lmax": 2
  },
  "basis": {
    "type": "GEM",
    "nmax": 6,
    "rmin": 0.2,
    "rmax": 2.0
  }
}
```

The exact quantum-number schema may change, but it should reflect the real baryon basis generator rather than pretending baryons use the meson `(S, L, J)` input shape unchanged.

## Specific Porting Advice From `recycle/`

### Reuse conceptually

- `getlsj_sl()` basis enumeration logic
- `basis_mlsj_sl()` spin-orbital coefficient generation
- pair remapping helpers in `vtype.c`
- pairwise operator kernels `vcent`, `vcont`, `vtens`, `vsorp`, `vsorn`
- matrix-element accumulation pattern in `getmfi()`

### Do not reuse verbatim

- `task_baryon()` as the runtime entry point
- hardcoded parameters in `task.c`
- raw multi-level allocation patterns exposed across modules
- legacy non-namespaced headers
- minimal console-only output behavior

## Risks and Design Constraints

### 1. The current basis header is misleading

`include/gemstore/basis/basis.h` contains structures named `basis_meson_t` and `basis_baryon_t`, but `basis_meson_t` already includes many baryon-style fields such as `m3`, `lrho`, `llam`, `nrho`, and `nlam`.

This suggests leftover design drift from the older code. Before baryon integration, this header should be treated carefully and likely redesigned instead of extended blindly.

### 2. Baryon support is not just a new switch case

Adding `SYSTEM_BARYON` to parser and dispatch is easy. The difficult part is building a maintainable baryon basis and Hamiltonian assembly layer that matches the active architecture.

### 3. Parameter-model semantics may diverge

The old baryon code uses a different parameter set style from the active meson GI code. The project must decide whether baryons will:

- reuse `argsGIModel_t` directly,
- extend it for baryon-specific parameters,
- or introduce a baryon model struct derived from a shared base.

The cleanest near-term option is to extend `argsGIModel_t` only if required by actual baryon physics terms.

## Recommended First Implementation Milestone

The best first deliverable is not full baryon physics. It is:

1. parse baryon JSON input,
2. dispatch to `compute_spectra_baryon()`,
3. build a baryon basis object for a small test case,
4. print a basis summary,
5. stop before Hamiltonian solve.

That milestone forces the codebase to settle the input schema and basis representation early, which are the hardest API decisions. After that, the matrix-element port becomes much safer.

## Summary

The legacy baryon implementation in `recycle/` contains valuable physics and basis-generation logic, but it does not fit the architecture of the active codebase as-is.

The correct path is:

1. extend the active runtime to recognize baryons,
2. port the baryon basis generator into new active modules,
3. port the baryon matrix-element machinery behind cleaner interfaces,
4. add a dedicated baryon solver and output writer,
5. validate against the old implementation.

In short, reuse the baryon mathematics from `recycle/`, but do not reuse the `recycle/` application structure.
