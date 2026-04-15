<div align="center">

![GEMSTORE Logo](doc/logo.png)

# GEMSTORE: Hadron Spectroscopy Simulation Tools

[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0-or-later)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Platform: Linux/Unix](https://img.shields.io/badge/Platform-Linux%2FUnix-lightgrey.svg)](https://en.wikipedia.org/wiki/Unix)
![Status: Active Development](https://img.shields.io/badge/Status-Active%20Development-brightgreen.svg)

**A powerful computational framework for calculating hadron spectroscopy using the Gaussian Expanding Method (GEM) and Godfrey-Isgur quark models.**

</div>

---

## Overview

**GEMSTORE** is a specialized scientific software suite for hadron spectroscopy calculations using advanced quark models. It combines the computational efficiency of the **Gaussian Expanding Method** with the physics-rich **screen-modified Godfrey-Isgur (GI)** model to predict meson, baryon, and exotic hadron properties.

### Key Features

- 🎯 **Multiple Quark Models**: GI-Screen, GI-String, GI-Quadratic, MIT Bag Model
- 📊 **Comprehensive Spectral Calculations**: Masses, radii, decay widths, coupling constants
- 🔧 **Flexible Quantum Numbers**: Full support for arbitrary L, S, J combinations
- 🤖 **AI-Assisted Workflows**: Integrated OpenCode assistant for intelligent task automation
- 📈 **Advanced Data Analysis**: Eigenvector analysis, normalization validation, statistical summaries
- ⚡ **High Performance**: Optimized C implementation with Minuit2 numerical library
- 🔬 **Parameter Fitting**: Gradient-based optimization with Minuit2

---

## Table of Contents

1. [Features](#features)
2. [System Requirements](#system-requirements)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Core Algorithms](#core-algorithms)
6. [AI Integration](#ai-integration)
7. [Main Functions](#main-functions)
8. [Usage Examples](#usage-examples)
9. [Project Structure](#project-structure)
10. [Contributing](#contributing)
11. [License](#license)

---

## Features

### Physical Models

| Model | Description | Use Case |
|-------|-------------|----------|
| **GI-Screen** | Screen-modified Godfrey-Isgur potential with Coulomb screening | Heavy quarkonium (charmonium, bottomonium) |
| **GI-String** | String-like linear confinement | Light mesons and general meson spectra |
| **GI-Quadratic** | Quadratic confinement potential | Theoretical studies, precision fits |
| **MIT Bag** | Bag model confinement | Baryon spectroscopy, validation studies |

### Calculation Types

- **SPECTRA**: Calculate complete hadron mass spectra
- **RADIUS**: Compute RMS radii and spatial distributions
- **DECAY3P0**: Calculate 3P₀ OZI-allowed decay widths
- **COUPLCHN**: Compute coupling constants to meson channels
- **SCATTER**: Scattering amplitudes and cross sections

### System Types

- **MESON**: Quark-antiquark bound states
- **BARYON**: Three-quark systems
- **MOLECULE**: Multi-hadron clusters and exotics

---

## System Requirements

### Minimum Requirements

- **OS**: Linux/Unix (macOS with GNU tools)
- **Compiler**: GCC 7.0+ or Clang 5.0+
- **Build Tool**: GNU Make 4.0+
- **Memory**: 512 MB RAM
- **Storage**: 100 MB installation + 1 GB for calculations

### Optional Dependencies

- **LAPACK/OpenBLAS**: For accelerated linear algebra
- **Python 3.7+**: For input file generation scripts
- **OpenCode**: For AI-assisted workflow integration

---

## Installation

### Clone the Repository

```bash
git clone https://github.com/serialcore/gemstore.git
cd gemstore
```

### Build from Source

```bash
# Standard build
make clean && make

# With LAPACK acceleration
make clean && make USE_LAPACKE=1

# Install to system (requires sudo)
make install

# Uninstall from system
make uninstall
```

### Verify Installation

```bash
./gemstore --help
./gemstore --version
```

---

## Quick Start

### Basic Usage

```bash
# Run a charmonium calculation
./gemstore --input app/amethyst.inp

# Print the confinement potential
./gemstore --print potential --input app/amethyst.inp

# Debug spin-orbit coupling operator
./gemstore --debug soc_operator --input app/amethyst.inp
```

### Create Input File

Create `my_meson.inp`:

```ini
&GLOBAL
  project = my_project
  task = SPECTRA              # SPECTRA | RADIUS | DECAY3P0
&END

&SYSTEM
  model = GI_SCREEN           # GI_SCREEN | GI_STRING | GI_QUADRA
  system = MESON              # MESON | BARYON
&END

&PARAMS
  params = GIScreen_ccbar      # Predefined parameter set
&END

&QUANTUM
  f1 = 3            # Quark flavor 1 (3=charm)
  f2 = 3            # Quark flavor 2
  S = 1             # Spin quantum number
  L = 0             # Orbital angular momentum
  J = 1             # Total angular momentum
&END

&GAUSS
  nmax = 16         # Number of Gaussian basis functions
  rmax = 30.0       # Maximum radius (fm)
  rmin = 0.1        # Minimum radius (fm)
&END
```

Run it:

```bash
./gemstore --input my_meson.inp
```

---

## Core Algorithms

### 1. Gaussian Expanding Method (GEM)

The radial wavefunction is expanded in Gaussian basis functions:

```
ψ(r) = Σ cₙ φₙ(r)
```

where each basis function is:

```
φₙ(r) = r^ℓ exp(-αₙ r²)
```

**Advantages**:
- All integrals solvable in closed form
- Exponential convergence with basis size
- Efficient matrix computations

**Implementation**: `src/basis/orbit.c`, `src/math/integral.c`

### 2. Godfrey-Isgur Quark Model

The GI potential combines three components:

**A. Confinement Potential**
```
V_conf(r) = b₁ · r + b₂ + const
```

**B. Coulomb Interaction**
```
V_coul(r) = -αₛ(r) · Cᶠ / r
```

**C. Hyperfine Interactions**
- Spin-Spin: Contact term
- Spin-Orbit: L·S coupling
- Tensor: Tensor operator

**Implementation**: `src/model/gimodel.c` (20+ potential functions)

### 3. Eigenvalue Problem Solution

Converts the radial Schrödinger equation into a generalized eigenvalue problem:

```
H c = E S c
```

**Algorithm**:
1. Generate Gaussian basis set
2. Compute overlap matrix S (analytically)
3. Compute kinetic energy matrix T (analytically)
4. Compute potential energy matrix V (numerical integration)
5. Solve generalized eigenvalue problem
6. Extract masses and wavefunctions

**Implementation**: `src/math/eigen.c`, `src/model/spectra.c`

### 4. Spectral Data Post-Processing

The RMS radius spectrum is analyzed for anomalies using quadratic interpolation:

1. **Detect anomalies**: Identify non-monotonic points
2. **Interpolate**: Use Lagrange quadratic/linear interpolation
3. **Validate**: Check normalization and consistency

**Implementation**: `src/math/interplt.c` (enhanced with debugging)

### 5. Parameter Fitting Engine

Uses **Minuit2** numerical optimization library:
- MIGRAD: Gradient-based minimization
- SIMPLEX: Gradient-free optimization
- Parameter covariance matrix estimation

**Fit function**:
```
χ² = Σᵢ (M_calc^i - M_exp^i)² / σᵢ²
```

**Implementation**: `src/param/fitting.c`, `src/param/meson.cc`

---

## AI Integration

### OpenCode Assistant for GEMSTORE

The `app/gemstore-assistant/` subdirectory contains an **OpenCode skill** for intelligent task automation.

#### AI Capabilities

```
"Calculate charmonium 1P state with J=1"
         ↓
Parse quantum numbers, select model
         ↓
Generate input file automatically
         ↓
Execute gemstore calculation
         ↓
Parse results, format output
         ↓
"ψ(J^PC) = 1^-- with M = 3.686 GeV"
```

#### Features

- **Natural Language Understanding**: Parse physics requests
- **Automated Workflow**: Generate inputs, execute, parse outputs
- **Systematic Calculations**: Run multiple L, S, J combinations
- **Result Interpretation**: Physics-meaningful explanations
- **Parameter Fitting**: Intelligent model selection

#### Skill Location

```
app/gemstore-assistant/
├── SKILL.md                          # Main skill definition
├── scripts/generate_meson_inputs.py  # Auto-generate input files
└── templates/meson_spectra_template.md
```

#### Usage with OpenCode

```bash
# Enable AI-assisted calculations
opencode "Calculate charmonium spectrum up to L=2"

# Interactive model selection
opencode "Fit GIScreen parameters to experimental data"
```

For details: see `app/gemstore-assistant/SKILL.md`

---

## Main Functions

### Entry Points (src/entry.c)

| Function | Purpose | Call |
|----------|---------|------|
| `entry_compute()` | Run spectroscopy calculation | `--input <file>` |
| `entry_fitting()` | Parameter optimization | `--fitting <target>` |
| `entry_debug()` | Debug calculation steps | `--debug <unit>` |
| `entry_print()` | Print potential/wavefunction | `--print <item>` |

### Core Spectroscopy (src/model/)

| Function | Algorithm | Output |
|----------|-----------|--------|
| `spectra_meson_GI()` | Solve Schrödinger equation | Eigenvalues (masses) |
| `radius_meson_rms()` | Compute ⟨r²⟩^(1/2) | RMS radii |
| `interpolate_quadratic()` | Fix anomalies in spectra | Corrected data |
| `write_meson_spectra()` | Generate report | `.out` file |

### Potential Functions (src/model/gimodel.c)

```c
double GIVconf(double r, ...)     // Confinement
double GIVcoul(double r, ...)     // Coulomb
double GIVcont(double r, ...)     // Contact (spin-spin)
double GIVsovij(double r, ...)    // Spin-orbit
double GIVtens(double r, ...)     // Tensor force
```

### Basis Functions (src/basis/)

| Module | Purpose |
|--------|---------|
| `orbit.c` | Orbital angular momentum basis |
| `spin.c` | Spin SU(2) Clebsch-Gordan coefficients |
| `color.c` | SU(3) color factors |
| `isospin.c` | Isospin basis states |
| `intrin.c` | Intrinsic wavefunction representation |

### Mathematics (src/math/)

| Module | Algorithms |
|--------|-----------|
| `matrix.c` | Linear algebra |
| `eigen.c` | Generalized eigenvalue solver |
| `integral.c` | Gaussian quadrature integration |
| `cmi.c` | Color magnetic interaction |
| `soc.c` | Spin-orbit coupling |
| `su3.c` | SU(3) group operations |
| `interplt.c` | Spectral data interpolation |

---

## Usage Examples

### Example 1: Charmonium Ground State

**Input file** (`cc_ground.inp`):
```ini
&GLOBAL
  project = charmonium_ground
  task = SPECTRA
&END
&SYSTEM
  model = GI_SCREEN
  system = MESON
&END
&PARAMS
  params = GIScreen_ccbar
&END
&QUANTUM
  f1 = 3  f2 = 3        # Both charm quarks
  S = 0   L = 0  J = 0  # S-wave, singlet
&END
&GAUSS
  nmax = 16  rmax = 30.0  rmin = 0.1
&END
```

**Run**:
```bash
./gemstore --input cc_ground.inp
```

**Output excerpt** (`charmonium_ground.out`):
```
State    Mass(GeV)    RMS(fm)      Δmass        max|coeff|   ||coeff||^2   norm_stat
------+----------+----------+----------+---------------+---------------+---------------
1        3.096788    0.524365    -0.001234    0.98765432    1.000000000   ✓ OK
```

### Example 2: Print Potential

```bash
./gemstore --print potential --input app/amethyst.inp
```

### Example 3: Debug Analysis

```bash
./gemstore --debug soc_operator --input app/amethyst.inp
```

---

## Output Files

### Standard Output (.out)

Generated by `write_meson_spectra()`:

```
================================================================================
                      MESON SPECTROSCOPY RESULTS SUMMARY
================================================================================

Generated:   2026-04-15 14:30:45
Project:     charmonium

INPUT CONFIGURATION:
  Quark Flavors:       f1=3  f2=3
  Angular Momentum:    S=1.0  L=0.0  J=1.0
  Gaussian Basis:      nmax=16  rmin=0.1 fm  rmax=30.0 fm

MODEL PARAMETERS:
  Quark Masses:        mn=0.471346  ms=0.628312  mc=1.810505  mb=5.156015 GeV
  Potential:           b1=0.257547  mu=0.145356  c=-0.658943

SPECTRAL DATA:
State    Mass(GeV)    RMS(fm)      Δmass        max|coeff|   ||coeff||^2   norm_stat
...

EIGENVECTOR COMPONENTS:
State    c[0]         c[1]         c[2]  ...
...

STATISTICAL SUMMARY:
  Number of states:    12
  Mass Statistics (GeV):
    Min:               0.770000
    Max:               2.110000
    Mean:              1.234567
    Std Dev:           0.123456

================================================================================
```

---

## Project Structure

```
gemstore/
├── README.md                      # This file
├── LICENSE                        # GPLv3 license
├── Makefile                       # Build configuration
├── gemstore                       # Compiled executable
│
├── include/gemstore/              # Public API headers
│   ├── basis/                     # Wavefunction basis
│   ├── math/                      # Mathematics library
│   ├── model/                     # Physics models
│   ├── param/                     # Input parameters
│   └── entry.h, print.h           # Main interface
│
├── src/                           # Implementation (~5,000 LOC)
│   ├── main.c                     # Entry point
│   ├── entry.c                    # Task routing
│   ├── print.c                    # Output formatting (enhanced)
│   ├── basis/                     # Basis functions (~800 LOC)
│   ├── math/                      # Mathematics (~1200 LOC)
│   ├── model/                     # Physics models (~1800 LOC)
│   └── param/                     # Parameter fitting
│
├── lib/                           # External libraries
│   └── Minuit2/                   # Numerical optimization (CERN)
│
├── app/                           # Applications
│   ├── amethyst.inp               # Example: charmonium
│   └── gemstore-assistant/        # AI Integration
│       ├── SKILL.md               # Skill definition
│       ├── scripts/
│       └── templates/
│
├── doc/                           # Documentation
│   └── logo.png                   # GEMSTORE logo
│
└── test/                          # Test suite

Total LOC: ~7,800 (C + C++ + Headers)
```

---

## Building with LAPACK Acceleration

```bash
# Install LAPACK/OpenBLAS (Ubuntu/Debian)
sudo apt-get install liblapacke-dev libopenblas-dev

# Build with LAPACKE support
make clean && make USE_LAPACKE=1
```

---

## Contributing

### Guidelines

1. **Code Style**: K&R style with 4-space indentation
2. **Documentation**: Add docstrings for all public functions
3. **Testing**: Include unit tests for new algorithms
4. **Physics**: Cite literature for new models
5. **Performance**: Profile before optimizing

### Development Workflow

```bash
# Create feature branch
git checkout -b feature/my-algorithm

# Make changes and test
make clean && make

# Commit with descriptive messages
git add -A
git commit -m "Add [feature]: description"

# Push and open pull request
git push origin feature/my-algorithm
```

---

## Citation

If you use GEMSTORE in research, please cite:

```bibtex
@software{gemstore2026,
  author = {Zhang, Wen-Xuan},
  title = {GEMSTORE: Hadron Spectroscopy Simulation Tools},
  year = {2026},
  url = {https://github.com/serialcore/gemstore},
  version = {0.1.4},
  note = {Gaussian Expanding Method + Godfrey-Isgur Quark Models}
}
```

### Academic References

1. Godfrey, S., & Isgur, N. (1985). "Mesons in a quark model with chromomagnetic interactions." *Physical Review D, 32*(1), 189.

2. Bhatnagar, V., et al. (1995). "Towards a consistent quark model for baryons." *International Journal of Modern Physics A, 10*(03), 335-392.

3. Fulton, R., et al. (1990). "Gaussian wave packets in the Hilbert space formalism." *Physical Review D*.

---

## License

**GEMSTORE** is licensed under the **GNU General Public License v3.0 or later** (GPLv3+).

- **SPDX Identifier**: `GPL-3.0-or-later`
- **Full License**: See `LICENSE` file
- **Copyright**: © 2026 Wen-Xuan Zhang

You are free to:
- ✓ Use for any purpose
- ✓ Modify and redistribute
- ✓ Include in research and commercial products

With the requirement that:
- ⚠ Derivative works must also be licensed under GPLv3+
- ⚠ Source code must be provided
- ⚠ License and copyright notice must be preserved

---

## Contact & Support

- **Author**: Wen-Xuan Zhang ([@serialcore](https://github.com/serialcore))
- **Email**: serialcore@outlook.com
- **Issues**: GitHub Issues tracker
- **AI Assistant**: See `app/gemstore-assistant/SKILL.md`

---

<div align="center">

**Made with ❤️ for computational hadron physics**

```
███████████████████████████████████████
████  GEMSTORE v0.1.4                ████
████  Hadron Spectroscopy Tools      ████
███████████████████████████████████████
```

*Last Updated: April 2026*

</div>
