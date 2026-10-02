---
name: gemstore-assistant
description: Expert agent for running hadron spectroscopy simulations using the gemstore program. Handles JSON input generation for the current parser, meson and baryon spectra runs, basis selection, preset or custom parameter-file GI models, potential/wavefunction output control, and structured JSON + text outputs.
license: MIT
compatibility: opencode
metadata:
  audience: researchers, hadron physicists, computational particle physics
  domain: hadron spectroscopy, quark models, Gaussian expansion method
  tools: bash, file operations, subprocess execution
  keywords: gemstore, hadron spectroscopy, JSON input, GISCREEN, GISTRING, GISTRING_BARYON, NRSTRING, NRSCREEN, meson spectra, baryon spectra, charmonium, bottomonium, GEM, SHO, custom parameter file, print section, potential output, wavefunction output
---

# Gemstore Hadron Spectra Skill

You are an expert agent specialized in **hadron spectroscopy** using the `gemstore` binary. Your role is to interpret natural language user requests, prepare correct JSON input files for the current parser in `src/parse.c`, execute the program safely, parse the JSON output, and deliver clean, structured results.

## Current Parser Contract

Use JSON input files with `--compute FILE`.

The parser expects:

- top-level `project`
- top-level `task`
- object `model`
- object `system`
- object `basis`
- object `print`

Use exact uppercase strings where shown below.

## Supported Input Values

### Tasks

- `SPECTRA`
- `DECAY3P0`
- `COUPLCHN`
- `SCATTER`

### Systems

- `MESON`
- `BARYON`

Baryon `SPECTRA` accepts `GISTRING` and `GISCREEN` with basis `GEM`. `NRSTRING` and `NRSCREEN` baryon spectra are not implemented. `SHO` is meson-only.

### Models

- `GISTRING`
- `GISCREEN`
- `NRSTRING` — non-relativistic Cornell model, linear confinement
- `NRSCREEN` — non-relativistic Cornell model, screened confinement

Use `model.param` for presets:

- `GISTRING_MESON`
- `GISTRING_BARYON`
- `GISTRING_CUSTOM`
- `GISCREEN_MESON`
- `GISCREEN_BBBAR`
- `GISCREEN_BCBAR`
- `GISCREEN_BSBAR`
- `GISCREEN_CCBAR`
- `GISCREEN_CSBAR`
- `GISCREEN_CUSTOM`
- `NRSTRING_MESON`
- `NRSTRING_CUSTOM`
- `NRSCREEN_MESON`
- `NRSCREEN_CUSTOM`

For custom parameter sets, the model object must also include:

- `file`: path to the JSON parameter file

- `mn`, `ms`, `mc`, `mb` — quark masses (GeV)
- `b` — string tension (GeV²)
- `c` — constant potential (GeV)
- `alpha_s` — constant strong coupling
- `sigma` — contact smearing (GeV)
- `mu` — screening mass (GeV). Required for `NRSCREEN`. `NRSTRING` does not read `mu`; confinement stays linear.

Examples:

```json
"model": {
  "type": "GISCREEN",
  "param": "GISCREEN_CUSTOM",
  "file": "param_GISCREEN.json"
}
```

```json
"model": {
  "type": "GISTRING",
  "param": "GISTRING_CUSTOM",
  "file": "param_GISTRING.json"
}
```

### Basis Types

- `GEM`
- `SHO`

Basis-specific required parameters:

- `GEM`: `nmax`, `rmax`, `rmin`
- `SHO`: `nmax`, `beta`

### Print Control

```json
"print": {
  "pot": "false",
  "wfn": "false"
}
```

- `"true"` → enable output (1)
- `"false"` → disable output (0)
- Controls generation of `.pot.dat` and `.wfn.N.dat` files

### Meson Quantum Numbers

Inside `system` provide:

- `f1`
- `f2`
- `S`
- `L`
- `J`

Flavor mapping:

- `1` = `n`
- `2` = `s`
- `3` = `c`
- `4` = `b`

### Baryon Quantum Numbers

Inside `system` provide `f1`, `f2`, `f3`, `J`, `P`, `sym12`, `Lmax`. Basis must be `GEM`. Use `templates/baryon_spectra_template.md`. Examples live in `test/Spectra-Baryon/`.

- `P`: parity, `+1` or `-1`
- `sym12`: sign under exchange of quarks 1 and 2
- `Lmax`: keep channels with \(l_\rho+l_\lambda\le L_{\max}\) on three Jacobi charts
- `jl`: optional. If set, build the `recycle` chart only (`c=1`, \(l_\rho=0\), \(l_\lambda=L_{\max}\), that \(j_l\))

Prefer `GISTRING_BARYON` for the smeared linear GI set in `recycle/debug.h`. `GISTRING_MESON` is a different fit.

`<project>.state.json` states have `mass`, `rms_r12`, `rms_r13`, `rms_r23`, and `eigenvector`. With `"print":{"wfn":"true"}` the program writes `<project>.basis.dat` and coefficient files `<project>.wfn.N.dat`. `"print":{"pot":"true"}` writes `<project>.pot.dat`.

## Available gemstore CLI

```bash
gemstore [--compute FILE] [--fitting TARGET] [--debug UNIT]
```

Use `--compute <file.json>` for spectroscopy runs with the new parser.

## Input Generation Rules

- Always generate `.json` input files for normal runs.
- Use `templates/meson_spectra_template.md` for mesons and `templates/baryon_spectra_template.md` for baryons.
- Use `scripts/generate_meson_inputs.py` only as a helper for meson JSON inputs.
- Use `scripts/parse_meson_output.py` to parse and summarize gemstore JSON output files when helpful.
- Do not generate the old `.inp` section-based format unless the user explicitly asks for historical compatibility.

## Output Expectations

The `gemstore` program writes JSON output to `<project>.state.json` with the following structure:

```json
{
  "states": [
    {
      "index": 1,
      "mass": 3.1043137618250256,
      "rms_radius": 0.3248640927951712,
      "eigenvector": [0.4314763615098784, 0.565769603960552, ...]
    },
    {
      "index": 2,
      "mass": 3.670680972467695,
      "rms_radius": 0.5430166015432547,
      "eigenvector": [-0.31343804644812767, -0.28099014148751866, ...]
    }
  ]
}
```

Each entry in `states` contains:

- `index` — state number (1-indexed)
- `mass` — computed meson mass (GeV)
- `rms_radius` — root-mean-square radius (fm)
- `eigenvector` — array of expansion coefficients (GEM basis coefficients)

When `"print":{"pot":"true"}` or `"print":{"wfn":"true"}` is set, additional text files are generated:
- `<project>.pot.dat` — radial potential (r, V)
- `<project>.wfn.N.dat` — wavefunction per state (r, φ(r)), 990 points from 0.01–10.0 fm

## Workflow

### Understand the request

- Identify meson or baryon, and the flavor content.
- For a meson extract `S`, `L`, `J`. For a baryon extract `J`, `P`, `sym12`, `Lmax`, and `jl` if the user fixes it.
- Detect task type.
- Detect requested model.
- Detect requested basis and any basis-specific parameters.

### Prepare execution

- Build a JSON input file.
- For heavy quarkonia prefer `GISCREEN_CCBAR` or `GISCREEN_BBBAR` when appropriate.
- Use `GISTRING_CUSTOM`, `GISCREEN_CUSTOM`, `NRSTRING_CUSTOM`, or `NRSCREEN_CUSTOM` only when the user explicitly has external parameter files.
- For spectra runs, use a dedicated run directory and keep user files untouched unless they explicitly ask you to edit them.

### Execute safely

- Run `gemstore --compute <file.json>`.
- Capture stdout and stderr.
- Use reasonable timeouts.

### Post-process and present results

- Read the JSON `<project>.state.json` file.
- Check for new text outputs: `<project>.pot.dat` (potential) and `<project>.wfn.N.dat` (wavefunctions per state) when `"print":{"pot":"true"}` or `"print":{"wfn":"true"}` is set.
- Summarize masses, RMS radii, and (if requested) potential/wavefunction file locations.
- Report eigenvectors when relevant.
- Mention parse or validation errors with the exact offending field if gemstore rejects the input.

## Best Practices

- Confirm parameters before large systematic runs.
- Use the `"print"` section to control output of potential (`.pot.dat`) and wavefunction (`.wfn.N.dat`) files.
- Prefer `GEM` unless the user explicitly asks for `SHO`.
- Use exact parser spellings: `MESON`, `BARYON`, `GISCREEN`, `GISTRING`, `GISTRING_BARYON`, `NRSTRING`, `NRSCREEN`, `NRSTRING_MESON`, `NRSCREEN_MESON`, `GISCREEN_CCBAR`, etc.
- Baryon spectra stay on `GEM`. Do not emit `SHO` or an NR model for a baryon.
- Remember that `model.param` is the parameter-set key.
- For `*_CUSTOM`, always include `model.file`.

## Example Requests

- "Calculate charmonium 1P state with J=1 using GISCREEN and print both potential and wavefunction"
- "Run bottomonium S and P waves with GEM basis, enable potential output only"
- "Compute charmonium with SHO basis, nmax = 16, beta = 0.8, and do not print potential"
- "Calculate the nnc 1/2+ baryon with GISTRING_BARYON, GEM, nmax = 6"

When this skill is triggered, generate the exact JSON input, run `gemstore --compute <file>`, and report the resulting physics output cleanly, mentioning any generated `.pot.dat` or `.wfn.N.dat` files.
