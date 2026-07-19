# H₂⁺ TDSE 2D

**Two-Dimensional Time-Dependent Schrödinger Equation Solver for H₂⁺**

A Fortran-based code for solving the time-dependent Schrödinger equation (TDSE) for the H₂⁺ molecular ion in one and two dimensions (internuclear separation R + electron coordinate x). The code handles nuclear and electronic dynamics under the non-Born-Oppenheimer framework in the presence of intense, ultrafast laser fields.

---

## Capabilities / Features

### Dimensionality
- **1D Propagation** — Nuclear dynamics only (R coordinate), with multi-state electronic coupling (Born-Oppenheimer)
- **2D Propagation** — Nuclear + electronic dynamics (R + x coordinates), with full 2D wavefunction evolution

### Electronic Structure
- Born-Oppenheimer potential energy curves and transition dipole moments read from external data files
- Adiabatic electronic state calculation via the `adiabatic` module

### Vibrational States
- Imaginary Time Propagation (ITP) to compute bound vibrational eigenstates for each electronic surface
- Configurable number of vibrational states per electronic state

### Propagation Methods
- **Split-Operator** (Strang splitting) — 2nd-order symplectic integrator
- **4th-order Runge-Kutta (RK4)** — Higher-order explicit integrator

### Laser Pulse Models
- **Envelope shapes**: cos², Gaussian, trapezoidal (with configurable rise time), and CW (continuous wave)
- **Arbitrary number of laser pulses**: Define one, two, or many `[[laser.pulses]]` blocks
- **Carrier-envelope phase** control (in units of π)

### Calculation Modes
- **Lab frame** — Standard length-gauge propagation
- **Kramers-Henneberger (KH) frame** — Cycle-averaged static KH potential
- **Time-dependent KH frame** — Instantaneous KH potential recalculated at each time step

### Initial Wavefunction Distributions
- Single vibrational eigenstate on a chosen electronic surface
- Gaussian distribution (centered at specified R with given width)
- Boltzmann distribution (thermal population of vibrational states)

### Boundary Conditions
- **Complex Absorbing Potential (CAP)** — Complex exponential absorber tuned to an optimal momentum
- **Mask function** — Smooth exponential mask function on one or both grid boundaries

### Parallelization
- **FFTW3 multi-threading** for FFT operations
- **OpenMP parallelization** for interstate coupling matrix operations
- Configurable number of threads via input file

### Observable Outputs
- Time-dependent wavefunction norm and densities (R and x grids)
- Expectation values of position (⟨R⟩, ⟨x⟩)
- Localized-state population analysis (gerade/ungerade)
- Vibrational population analysis on each electronic surface
- Kinetic Energy Release (KER) spectra
- Momentum spectra for continuum/dissociated wavepackets
- Absorbed wavepacket tracking and analysis
- Electric field and vector potential time-histories

---

## Directory Structure

```
.
├── input.toml                  # Sample TOML input file with all simulation parameters
├── input.ini                   # Legacy namelist input (for reference; no longer used)
├── meson.build                 # Top-level Meson build configuration
├── input_data/                 # Input potential curves and data files
│   ├── H2+_BO.dat              # Born-Oppenheimer potential curves
│   ├── 12.dat, 13.dat, ...     # Transition dipole moments between states
│   └── Z_interpolated_...      # Soft-core parameter files
└── src/
    ├── main.f08                # Program entry point
    ├── propagation.f08         # 1D propagation module
    ├── propagation_2d.f08      # 2D propagation module
    ├── adiabatic.f08           # Adiabatic surface calculation (ITP)
    ├── nuclear_wv.f08          # Vibrational eigenstate calculation (ITP)
    ├── IO_modules/             # Input/Output and utility modules
    │   ├── input_vars.f08      # Input variable definitions
    │   ├── readinputmodule.f08 # TOML input file parser
    │   ├── toml_parser.f08     # TOML subset parser (section, array-of-tables)
    │   ├── global_vars.f08     # Shared simulation variables and arrays
    │   ├── data_au.f08         # Atomic units and conversion constants
    │   ├── commandlinemodule.f08 # Command-line argument parsing
    │   ├── pot_param.f08       # Potential parameters
    │   ├── output_dir.f08      # Output directory management
    │   ├── printinput.f08      # Input parameter printing
    │   └── varprecision.f08    # Precision definitions
    ├── processes/              # Core physics modules
    │   ├── split_operator.f08      # 1D split-operator propagator
    │   ├── split_operator_2d.f08   # 2D split-operator propagator
    │   ├── rk4_operator.f08        # 1D RK4 propagator
    │   ├── rk4_operator_2d.f08     # 2D RK4 propagator
    │   ├── pulse_gen.f08           # Laser pulse generation (N-laser support)
    │   ├── setpot.f08              # Potential builder (including KH)
    │   ├── initializer.f08         # Grid and array setup
    │   └── continuum_1d.f08        # Continuum/dissociation analysis
    └── libs/                   # External library interfaces
        ├── fftw3.f08           # FFTW3 Fortran interface
        ├── blas_interface.f08  # BLAS/LAPACK interface
        ├── differentiation.f08 # Numerical differentiation
        └── timeit.f08          # Timing utility
```

---

## Dependencies

| Dependency | Version | Purpose |
|-----------|---------|---------|
| Fortran compiler | GCC ≥ 10 (gfortran) | Compilation; requires LTO and OpenMP support |
| [Meson](https://mesonbuild.com/) | ≥ 1.0 | Build system |
| [FFTW3](http://www.fftw.org/) | ≥ 3.3 | Fast Fourier Transforms |
| OpenMP | Compiler built-in | Shared-memory parallelization |
| BLAS / LAPACK | System-provided | Matrix operations for interstate coupling |

---

## Installation & Build Guide

### 1. Install Dependencies

**Ubuntu / Debian:**
```bash
sudo apt update
sudo apt install gfortran meson ninja-build libfftw3-dev libopenblas-dev liblapack-dev
```

**macOS (Homebrew):**
```bash
brew install gcc meson fftw openblas lapack
```

**HPC Clusters (module-based):**
```bash
module load gcc fftw openblas meson
```

### 2. Clone the Repository

```bash
git clone https://github.com/saurabbh14/H2p_TDSE_2D.git
cd H2p_TDSE_2D
```

### 3. Build with Meson

```bash
meson setup builddir
meson compile -C builddir
```

The compiled executable will be located at `builddir/src/TDSE-2D`.

### 4. Build Options

To build with debugging symbols and runtime checks (slower but useful for development):
```bash
meson setup builddir_debug -Doptimization=0 -Ddebug=true
meson compile -C builddir_debug
```

---

## Usage

Run the simulation by passing a TOML input file to the executable:

```bash
./builddir/src/TDSE-2D input.toml
```

All simulation parameters are specified in a [TOML](https://toml.io/) input file. A sample `input.toml` is provided in the repository root. The legacy `input.ini` (Fortran namelist format) is kept for reference but is no longer used.

### Input File Format (`input.toml`)

The input file is organized into **scalar sections** (`[section]`) and **array-of-tables laser blocks** (`[[laser.pulses]]`). Values can be integers, floats, or quoted strings. Comments use `#`.

#### `[grid]` — Spatial Grid

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `R_points` | int | — | Number of points on the nuclear (R) grid |
| `R_min` | float | Å | Minimum R coordinate |
| `R_max` | float | Å | Maximum R coordinate |
| `x_points` | int | — | Number of points on the electronic (x) grid |
| `x_min` | float | Å | Minimum x coordinate |
| `x_max` | float | Å | Maximum x coordinate |

#### `[system]` — Physical System

| Key | Type | Description |
|-----|------|-------------|
| `m1` | float | Mass of nucleus 1 (amu) |
| `m2` | float | Mass of nucleus 2 (amu) |
| `Nstates` | int | Number of electronic BO states |
| `sc_kind` | string | Soft-core potential: `"on_grid"` or `"static"` |
| `sc_params` | string | Soft-core parameter filename (if `sc_kind = "on_grid"`) |
| `bo_pot_kind` | string | BO potential: `"on_nuclr_grid"` or `"Morse"` |
| `CalcMode` | string | `"Lab"`, `"KH"`, or `"KH_td"` |
| `guess_vstates` | int | Max vibrational states to compute |

#### `[time]` — Time Grid

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `dt` | float | a.u. | Time step |
| `Nt` | int | — | Number of time steps |

#### `[initial_guess]` — ITP Initial Guess

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `RI` | float | Å | Center of initial Gaussian |
| `kappa` | float | — | Width of initial Gaussian (negative) |

#### `[initial_state]` — TDSE Initial Wavefunction

| Key | Type | Description |
|-----|------|-------------|
| `distribution` | string | `"single_vibrational"`, `"gaussian"`, or `"boltzmann"` |
| `N_ini` | int | Initial electronic state index |
| `v_ini` | int | Initial vibrational state index |
| `RI_tdse` | float | Gaussian center for TDSE (Å) |
| `kappa_tdse` | float | Gaussian std dev for TDSE |

#### `[io]` — Input/Output Paths

| Key | Type | Description |
|-----|------|-------------|
| `input_data_dir` | string | Directory with input potential/data files |
| `output_data_dir` | string | Directory for output files |
| `adb_pot` | string | Filename for BO potential curves |
| `trans_dip_prefix` | string | Optional prefix for dipole files |

#### `[methods]` — Propagation & Numerics

| Key | Type | Description |
|-----|------|-------------|
| `absorber` | string | `"mask"` or `"CAP"` |
| `propagator` | string | `"split_operator"` or `"rk4"` |
| `gauge` | string | `"length"` or `"velocity"` |
| `total_trans_off` | int | Number of transitions to switch off |
| `trans_off` | string | Space-separated transition pairs (e.g., `"12 23"`) |

#### `[parallel]` — Parallelization

| Key | Type | Description |
|-----|------|-------------|
| `prop_fftw` | string | TDSE FFTW: `"parallel"` or `""` |
| `itp_fftw` | string | ITP FFTW: `"parallel"` or `""` |
| `omp_nthreads` | int | Number of OpenMP threads (`0` = auto) |

#### `[[laser.pulses]]` — Laser Parameters (Array of Tables)

Add one `[[laser.pulses]]` block per laser pulse. Any number of pulses is supported.

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `envelope` | string | — | `"cos2"`, `"sin2"`, `"gaussian"`, or `"trapezoidal"` |
| `lambda` | float | nm | Wavelength |
| `tp` | float | fs | Pulse duration (FWHM for gaussian; total width for sin²/cos²; flat-top width for trapezoidal) |
| `t_mid` | float | fs | Pulse midpoint time |
| `alpha0` | float | a.u. | Quiver amplitude |
| `phi` | float | π | Carrier-envelope phase (0.0 to 2.0) |
| `rise_time` | float | fs | Rise/fall time (trapezoidal only) |

---

## Output Files

All output is written to the specified output directory, organized into subdirectories:

```
<output_data_dir>/
├── pulse_data/                  # Laser pulse electric field & vector potential
├── nuclear_wavepacket_data/     # Vibrational eigenstates & energies for each electronic state
├── adiabatic_data/              # Adiabatic electronic wavefunction data
└── time_prop/
    ├── 1d/                      # 1D propagation output
    │   ├── norm_1d.out          # Wavefunction norm vs. time
    │   ├── avgR_1d.out          # Expectation value ⟨R⟩(t)
    │   ├── density_1d_pm3d.out  # Ground-state density map (pm3d)
    │   ├── vibpop1D_lambda.out  # Vibrational population time-evolution
    │   └── KER_spectra_from_state_g*.out  # Kinetic Energy Release spectra
    └── 2d/                      # 2D propagation output
        ├── norm_2d.out          # Wavefunction norm vs. time
        ├── avgR_2d.out          # Expectation value ⟨R⟩(t)
        ├── avgx_2d.out          # Expectation value ⟨x⟩(t)
        ├── td-density_R.out     # Time-dependent R density map
        ├── td-density_x.out     # Time-dependent x density map
        └── field_2d.out         # Electric field & vector potential
```

---

## References

TODO

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.