# H₂⁺ TDSE 2D

**Two-Dimensional Time-Dependent Schrödinger Equation Solver for H₂⁺**

A Fortran-based code for solving the time-dependent Schrödinger equation (TDSE) for the H₂⁺ molecular ion in one and two dimensions (internuclear separation R + electron coordinate x). The code handles nuclear and electronic dynamics under the Born-Oppenheimer framework in the presence of intense, ultrafast laser fields.

---

## Capabilities / Features

### Dimensionality
- **1D Propagation** — Nuclear dynamics only (R coordinate), with multi-state electronic coupling
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
- **Envelope shapes**: cos², Gaussian, and trapezoidal (with configurable rise time)
- **Two-color fields**: Two independent laser pulses with separate wavelengths, intensities, durations, and phases
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
├── input.ini                  # Sample input file with all simulation parameters
├── meson.build                # Top-level Meson build configuration
├── input_data/                # Input potential curves and data files
│   ├── H2+_BO.dat             # Born-Oppenheimer potential curves
│   ├── 12.dat, 13.dat, ...    # Transition dipole moments between states
│   └── Z_interpolated_...     # Soft-core parameter files
└── src/
    ├── main.f90               # Program entry point
    ├── propagation.f90        # 1D propagation module
    ├── propagation_2d.f90     # 2D propagation module
    ├── adiabatic.f90          # Adiabatic surface calculation (ITP)
    ├── nuclear_wv.f90         # Vibrational eigenstate calculation (ITP)
    ├── IO_modules/            # Input/Output and utility modules
    │   ├── input_vars.f90     # Input variable definitions
    │   ├── readinputmodule.f90# Input file parser
    │   ├── global_vars.f90    # Shared simulation variables and arrays
    │   ├── data_au.f90        # Atomic units and conversion constants
    │   ├── commandlinemodule.f90 # Command-line argument parsing
    │   ├── pot_param.f90      # Potential parameters
    │   ├── output_dir.f03     # Output directory management
    │   ├── printinput.f90     # Input parameter printing
    │   └── varprecision.f90   # Precision definitions
    ├── processes/             # Core physics modules
    │   ├── split_operator.f03     # 1D split-operator propagator
    │   ├── split_operator_2d.f03  # 2D split-operator propagator
    │   ├── rk4_operator.f03       # 1D RK4 propagator
    │   ├── rk4_operator_2d.f03    # 2D RK4 propagator
    │   ├── pulse_gen.f90          # Laser pulse generation
    │   ├── setpot.f90             # Potential builder (including KH)
    │   ├── initializer.f90        # Grid and array setup
    │   └── continuum_1d.f03       # Continuum/dissociation analysis
    └── libs/                  # External library interfaces
        ├── fftw3.f90          # FFTW3 Fortran interface
        ├── blas_interface.f03 # BLAS/LAPACK interface
        ├── differentiation.f03# Numerical differentiation
        └── timeit.f90         # Timing utility
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

The compiled executable will be located at `builddir/TDSE-2D` (or `builddir/src/TDSE-2D` depending on the Meson layout).

### 4. Build Options

To build with debugging symbols and runtime checks (slower but useful for development):
```bash
meson setup builddir_debug -Doptimization=0 -Ddebug=true
meson compile -C builddir_debug
```

---

## Usage

Run the simulation by passing an input file to the executable:

```bash
./builddir/TDSE-2D input.ini
```

All simulation parameters are specified in the input file via Fortran namelists. A sample `input.ini` is provided in the repository root.

### Key Input Parameters

#### Grid Definitions
- **`&R_grid`** — `NR`, `Rmin`, `Rmax`: Number of points and bounds for the nuclear (R) grid (in Å)
- **`&x_grid`** — `Nx`, `xmin`, `xmax`: Number of points and bounds for the electronic (x) grid (in Å)
- **`&time_grid`** — `dt`, `Nt`: Time step (a.u.) and total number of time steps

#### Physical Parameters
- **`&nucl_masses`** — `m1`, `m2`: Masses of the two nuclei (in units of proton mass)
- **`&elec_states`** — `Nstates`: Number of electronic states; `sc_kind`: Potential type (`"on_grid"` or `"static"`)
- **`&softcore_params`** — `sc_params`: Soft-core parameter file; `CalcMode`: `"Lab"`, `"KH"`, or `"KH_td"`
- **`&vib_states`** — `guess_vstates`: Number of vibrational states to compute

#### Initial Conditions
- **`&ini_guess_wf`** — `RI`, `kappa`: Center and width of initial Gaussian guess for ITP (in Å)
- **`&ini_state`** — `initial_distribution`: `"single vibrational state"`, `"gaussian distribution"`, or `"Boltzmann distribution"`; `N_ini`, `v_ini`: Initial electronic and vibrational state indices; `RI_tdse`, `kappa_tdse`: Gaussian parameters for TDSE

#### Laser Parameters
- **`&laser_param`** — Two independent laser fields with:
  - `envelope_shape_laser1/2`: `"cos2"`, `"gaussian"`, or `"trapezoidal"`
  - `lambda1/2`: Wavelength (nm)
  - `tp1/2`: Pulse duration (fs)
  - `t_mid1/2`: Center time (fs)
  - `E01/2`: Peak electric field amplitude (a.u.)
  - `phi1/2`: Carrier-envelope phase (in units of π)
  - `rise_time1/2`: Rise time for trapezoidal pulses (fs)

#### Propagation & Parallelization
- **`&propagation_method`** — `propagator_method`: `"split_operator"` or `"rk4"`
- **`&absorber_choice`** — `absorber`: `"mask"` or `"CAP"`
- **`&parallelization`** — `prop_par_FFTW`, `ITP_par_FFTW`: FFTW threading (`"parallel"` or `""`)
- **`&openmp_threads`** — `omp_nthreads`: Number of OpenMP threads (0 = auto)

#### File Paths
- **`&input_files`** — `input_data_dir`: Directory containing input potential/data files
- **`&output_files`** — `output_data_dir`: Directory for output files

---

## Output Files

All output is written to the specified output directory, organized into subdirectories:

- **`time_propagation_data_1d/`** (1D propagation):
  - `norm_1d.out` — Time-dependent wavefunction norm
  - `density_1d_pm3d.out` — Ground-state density map (pm3d format)
  - `ex_density_1d_pm3d.out` — Excited-state density map
  - `avgR_1d.out` — Expectation value ⟨R⟩(t)
  - `vibpop1D_lambda.out` — Vibrational population time-evolution
  - `KER_spectra_from_state_g*.out` — Kinetic Energy Release spectra (normalized and unnormalized)
  - `momt_spectra_from_state_g*.out` — Momentum spectra
  - `psi_outR_*` — Absorbed wavepacket analysis

- **`time_propagation_data_2d/`** (2D propagation):
  - `norm_2d.out` — Time-dependent norm
  - `td-density_R.out`, `td-density_x.out` — Time-dependent R and x density maps
  - `avgR_2d.out`, `avgx_2d.out` — Expectation values ⟨R⟩(t), ⟨x⟩(t)
  - `field_2d.out` — Electric field and vector potential time-history
  - `psi_outR_*` — Absorbed/dissociated wavepacket analysis

- **`nuclear_wavepacket_data/`** — Vibrational eigenstates and energies for each electronic state

---

## References

If you use this code in your research, please cite the relevant publications. *[Add citations / DOIs here]*

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Author

Saurabh Bhatta