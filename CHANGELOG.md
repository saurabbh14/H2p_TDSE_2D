# Changelog

All notable changes to the TDSE-2D solver will be documented in this file.

---

## [0.2.0] — 2026-07-19

### Breaking Change
- **Input format migrated from Fortran namelist to TOML.** The old `input.ini` (17 `&block /` namelists) is replaced by `input.toml`. Parameters are organized into logical sections (`[grid]`, `[system]`, `[time]`, `[initial_guess]`, `[initial_state]`, `[io]`, `[methods]`, `[parallel]`) with laser pulses in a `[[laser.pulses]]` array-of-tables. The legacy `input.ini` is kept for reference but is no longer read by the code.

### Added
- **N-laser support** — Pulse generation now uses an allocatable array of per-pulse parameters. Any number of `[[laser.pulses]]` blocks can be defined, removing the hardcoded 2-laser limit.
- **TOML parser** (`src/IO_modules/toml_parser.f08`) — Self-contained, pure-Fortran TOML subset parser supporting `[section]` headers, `[[array_of_tables]]`, key = value pairs (string/integer/float), `#` comments, and quoted strings.
- **Orphan entry handling** — Key-value pairs placed before the first `[section]` header are collected and reported as a non-fatal warning, with each orphan listed by line number.

### Changed
- **Single-pass input reading** — `readinputmodule.f08` parses the entire TOML file once, extracting all parameters including laser pulses. `pulse_gen.f08` no longer opens or reads the input file independently.
- **`pulse_gen.f08` refactored** — `pulse_param` type holds `single_pulse_data` array instead of hardcoded `laser1`/`laser2` variables. `read_pulse_params` removed; replaced by `initialize_from_lasers(lasers, N_lasers)`.
- **`input_vars.f08`** — Added `LaserParams` derived type and `lasers(:)` allocatable array, plus `N_lasers` count.
- **Source file extension** renamed from `.f90`/`.f03` to `.f08`.
- **`meson.build`** — Added `toml_parser.f08` to IO module sources.
- **`README.md`** — Updated with full TOML input reference, correct usage, and updated directory tree.

### Fixed
- **Parser: section-entry counter bug** — `n_entries` was reset to 0 at every section header, causing entries to be overwritten. Fixed by introducing `entries_this_block` as a separate per-section counter.
- **Parser: array-of-tables deduplication** — Multiple `[[laser.pulses]]` blocks were merged into one table. Fixed by stamping block IDs on entries and grouping by block boundary instead of section name.
- **Input: initial_distribution string mismatch** — TOML default was `"single_vibrational"` but propagation code expects `"single vibrational state"`. Fixed in both `input.toml` and `readinputmodule.f08`.

---

## [0.1] — Initial release
- Fortran namelist input (`input.ini`)
- 1D and 2D propagation (split-operator and RK4)
- Two laser pulses with sin², cos², Gaussian, and trapezoidal envelopes
- Lab, KH, and time-dependent KH frames
- Imaginary Time Propagation for vibrational eigenstates
- Mask and CAP absorbers
- FFTW3 + OpenMP parallelization