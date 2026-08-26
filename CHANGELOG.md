# Changelog

All notable changes to the TDSE-2D solver will be documented in this file.

---

## [Unreleased]

### Added
- **RK4 evolution mode for the continuum channels** — the absorbed wave packets
  are now integrated with the *same* scheme as the bound wave packet: setting
  `[methods] propagator = "rk4"` switches the 2D ionization and dissociation
  channels (`continuum_2d.f08`) and the 1D absorbed packet (`continuum_1d.f08`)
  from the Strang split-operator to a 4-stage RK4 step, using the vector
  potential at `t`, `t+dt/2` and `t+dt` exactly like `rk4_operator_2d%rk4_step`.
  Each channel prints its evolution mode at start-up. The two doubly-continuum
  channels (dissociation-after-ionization, ionization-after-dissociation) have a
  purely diagonal Hamiltonian in `(PR,Px)`, so their analytic phase is exact and
  is kept in both modes (documented in the code). Note that RK4 is not unitary:
  the channel yields acquire an O(dt⁵)/step drift that the split-operator scheme
  does not have, and each RK4 channel costs ~4 RHS evaluations per step.
- **Explicit unit selection for pulse durations** — new `time_unit` key for laser
  pulses, accepting `"fs"` (default) or `"cycles"` (optical cycles of the pulse
  carrier; aliases `cycle`, `oc`, `optical_cycle`, `optical_cycles`, case
  insensitive). It applies to **both** `tp` and `rise_time`, so the two can no
  longer be given in different units. `t_mid` remains in fs in all cases.
  An optional `[laser]` table sets the default unit for all pulses, which any
  `[[laser.pulses]]` block may override.
- **Effective-duration reporting** — `pulse%param_print()` now prints the optical
  period, the duration unit, and each duration both *as given* and *as actually
  used*, in fs **and** in optical cycles (previously only the requested value was
  printed, which hid the internal cycle rounding).
- **Pulse sanity checks in `initialize`** — non-positive `tp` (cos²/sin²/gaussian)
  or `rise_time` (trapezoidal) and non-positive `lambda` are now detected before
  the (expensive) ITP stage: fatal for a pulse carrying field, and a skipped
  pulse (with a warning) for a zero-amplitude placeholder pulse.

### Changed
- **Continuum API** — `continuum_2d%initialize` and `continuum_1d%initialize`
  gained an *optional* trailing `propagator` argument (the main propagator name),
  and `continuum_2d%propagate` gained optional `A_half` / `A_next` arguments
  (used only in RK4 mode; they default to `A`). Existing call sites without the
  new arguments keep the previous split-operator behaviour.
- **`propagation_2d.f08`** — the half-step and next-step field values
  (`E_half`, `A_half`, `E_next`, `A_next`) are now computed once at the top of the
  time loop and shared by the main RK4 step and the continuum RK4 step, instead of
  being recomputed inside the `case("rk4")` branch (numerically identical).
- **Electron kinetic energy in the 2D main propagators** — see *Fixed*.
- **Duration handling centralised in `initialize_from_lasers`** — the whole-cycle
  rounding that used to be hidden inside `generate_single` is now applied once
  during initialisation and stored in `tp_eff` / `rise_eff`, which the field
  generator uses. Behaviour in `"fs"` mode is unchanged (`n = int(t/T) + 1`
  whole cycles for the sin² width and the trapezoidal rise/fall); in `"cycles"`
  mode the requested duration is used exactly, with no rounding.
- **`single_pulse_data`** — added `time_unit`, `T_cycle`, `tp_eff`, `rise_eff`
  and `skip`; the `cycles_rise/flat/fall/total` counters became `real(dp)` so
  fractional cycle counts are reported honestly.
- **`input_vars.f08`** — `LaserParams` gained `time_unit`; the `tp`/`rise_time`
  comments now state the unit is selected by `time_unit`.
- **`readinputmodule.f08`** — reads `[laser] time_unit` and the per-pulse
  `time_unit`, normalises the spelling and warns (falling back to the default)
  on an unknown unit.
- **`input.toml` / `README.md`** — documented the unit modes; the second
  (zero-amplitude) example pulse now demonstrates `time_unit = "cycles"`.

### Fixed
- **Electron mass in the 2D main propagators** — the electron kinetic term was
  built with mass 1 (`0.5*px²`) in `split_operator_2d%kprop_gen_len/kprop_gen_vel`
  and `rk4_operator_2d%kin_energy_gen/rhs_2d`, while the rest of the 2D code
  (continuum channels, dipole velocity `⟨px⟩/m_eff`, dipole acceleration, the KH
  displacement `kap/m_eff·α`, `setpot`) and the reference implementation
  `check/prop2d.f90` all use the effective electron mass
  `m_eff = (m1+m2)/(m1+m2+1)`. Both propagators now use `px²/(2*m_eff)`
  (velocity gauge: `(px+kap*A)²/(2*m_eff)`). For H₂⁺ this shifts the electron
  kinetic energy by ~0.03 %, so absolute phases differ slightly from earlier runs.
  Note the electronic ITP in `adiabatic.f08` still uses mass 1, so the BO curves
  are (as before) computed with a marginally different Hamiltonian.
- **Out-of-bounds read in the time-dependent KH + RK4 branch** —
  `build_kh_potential_at_time(pot_kh_next, alpha_t(k+1))` read `alpha_t(Nt+1)` on
  the last time step; the next-step quiver displacement is now clamped to
  `alpha_t(Nt)` (mirroring the existing `alpha_half` guard).
- **Misdocumented rise time** — `rise_time` was labelled "fs" in `input.toml`
  while the code silently converted it to a whole number of optical cycles
  (e.g. 4 fs at 228 nm became 6 cycles = 4.53 fs, and at 800 nm 2 cycles =
  5.30 fs). The rounding is unchanged in `"fs"` mode but is now documented and
  reported; `"cycles"` mode gives direct control over the cycle count.
- **Division by zero for `rise_time = 0`** — a trapezoidal pulse with zero
  rise time no longer relies on the accidental `max(n_cycles, 1)` floor; it is
  reported and either skipped (zero amplitude) or rejected.

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