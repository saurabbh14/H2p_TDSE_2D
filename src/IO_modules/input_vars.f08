!> High-level input variables read from input file.
!> These are the *declarative* variables describing the simulation setup.
module InputVars
    use VarPrecision, only: dp, idp
    use, intrinsic :: iso_c_binding
    implicit none

    !> Per-laser parameters — one instance per [[laser.pulses]] entry
    type :: LaserParams
        character(150) :: envelope = ""    ! "sin2" | "cos2" | "gaussian" | "trapezoidal"
        real(dp) :: lambda    = 0._dp      ! wavelength (nm)
        real(dp) :: tp        = 0._dp      ! pulse duration (fs)
        real(dp) :: t_mid     = 0._dp      ! pulse midpoint (fs)
        real(dp) :: alpha0    = 0._dp      ! quiver amplitude (a.u.)
        real(dp) :: phi       = 0._dp      ! carrier-envelope phase (units of pi)
        real(dp) :: rise_time = 0._dp      ! rise/fall time (fs) — trapezoidal only
    end type LaserParams

    ! R-grid
    integer(C_INT):: NR, Nx                    ! number of grid points (coordinate space)
    real(dp) :: Rmin, xmin                       ! grid minimum in a.u.
    real(dp) :: Rmax, xmax                       ! grid maximum in a.u.
    
    ! electronic states
    integer:: Nstates                      ! number of electronic BO states
    character(200):: sc_kind         ! "on_grid" | "static" (select potential source)
    character(2000):: sc_params     ! filename for soft-core parameters if sc_kind = "on_grid"
    character(2000):: CalcMode = "Lab"      ! "Lab" | "KH" (cycle-avg) | "KH_td" (time-dep KH)
    character(200):: bo_pot_kind         ! "on_nuclr_grid" | "Morse" (select potential source)
    real(dp):: alpha0                     ! alpha0 for KH potential (in a.u.)

    ! vibrational states 
    integer:: guess_vstates                ! number of vibrational eigenstates to compute
    integer, allocatable:: Vstates(:)      ! storage for computed vibrational energies
    
    ! time grid 
    integer:: Nt                           ! number of time steps for time propagation
    
    ! masses (input in atomic mass units or as specified; converted later)
    real(dp):: m1, m2                      ! masses of particle 1 and 2 (in code units before conversion)
    
    ! guess initial wavefunction
    real(dp):: RI, kappa                   ! RI: center of initial Gaussian (in a.u.)
    
    ! initial TDSE state (how to prepare the initial wavefunction for real-time propagation)
    integer:: N_ini, v_ini                 ! N_ini: electronic state index; v_ini: vibrational quantum number
    integer, allocatable:: v_dist_ini(:)   ! optional explicit vibrational-population vector
    real(dp):: temperature, kappa_tdse, RI_tdse ! TDSE initial distribution parameters (a.u.)
    character(2000):: initial_distribution ! string selecting initial distribution type ("single vibrational state", "Boltzmann distribution", etc.)
    
    ! input / output file paths and prefixes
    character(2000):: input_data_dir       ! directory with input grids, dipoles, potentials
    character(2000):: adb_pot, trans_dip_prefix ! adb_pot: filename for BO surfaces; trans_dip_prefix: optional prefix for dipole files
    character(2000):: output_data_dir      ! directory to write outputs
    
    ! transitions to be switched off (e.g. "12 23")
    integer:: total_trans_off
    character(2000):: trans_off
    
    ! Absorber choice for propagation (mask or CAP)
    character(5):: absorber                ! "mask" | "CAP"
    
    ! FFTW parallelization flags read from input
    character(10):: prop_par_FFTW
    character(10):: ITP_par_FFTW

    ! OMP threads (optional)
    integer :: omp_nthreads

    ! Wavefunction density snapshots (optional 2D dumps)
    logical :: snapshots_enabled = .false.   ! whether to write 2D density snapshots
    integer :: snapshot_frames = 50          ! total number of frames across propagation

    ! Time-dependent density output controls
    integer :: td_density_points = 250       ! number of time points in td density output

    ! Output grid resolution & cropping
    integer :: output_R_stride = 4           ! stride for R-grid in density/snapshot output
    integer :: output_x_stride = 4           ! stride for x-grid in density/snapshot output
    real(dp) :: output_R_start = 0.0_dp      ! fraction of R grid where output starts
    real(dp) :: output_R_end   = 1.0_dp      ! fraction of R grid where output ends
    real(dp) :: output_x_start = 0.25_dp     ! fraction of x grid where output starts
    real(dp) :: output_x_end   = 0.75_dp     ! fraction of x grid where output ends

    ! Propagation method selection
    character(2000):: propagator_method  ! "split_operator" | "rk4"

    ! Gauge selection for propagation
    character(20):: gauge_2d  ! "length" | "velocity"

    ! Laser pulses — allocatable array populated from TOML [[laser.pulses]]
    integer :: N_lasers = 0
    type(LaserParams), allocatable :: lasers(:)
    
end module InputVars