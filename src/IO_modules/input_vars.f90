!> High-level input variables read from input file.
!> These are the *declarative* variables describing the simulation setup.
!> Keep names and comments in sync with readinputmodule.f90.  
module InputVars
    use CommandLineModule
    use VarPrecision, only: dp, idp
    use, intrinsic :: iso_c_binding

    ! R-grid
    integer(C_INT):: NR, Nx                    ! number of grid points (coordinate space)
    real(dp) :: Rmin, xmin                       ! grid minimum (positive) in Angs
    real(dp) :: Rmax, xmax                       ! grid maximum in Angs
    
    ! electronic states
    integer:: Nstates                      ! number of electronic BO states
    character(200):: sc_kind         ! "on_grid" | "static" (select potential source)
    character(2000):: sc_params     ! filename for soft-core parameters if sc_kind = "on_grid"
    character(2000):: CalcMode = "Lab"      ! "Lab" | "KH" (cycle-avg) | "KH_td" (time-dep KH)
    real(dp):: alpha0                     ! alpha0 for KH potential (in Angstrom)

    ! vibrational states 
    integer:: guess_vstates                ! number of vibrational eigenstates to compute
    integer, allocatable:: Vstates(:)      ! storage for computed vibrational energies
    
    ! time grid 
    integer:: Nt                           ! number of time steps for time propagation
    
    ! masses (input in atomic mass units or as specified; converted later)
    real(dp):: m1, m2                      ! masses of particle 1 and 2 (in code units before conversion)
    
    ! guess initial wavefunction
    real(dp):: RI, kappa                   ! RI: center of initial Gaussian (units: Angstrom unless converted)
    
    ! initial TDSE state (how to prepare the initial wavefunction for real-time propagation)
    integer:: N_ini, v_ini                 ! N_ini: electronic state index; v_ini: vibrational quantum number
    integer, allocatable:: v_dist_ini(:)   ! optional explicit vibrational-population vector
    real(dp):: temperature, kappa_tdse, RI_tdse ! parameters used for Boltzmann/Gaussian TDSE initial distributions
    character(2000):: initial_distribution ! string selecting initial distribution type ("single vibrational state", "Boltzmann distribution", etc.)
    
    ! input / output file paths and prefixes
    character(2000):: input_data_dir       ! directory with input grids, dipoles, potentials
    character(2000):: trans_dip_prefix ! trans_dip_prefix: optional prefix for dipole files
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

    ! Propagation method selection
    character(2000):: propagator_method  ! "split_operator" | "rk4"
    
end module InputVars
