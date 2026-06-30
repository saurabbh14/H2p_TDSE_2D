module global_vars
! Global arrays and derived variables used across modules.
    use InputVars
    real(dp) :: dR, dx                          ! grid spacing in coordinate space
    real(dp), allocatable :: R(:), x(:)         ! coordinate grid
    real(dp), allocatable :: alpha2(:), zeff(:) ! softcore and effective charge parameters
    real(dp), allocatable :: PR(:), Px(:)       ! momentum grid
    real(dp), allocatable :: Pot(:,:), ewf(:,:,:) !2D potential and electronic wavefunction
    real(dp), allocatable :: chi0(:,:,:)        !1D vibrational potential and wavefunctions
    real(dp), allocatable, dimension(:,:,:) :: mu_all ! transition dipole arrays
    real(dp), allocatable, dimension(:,:) :: adb ! adiabatic BO potentials
    real(dp) :: kap, lam                      ! derived coefficients used in dipole / kinetic expressions
    real(dp) :: dt                            ! time step 
    real(dp) :: dpr, dpx                      ! momentum-grid spacing
    real(dp) :: m_eff, m_red                  ! effective and reduced mass
    real(dp) :: mn, mn1, mn2                  ! total mass and individual mass ratios
    real(dp), allocatable:: time(:)         ! time grid
    character(len=2000) :: pulse_data_dir   ! output directory path for pulse data
    character(len=2000) :: ewf_dir   ! output directory path for electronic wavefunction data
    character(len=2000) :: nucl_wf_dir      ! output directory path for nuclear wavefunction data
    character(len=2000) :: time_prop_dir    ! output directory path for time propagation data
    character(len=2000) :: time_prop_dir_1d    ! output directory path for 1d time propagation data
    character(len=2000) :: time_prop_dir_2d    ! output directory path for 2d time propagation data 
end module global_vars
