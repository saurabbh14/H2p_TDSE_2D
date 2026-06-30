!> Physical constants and unit-conversion factors used throughout the code.
!> All constants are in double precision (dp).
module data_au
    use VarPrecision, only: dp
    real(dp),parameter:: au2a=0.52917706d0  ! length: atomic units -> Angstrom
    real(dp),parameter:: cm2au=4.5554927d-6 ! energy: wavenumbers (cm^-1) -> a.u.
    real(dp),parameter:: au2fs=0.024        ! time: atomic units -> femtoseconds
    real(dp),parameter:: j2eV = 6.242D18    ! energy: Joule -> eV (approx.)
    real(dp),parameter:: au2eV=27.2116d0    ! energy: a.u. -> eV
    real(dp),parameter:: eV2nm=1239.84      ! conversion: eV -> nm (lambda)
    real(dp),parameter:: kB = 3.167d-6      ! Boltzmann constant in hartree/K
    real(dp),parameter:: i2au=2.0997496D-9  ! intensity conversion (W/cm^2 -> a.u.)
    real(dp),parameter:: e02au=5.142206707e11 ! electric field conversion (V/m -> a.u.)
    real(dp),parameter:: pi=3.141592653589793d0
    real(dp),parameter:: mass=1836.15d0     ! electron mass unit reference (e.g. proton mass / electron mass)
    real(dp),parameter:: me =1.d0           ! electron mass 1 a.u.
    real(dp),parameter:: c_speed =137.03604  ! speed of light in a.u. (approx.)
    complex(dp),parameter:: im=(0.d0,1.d0)   ! imaginary unit
end module