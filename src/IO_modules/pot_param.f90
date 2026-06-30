!> Parameters related to potentials and absorbers.
module pot_param
    use varprecision, only: dp
    implicit none
    real(dp):: R0, x0     ! Grid-Parameter, start..
    real(dp)::Rend, xend   !..and end
    real(dp),parameter:: cpmR= 3.2_dp*2*2 !*2 !absorber position from the end of R-grid
    real(dp),parameter:: cpmx= 25.6_dp !*2 !absorber position from the end of x-grid

    contains
        function morse_potential(de,a,re,r) result(pot)
            use varprecision, only: dp
            ! Simple Morse potential generator used when Elec_pot_kind="Morse".
            ! Parameters:
            !  - de : dissociation energy (in same units as returned pot)
            !  - a  : range parameter controlling width
            !  - re : equilibrium bond length (in same units as r)
            !  - r  : coordinate at which to evaluate potential
            real(dp), intent(in):: de, a, re, r
            real(dp) :: pot

            pot = de * (1._dp - exp(-a * (r-re)))**2

        end function
end module pot_param
