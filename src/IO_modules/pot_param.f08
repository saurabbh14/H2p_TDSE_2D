!> Parameters related to potentials and absorbers.
module pot_param
    use varprecision, only: dp
    implicit none
    real(dp):: R0, x0     ! Grid-Parameter, start..
    real(dp)::Rend, xend   !..and end

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

        function sc_alpha_parameter_func(r) result(alpha_2)
            use varprecision, only: dp
            real(dp), intent(in):: r
            real(dp) :: c0, c1, r_int, alpha, r_asym, r_asym_0 
            real(dp) :: a1, n1, b1
            real(dp) :: a2, n2, b2
            real(dp) :: a3, n3, b3
            real(dp) :: a4, n4, b4
            real(dp) :: a5, n5, b5
            real(dp) :: c, omega, phi, gamma

            real(dp) :: alpha_2

            c0     = 2.49951040510490136e-01
            c1     = -6.64610412316007260e-01
            r_int  = 4.25306510606404842e-01
            alpha  = 4.10165643142133884e+00
            a1     = 5.71797590213225737e-01
            n1     = 4.35923616858636842e-01
            b1     = 1.38514913001198808e-01
            a2     = 2.22577871377661696e+00
            n2     = 1.21390743555773839e+00
            b2     = 2.97526824296839143e+00
            a3     = -8.24865763318051993e-06
            n3     = 1.99999999999995488e+01
            b3     = 5.41359481166218792e+00
            a4     = -1.45833627482172788e-02
            n4     = 1.14940752709907805e+01
            b4     = 4.31296712754159906e+00
            a5     = 1.43381482600454357e+00
            n5     = 4.43953623254223650e+00
            b5     = 4.09530376115672912e+00
            C      = -9.22940566493708037e-01
            omega  = 1.77316954648702091e-01
            phi    = 5.00863740105063626e-01
            gamma  = 3.39790958199151283e-01

            r_asym = c0 + c1 /(r + 1.0_dp)
            r_asym_0 = c0 + c1

            alpha_2 = r_asym - (r_asym_0 - r_int) * exp(-alpha * (r)) &
                    & + a1 * (r)**n1 * exp(-b1 * (r)) &
                    & + a2 * (r)**n2 * exp(-b2 * (r)) &
                    & + a3 * (r)**n3 * exp(-b3 * (r)) &
                    & + a4 * (r)**n4 * exp(-b4 * (r)) &
                    & + a5 * (r)**n5 * exp(-b5 * (r)) &
                    & + c * exp(-gamma * (r)) * sin(omega * (r) + phi)
        end function

        function sc_zeff_parameter_func(r) result(zeff)
            use varprecision, only: dp
            real(dp), intent(in):: r
            real(dp) :: a, alpha 
            real(dp) :: b, beta
            real(dp) :: zeff

            a = -1.88372604687259493_dp
            alpha = 0.180948402509551687_dp
            b = 1.88652658316634692_dp
            beta = 0.136377938456432868_dp

            zeff = 0.5_dp + a * exp(-alpha * (r)) + b * exp(-beta * (r))

        end function

end module pot_param
