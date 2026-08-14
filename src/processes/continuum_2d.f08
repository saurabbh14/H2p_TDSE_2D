module continuum_2d_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: continuum_2d_type, ionization_prop_type

    !> Propagates the ionized wavefunction in 2D.
    !! Keeps x-coordinate in momentum space (Px), propagates R with split-operator.
    !! The accumulated wavefunction psi_out_acc is stored in (R, Px) representation.
    type :: ionization_prop_type
        private
        ! enable/disable toggle
        logical :: enabled = .false.

        ! absorber positions
        integer :: ix_absorber, iR_absorber
        ! absorber masks
        complex(dp), allocatable :: absorber_x(:)  ! x-absorber mask
        complex(dp), allocatable :: absorber_R(:)  ! R-absorber mask
        ! kinetic propagators
        complex(dp), allocatable :: kprop_x(:)           ! exp(-i*dt*Px²/(2*m_eff))
        complex(dp), allocatable :: kprop_R(:)           ! exp(-i*dt*PR²/(2*m_red))
        ! Coulomb half-step propagator at the x-absorber boundary
        complex(dp), allocatable :: vprop_coul_halfR(:)   ! exp(-i*0.5*dt*pot2D(:,Nx-ix_absorber))
        complex(dp), allocatable :: vprop_coul_halfL(:)   ! exp(-i*0.5*dt*pot2D(:,ix_absorber))

        ! gauge transform for x- (to velocity gauge)
        character(20) :: gauge ! Guage Type: "length" or "velocity" or "KH"
        complex(dp), allocatable :: gauge_transform(:,:)  ! exp(i*A*(x+R)) for velocity gauge

        ! FFTW plans (1D C2C) and work buffers
        type(C_PTR) :: planFx, planBx, planFR, planBR
        type(C_PTR) :: px_in, px_out, pr_in, pr_out
        complex(dp), allocatable :: psi_out_new(:,:)  ! newly ionized wavefunction in (R,Px)
        complex(C_DOUBLE_COMPLEX), pointer :: psi1d_x_in(:), psi1d_x_out(:)
        complex(C_DOUBLE_COMPLEX), pointer :: psi1d_R_in(:), psi1d_R_out(:)
        ! grid sizes
        integer :: NR, Nx
        ! accumulated ionized wavefunction in (R, Px) space
        complex(dp), allocatable, public :: psi_out_acc(:,:)
    contains
        procedure :: initialize => ion_init
        procedure :: set_enabled => ion_set_enabled
        procedure :: extract => ion_extract
        procedure :: accumulate => ion_accum
        procedure :: propagate => ion_propagate
        procedure :: yield => yield
        procedure :: finalize   => ion_finalize
    end type ionization_prop_type

    !> Parent coordinator for all continuum sub-propagators.
    !! Orchestrates ionization, dissociation, and dissociative-ionization propagation.
    type :: continuum_2d_type
        type(ionization_prop_type) :: ion
        ! future: type(dissociation_prop_type)   :: diss
        ! future: type(dissoc_ion_prop_type)     :: diss_ion
    contains
        procedure :: initialize          => cont2d_init
        procedure :: ionization_enable   => cont2d_enable
        procedure :: ionization_disable  => cont2d_disable
        procedure :: ionization_extract  => cont2d_extract
        procedure :: ionization_propagate => cont2d_propagate
        procedure :: ionization_yield    => cont2d_yield
        procedure :: finalize            => cont2d_finalize
    end type continuum_2d_type

contains

    !===========================================================================
    ! ionization_prop_type methods
    !===========================================================================

    !> Initialize FFTW plans, kinetic/Coulomb propagators, absorber phase,
    !! and allocate the accumulated ionized wavefunction array.
    subroutine ion_init(this, ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        use global_vars, only: dt, m_eff, m_red, R, NR, Nx, pot, prop_par_FFTW
        use data_au, only: im
        use FFTW3
        class(ionization_prop_type), intent(inout) :: this
        integer, intent(in) :: ix_absorber       ! grid index of x-absorber boundary
        integer, intent(in) :: iR_absorber       ! grid index of R-absorber boundary
        complex(dp), intent(in) :: abs_x(Nx)      ! x-absorber mask
        complex(dp), intent(in) :: abs_R(NR)      ! R-absorber mask
        character(20), intent(in) :: gauge       ! Gauge Type: "length" or "velocity" or "KH"

        integer :: i

        this%NR = NR
        this%Nx = Nx
        this%ix_absorber = ix_absorber
        this%iR_absorber = iR_absorber
        allocate(this%absorber_x(Nx), this%absorber_R(NR))
        this%absorber_x = abs_x
        this%absorber_R = abs_R
        this%gauge = gauge

        print*
        print*, "Continuum 2D / Ionization: FFTW initialization ..."

        ! Allocate aligned memory for 1D FFTW work buffers
        this%px_in  = fftw_alloc_complex(int(this%Nx, C_SiZE_T))
        call c_f_pointer(this%px_in,  this%psi1d_x_in,  [this%Nx])
        this%px_out = fftw_alloc_complex(int(this%Nx, C_SiZE_T))
        call c_f_pointer(this%px_out, this%psi1d_x_out, [this%Nx])

        this%pr_in  = fftw_alloc_complex(int(this%NR, C_SiZE_T))
        call c_f_pointer(this%pr_in,  this%psi1d_R_in,  [this%NR])
        this%pr_out = fftw_alloc_complex(int(this%NR, C_SiZE_T))
        call c_f_pointer(this%pr_out, this%psi1d_R_out, [this%NR])

        call fftw_initialize_threads

        print*, "Continuum 2D / Ionization: creating 1D C2C plans ..."
        ! planFx: forward FFT along x (row-by-row, reused NR times)
        call fftw_create_c2c_plans(this%psi1d_x_in, this%psi1d_x_out, this%Nx, &
            & this%planFx, this%planBx, prop_par_FFTW)
        ! planFR / planBR: forward/backward FFT along R (col-by-col, reused Nx times)
        call fftw_create_c2c_plans(this%psi1d_R_in, this%psi1d_R_out, this%NR, &
            & this%planFR, this%planBR, prop_par_FFTW)

        ! --- kinetic propagator in x-momentum space ---
        allocate(this%kprop_x(this%Nx))

        ! --- kinetic propagator in R-momentum space ---
        allocate(this%kprop_R(this%NR))

        ! --- Right hand side Coulomb half-step ---
        allocate(this%vprop_coul_halfR(this%NR), this%vprop_coul_halfL(this%NR))
        this%vprop_coul_halfR(:) = exp(-im * 0.5_dp * dt * pot(:, Nx - ix_absorber))
        this%vprop_coul_halfL(:) = exp(-im * 0.5_dp * dt * pot(:, ix_absorber))

        ! --- gauge transform for x- (to velocity gauge) ---
        allocate(this%gauge_transform(this%NR, this%Nx))
        this%gauge_transform(:,:) = 1._dp  ! default: no gauge transform

        ! --- newly ionized wavefunction in (R,Px) space ---
        allocate(this%psi_out_new(this%NR, this%Nx))
        this%psi_out_new = (0._dp, 0._dp)

        ! --- accumulated ionized wavefunction in (R,Px) space ---
        allocate(this%psi_out_acc(this%NR, this%Nx))
        this%psi_out_acc = (0._dp, 0._dp)

        print*, "Continuum 2D / Ionization: Done."

    end subroutine ion_init

    !> Enable or disable ionization tracking.
    subroutine ion_set_enabled(this, state)
        class(ionization_prop_type), intent(inout) :: this
        logical, intent(in) :: state
        
        this%enabled = state
        if (state) then
            print*, "Continuum 2D / Ionization: ENABLED"
        else
            print*, "Continuum 2D / Ionization: DISABLED"
        end if
    end subroutine ion_set_enabled

    !> Extract newly ionized part from psi, FFT along x, and accumulate.
    !! Does NOT modify psi or apply the x-absorber mask.
    subroutine ion_extract(this, psi, A)
        use global_vars, only: NR, Nx, x, R
        use data_au, only: im
        use FFTW3
        class(ionization_prop_type), intent(inout) :: this
        complex(dp), intent(in)    :: psi(NR, Nx)       ! total wavefunction (R,x)
        real(dp), intent(in) :: A

        integer :: i, j
        complex(dp), allocatable :: psi_tmp(:,:)
        real(dp) :: sqrt_Nx

        if (.not. this%enabled) return

        sqrt_Nx = sqrt(dble(Nx))

        ! To be applied to the newly ionized wavefunction before accumulation, to transform to velocity gauge
        do j = 1, Nx
            this%gauge_transform(:,j) = exp(-im * A * (x(j) + R(:)))  ! gauge transform to velocity gauge
        end do

        allocate(psi_tmp(NR, Nx))

        ! Extract: psi_tmp = psi * (1 - abs_x)
        do j = 1, Nx
            psi_tmp(:, j) = psi(:, j) * (1._dp - this%absorber_x(j))
        end do
        
        if (this%gauge == "length") then
            psi_tmp = psi_tmp * this%gauge_transform  ! apply gauge transform
        endif

        ! FFT along x for each row: (R,x) → (R,Px)
        do i = 1, NR
            this%psi1d_x_in(1:Nx)  = psi_tmp(i, 1:Nx)
            call fftw_execute_dft(this%planFx, this%psi1d_x_in, this%psi1d_x_out)
            psi_tmp(i, 1:Nx) = this%psi1d_x_out(1:Nx) / sqrt_Nx
        end do

        this%psi_out_new = psi_tmp

        deallocate(psi_tmp)

    end subroutine ion_extract

    !> Accumulate the newly ionized wavefunction into the total accumulated wavefunction.
    subroutine ion_accum(this)
        class(ionization_prop_type), intent(inout) :: this
         
        this%psi_out_acc = this%psi_out_acc + this%psi_out_new
        this%psi_out_new = (0._dp, 0._dp)
    end subroutine ion_accum

    !> Propagate the accumulated ionized wavefunction in velocity gauge:
    !! x-kinetic + Coulomb split-operator on R.
    subroutine ion_propagate(this, A)
        use global_vars, only: NR, Nx, dt, m_eff, m_red, Px, PR, kap, lam
        use data_au, only: im
        use FFTW3
        class(ionization_prop_type), intent(inout) :: this
        real(dp), intent(in) :: A

        integer :: j
        real(dp) :: sqrt_NR

        if (.not. this%enabled) then
            return
        end if

        ! velocity gauge kinetic propagators
        this%kprop_x = exp(-im * 0.5_dp * dt * (Px + kap * A) * (Px + kap * A) / m_eff) 
        this%kprop_R = exp(-im * 0.5_dp * dt * (PR + lam * A) * (PR + lam * A) / m_red)

        sqrt_NR = sqrt(dble(NR))

        ! 2. Split-operator on R coordinate (per Px column)
        do j = 1, Nx
            ! half-step Coulomb
            if (j > Nx/2) then
                this%psi1d_R_in(1:NR) = this%psi_out_acc(1:NR, j) &
                    & * this%vprop_coul_halfR(1:NR)
            else
                this%psi1d_R_in(1:NR) = this%psi_out_acc(1:NR, j) &
                    & * this%vprop_coul_halfL(1:NR)
            end if
            ! FFT R → PR
            call fftw_execute_dft(this%planFR, this%psi1d_R_in, this%psi1d_R_out)
            ! apply kinetic + normalize forward FFT
            this%psi1d_R_in(1:NR) = this%psi1d_R_out(1:NR) * this%kprop_x(j) &
                & * this%kprop_R(1:NR) / sqrt_NR
            ! iFFT PR → R
            call fftw_execute_dft(this%planBR, this%psi1d_R_in, this%psi1d_R_out)
            this%psi1d_R_in(1:NR) = this%psi1d_R_out(1:NR) / sqrt_NR
            ! half-step Coulomb
            if (j > Nx/2) then
                this%psi_out_acc(1:NR, j)= this%psi1d_R_in(1:NR) *  &
                    & this%vprop_coul_halfR(1:NR)
            else
                this%psi_out_acc(1:NR, j)= this%psi1d_R_in(1:NR) *  &
                    & this%vprop_coul_halfL(1:NR)
            end if
        end do

        ! Apply R-absorber mask (if any) to the accumulated wavefunction
        do j = 1, Nx
            this%psi_out_acc(:, j) = this%psi_out_acc(:, j) * this%absorber_R(:)
        end do

        call this%accumulate()  ! accumulate newly ionized wavefunction

    end subroutine ion_propagate

    !> Returns the ionization yield (norm of accumulated wavefunction).
    subroutine yield(this, ionization_yield)
        use global_vars, only: dR, dx
        class(ionization_prop_type), intent(inout) :: this
        real(dp), intent(out) :: ionization_yield

        if (.not. this%enabled) then
            ionization_yield = 0._dp
            return
        end if

        ! Compute the norm of the accumulated ionized wavefunction
        ionization_yield = sum(abs(this%psi_out_acc)**2) * dR * dx

    end subroutine yield

    !> Clean up FFTW plans, memory, and propagators.
    subroutine ion_finalize(this)
        use FFTW3
        class(ionization_prop_type), intent(inout) :: this

        call fftw_destroy_plan(this%planFx)
        call fftw_destroy_plan(this%planBx)
        call fftw_destroy_plan(this%planFR)
        call fftw_destroy_plan(this%planBR)
        call fftw_free(this%px_in)
        call fftw_free(this%px_out)
        call fftw_free(this%pr_in)
        call fftw_free(this%pr_out)
        if (allocated(this%absorber_x))     deallocate(this%absorber_x)
        if (allocated(this%absorber_R))     deallocate(this%absorber_R)
        if (allocated(this%gauge_transform)) deallocate(this%gauge_transform)
        if (allocated(this%kprop_x))         deallocate(this%kprop_x)
        if (allocated(this%kprop_R))         deallocate(this%kprop_R)
        if (allocated(this%vprop_coul_halfL)) deallocate(this%vprop_coul_halfL)
        if (allocated(this%vprop_coul_halfR)) deallocate(this%vprop_coul_halfR)
        if (allocated(this%psi_out_acc))     deallocate(this%psi_out_acc)
        if (allocated(this%psi_out_new))     deallocate(this%psi_out_new)

        print*, "Continuum 2D / Ionization: cleaned up."

    end subroutine ion_finalize

    !===========================================================================
    ! continuum_2d_type methods (parent coordinator)
    !===========================================================================

    !> Initialize all sub-propagators.
    subroutine cont2d_init(this, ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        use global_vars, only: NR, Nx, pot
        class(continuum_2d_type), intent(inout) :: this
        integer,  intent(in) :: ix_absorber, iR_absorber
        complex(dp), intent(in) :: abs_x(:), abs_R(:)
        character(*), intent(in) :: gauge

        call this%ion%initialize(ix_absorber, iR_absorber, abs_x, abs_R, gauge)

    end subroutine cont2d_init

    !> Enable ionization tracking.
    subroutine cont2d_enable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%ion%set_enabled(.true.)
    end subroutine cont2d_enable

    !> Disable ionization tracking.
    subroutine cont2d_disable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%ion%set_enabled(.false.)
    end subroutine cont2d_disable

    !> Extract newly ionized part.
    !! Must be called BEFORE applying the x-absorber mask on psi.
    subroutine cont2d_extract(this, psi, A)
        class(continuum_2d_type), intent(inout) :: this
        complex(dp), intent(in) :: psi(:,:)
        real(dp), intent(in) :: A
        call this%ion%extract(psi, A)
    end subroutine cont2d_extract

    !> Propagates the accumulated ionized wavefunction and the accumulates the newly ionized wavefunction.
    subroutine cont2d_propagate(this, A)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(in) :: A
        call this%ion%propagate(A)
    end subroutine cont2d_propagate

    !> Return the ionization yield (norm of accumulated wavefunction).
    subroutine cont2d_yield(this, ionization_yield)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(out) :: ionization_yield
        call this%ion%yield(ionization_yield)
    end subroutine cont2d_yield

    !> Finalize all sub-propagators.
    subroutine cont2d_finalize(this)
        class(continuum_2d_type), intent(inout) :: this

        call this%ion%finalize()

    end subroutine cont2d_finalize

end module continuum_2d_mod

