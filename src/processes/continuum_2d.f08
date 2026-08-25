module continuum_2d_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: continuum_2d_type, ionization_prop_type, dissociation_prop_type, free_continuum_prop_type

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

    !> Propagates the dissociated (nuclei flying apart) wavefunction in 2D.
    !! Keeps R-coordinate in momentum space (PR) and x in coordinate space.
    !! The electron remains bound: the propagation includes the electron-nuclear
    !! Coulomb potential at the R-absorber boundary.
    !! The accumulated wavefunction psi_out_acc is stored in (PR, x) representation.
    type :: dissociation_prop_type
        private
        ! enable/disable toggle
        logical :: enabled = .false.

        ! absorber positions
        integer :: ix_absorber, iR_absorber
        ! absorber masks
        complex(dp), allocatable :: absorber_x(:)  ! x-absorber mask (keeps electron bound)
        complex(dp), allocatable :: absorber_R(:)  ! R-absorber mask
        ! kinetic propagators
        complex(dp), allocatable :: kprop_x(:)     ! exp(-i*0.5*dt*(Px+kap*A)^2/m_eff)
        complex(dp), allocatable :: kprop_R(:)     ! exp(-i*0.5*dt*(PR+lam*A)^2/m_red)
        ! Coulomb half-step propagator at the R-dissociation boundary (as function of x)
        complex(dp), allocatable :: vprop_coul_halfR(:)  ! exp(-i*0.5*dt*pot(iR_absorber,:))

        ! gauge transform (to velocity gauge)
        character(20) :: gauge
        complex(dp), allocatable :: gauge_transform(:,:)  ! exp(-i*A*(x+R)) for length->velocity

        ! FFTW plans (1D C2C) and work buffers
        type(C_PTR) :: planFx, planBx, planFR, planBR
        type(C_PTR) :: px_in, px_out, pr_in, pr_out
        complex(dp), allocatable :: psi_out_new(:,:)  ! newly dissociated wavefunction in (PR,x)
        complex(C_DOUBLE_COMPLEX), pointer :: psi1d_x_in(:), psi1d_x_out(:)
        complex(C_DOUBLE_COMPLEX), pointer :: psi1d_R_in(:), psi1d_R_out(:)
        ! grid sizes
        integer :: NR, Nx
        ! accumulated dissociated wavefunction in (PR, x) space
        complex(dp), allocatable, public :: psi_out_acc(:,:)
    contains
        procedure :: initialize => diss_init
        procedure :: set_enabled => diss_set_enabled
        procedure :: extract => diss_extract
        procedure :: accumulate => diss_accum
        procedure :: propagate => diss_propagate
        procedure :: yield => diss_yield
        procedure :: finalize => diss_finalize
    end type dissociation_prop_type

    !> Field-only (Volkov-like) propagator for the doubly-continuum channels
    !! "dissociation-after-ionization" and "ionization-after-dissociation".
    !! No Coulomb potential: the wave packet evolves only under the laser field.
    !! The accumulated wavefunction psi_out_acc is stored in (PR, Px) representation.
    type :: free_continuum_prop_type
        private
        logical :: enabled = .false.
        character(20) :: gauge
        character(40) :: name
        complex(dp), allocatable :: kprop_R(:)     ! exp(-i*0.5*dt*(PR+lam*A)^2/m_red)
        complex(dp), allocatable :: kprop_x(:)     ! exp(-i*0.5*dt*(Px+kap*A)^2/m_eff)
        integer :: NR, Nx
        complex(dp), allocatable, public :: psi_out_acc(:,:)  ! (PR,Px)
    contains
        procedure :: initialize => free_init
        procedure :: set_enabled => free_set_enabled
        procedure :: accumulate => free_accum
        procedure :: propagate => free_propagate
        procedure :: yield => free_yield
        procedure :: finalize => free_finalize
    end type free_continuum_prop_type

    !> Parent coordinator for all continuum sub-propagators.
    !! Orchestrates ionization, dissociation, dissociation-after-ionization and
    !! ionization-after-dissociation propagation.
    type :: continuum_2d_type
        type(ionization_prop_type)     :: ion
        type(dissociation_prop_type)   :: diss
        type(free_continuum_prop_type) :: diss_after_ion
        type(free_continuum_prop_type) :: ion_after_diss
        complex(dp), allocatable :: flux(:,:)  ! scratch (PR,Px) buffer between channels
    contains
        procedure :: initialize            => cont2d_init
        procedure :: ionization_enable     => cont2d_enable
        procedure :: ionization_disable    => cont2d_disable
        procedure :: ionization_extract    => cont2d_extract
        procedure :: ionization_yield      => cont2d_yield
        procedure :: dissociation_enable   => cont2d_diss_enable
        procedure :: dissociation_disable  => cont2d_diss_disable
        procedure :: dissociation_extract  => cont2d_diss_extract
        procedure :: dissociation_yield    => cont2d_diss_yield
        procedure :: diss_after_ion_enable  => cont2d_diss_after_ion_enable
        procedure :: diss_after_ion_disable => cont2d_diss_after_ion_disable
        procedure :: diss_after_ion_yield   => cont2d_diss_after_ion_yield
        procedure :: ion_after_diss_enable  => cont2d_ion_after_diss_enable
        procedure :: ion_after_diss_disable => cont2d_ion_after_diss_disable
        procedure :: ion_after_diss_yield   => cont2d_ion_after_diss_yield
        procedure :: propagate             => cont2d_propagate
        procedure :: finalize              => cont2d_finalize
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
    subroutine ion_propagate(this, A, flux_out)
        use global_vars, only: NR, Nx, dt, m_eff, m_red, Px, PR, kap, lam
        use data_au, only: im
        use FFTW3
        class(ionization_prop_type), intent(inout) :: this
        real(dp), intent(in) :: A
        complex(dp), intent(out), optional :: flux_out(:,:)

        integer :: j
        real(dp) :: sqrt_NR

        if (.not. this%enabled) then
            if (present(flux_out)) flux_out = (0._dp, 0._dp)
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

        ! Harvest the large-R flux (dissociation-after-ionization) BEFORE applying abs_R
        if (present(flux_out)) then
            do j = 1, Nx
                this%psi1d_R_in(1:NR) = this%psi_out_acc(1:NR, j) &
                    & * (1._dp - this%absorber_R(1:NR))
                call fftw_execute_dft(this%planFR, this%psi1d_R_in, this%psi1d_R_out)
                flux_out(1:NR, j) = this%psi1d_R_out(1:NR) / sqrt_NR
            end do
        end if

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
    ! dissociation_prop_type methods
    !===========================================================================

    !> Initialize FFTW plans, kinetic/Coulomb propagators and allocate arrays.
    subroutine diss_init(this, ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        use global_vars, only: dt, NR, Nx, pot, prop_par_FFTW
        use data_au, only: im
        use FFTW3
        class(dissociation_prop_type), intent(inout) :: this
        integer, intent(in) :: ix_absorber, iR_absorber
        complex(dp), intent(in) :: abs_x(Nx), abs_R(NR)
        character(*), intent(in) :: gauge

        this%NR = NR
        this%Nx = Nx
        this%ix_absorber = ix_absorber
        this%iR_absorber = iR_absorber
        this%gauge = gauge
        allocate(this%absorber_x(Nx), this%absorber_R(NR))
        this%absorber_x = abs_x
        this%absorber_R = abs_R

        print*
        print*, "Continuum 2D / Dissociation: FFTW initialization ..."

        this%px_in = fftw_alloc_complex(int(this%Nx, C_SIZE_T))
        call c_f_pointer(this%px_in, this%psi1d_x_in, [this%Nx])
        this%px_out = fftw_alloc_complex(int(this%Nx, C_SIZE_T))
        call c_f_pointer(this%px_out, this%psi1d_x_out, [this%Nx])

        this%pr_in = fftw_alloc_complex(int(this%NR, C_SIZE_T))
        call c_f_pointer(this%pr_in, this%psi1d_R_in, [this%NR])
        this%pr_out = fftw_alloc_complex(int(this%NR, C_SIZE_T))
        call c_f_pointer(this%pr_out, this%psi1d_R_out, [this%NR])

        call fftw_initialize_threads
        print*, "Continuum 2D / Dissociation: FFTW plan creation ..."
        call fftw_create_c2c_plans(this%psi1d_x_in, this%psi1d_x_out, this%Nx, &
            & this%planFx, this%planBx, prop_par_FFTW)
        call fftw_create_c2c_plans(this%psi1d_R_in, this%psi1d_R_out, this%NR, &
            & this%planFR, this%planBR, prop_par_FFTW)

        allocate(this%kprop_x(this%Nx))
        allocate(this%kprop_R(this%NR))

        ! Coulomb half-step in x at the R-dissociation boundary (electron-nuclear)
        allocate(this%vprop_coul_halfR(this%Nx))
        this%vprop_coul_halfR(:) = exp(-im * 0.5_dp * dt * pot(iR_absorber, :))

        allocate(this%gauge_transform(this%NR, this%Nx))
        this%gauge_transform = (1._dp, 0._dp)

        allocate(this%psi_out_new(this%NR, this%Nx))
        this%psi_out_new = (0._dp, 0._dp)
        allocate(this%psi_out_acc(this%NR, this%Nx))
        this%psi_out_acc = (0._dp, 0._dp)

        print*, "Continuum 2D / Dissociation: Done."

    end subroutine diss_init

    !> Enable/disable the dissociation channel.
    subroutine diss_set_enabled(this, state)
        class(dissociation_prop_type), intent(inout) :: this
        logical, intent(in) :: state
        this%enabled = state
    end subroutine diss_set_enabled

    !> Extract the newly dissociated part psi*(1-abs_R) and transform to (PR,x).
    !! Must be called BEFORE applying the R-absorber mask on psi.
    subroutine diss_extract(this, psi, A)
        use global_vars, only: NR, Nx, R, x
        use data_au, only: im
        use FFTW3
        class(dissociation_prop_type), intent(inout) :: this
        complex(dp), intent(in) :: psi(NR, Nx)
        real(dp), intent(in) :: A
        integer :: i, j
        complex(dp), allocatable :: psi_tmp(:,:)
        real(dp) :: sqrt_NR

        if (.not. this%enabled) return

        sqrt_NR = sqrt(dble(NR))

        do j = 1, Nx
            this%gauge_transform(:, j) = exp(-im * A * (x(j) + R(:)))
        end do

        allocate(psi_tmp(NR, Nx))

        ! Extract the part beyond the R-absorber (dissociating wave packet)
        do i = 1, NR
            psi_tmp(i, :) = psi(i, :) * (1._dp - this%absorber_R(i))
        end do

        ! Length gauge -> velocity gauge transform
        if (this%gauge == "length") then
            psi_tmp = psi_tmp * this%gauge_transform
        end if

        ! FFT along R for each column: (R,x) -> (PR,x)
        do j = 1, Nx
            this%psi1d_R_in(1:NR) = psi_tmp(1:NR, j)
            call fftw_execute_dft(this%planFR, this%psi1d_R_in, this%psi1d_R_out)
            psi_tmp(1:NR, j) = this%psi1d_R_out(1:NR) / sqrt_NR
        end do

        this%psi_out_new = psi_tmp
        deallocate(psi_tmp)

    end subroutine diss_extract


    !> Add the newly dissociated wavefunction to the accumulated one.
    subroutine diss_accum(this)
        class(dissociation_prop_type), intent(inout) :: this
        if (.not. this%enabled) return
        this%psi_out_acc = this%psi_out_acc + this%psi_out_new
    end subroutine diss_accum

    !> Propagate the accumulated dissociated wavefunction in velocity gauge:
    !! R-kinetic + electron-nuclear Coulomb split-operator on x.
    !! Optionally harvests the large-|x| flux (ionization-after-dissociation)
    !! into flux_out, already transformed to (PR,Px).
    subroutine diss_propagate(this, A, flux_out)
        use global_vars, only: NR, Nx, dt, m_eff, m_red, Px, PR, kap, lam
        use data_au, only: im
        use FFTW3
        class(dissociation_prop_type), intent(inout) :: this
        real(dp), intent(in) :: A
        complex(dp), intent(out), optional :: flux_out(:,:)

        integer :: i
        real(dp) :: sqrt_Nx

        if (.not. this%enabled) then
            if (present(flux_out)) flux_out = (0._dp, 0._dp)
            return
        end if

        ! velocity gauge kinetic propagators
        this%kprop_x = exp(-im * 0.5_dp * dt * (Px + kap * A)**2 / m_eff)
        this%kprop_R = exp(-im * 0.5_dp * dt * (PR + lam * A)**2 / m_red)

        sqrt_Nx = sqrt(dble(Nx))

        ! Split-operator on x coordinate (per PR row)
        do i = 1, NR
            ! half-step Coulomb
            this%psi1d_x_in(1:Nx) = this%psi_out_acc(i, 1:Nx) * this%vprop_coul_halfR(1:Nx)
            ! FFT x -> Px
            call fftw_execute_dft(this%planFx, this%psi1d_x_in, this%psi1d_x_out)
            ! apply kinetic + normalize forward FFT
            this%psi1d_x_in(1:Nx) = this%psi1d_x_out(1:Nx) * this%kprop_R(i) &
                & * this%kprop_x(1:Nx) / sqrt_Nx
            ! iFFT Px -> x
            call fftw_execute_dft(this%planBx, this%psi1d_x_in, this%psi1d_x_out)
            this%psi1d_x_in(1:Nx) = this%psi1d_x_out(1:Nx) / sqrt_Nx
            ! half-step Coulomb
            this%psi_out_acc(i, 1:Nx) = this%psi1d_x_in(1:Nx) * this%vprop_coul_halfR(1:Nx)
        end do

        ! Harvest the large-|x| flux (ionization-after-dissociation) BEFORE applying abs_x
        if (present(flux_out)) then
            do i = 1, NR
                this%psi1d_x_in(1:Nx) = this%psi_out_acc(i, 1:Nx) * (1._dp - this%absorber_x(1:Nx))
                call fftw_execute_dft(this%planFx, this%psi1d_x_in, this%psi1d_x_out)
                flux_out(i, 1:Nx) = this%psi1d_x_out(1:Nx) / sqrt_Nx
            end do
        end if

        ! Apply x-absorber mask (keep the electron bound)
        do i = 1, NR
            this%psi_out_acc(i, :) = this%psi_out_acc(i, :) * this%absorber_x(:)
        end do

        call this%accumulate()

    end subroutine diss_propagate

    !> Returns the dissociation yield (norm of accumulated wavefunction).
    subroutine diss_yield(this, dissociation_yield)
        use global_vars, only: dR, dx
        class(dissociation_prop_type), intent(inout) :: this
        real(dp), intent(out) :: dissociation_yield
        if (.not. this%enabled) then
            dissociation_yield = 0._dp
            return
        end if
        dissociation_yield = sum(abs(this%psi_out_acc)**2) * dR * dx
    end subroutine diss_yield

    !> Clean up FFTW plans, memory, and propagators.
    subroutine diss_finalize(this)
        use FFTW3
        class(dissociation_prop_type), intent(inout) :: this

        call fftw_destroy_plan(this%planFx)
        call fftw_destroy_plan(this%planBx)
        call fftw_destroy_plan(this%planFR)
        call fftw_destroy_plan(this%planBR)
        call fftw_free(this%px_in)
        call fftw_free(this%px_out)
        call fftw_free(this%pr_in)
        call fftw_free(this%pr_out)
        if (allocated(this%absorber_x))       deallocate(this%absorber_x)
        if (allocated(this%absorber_R))       deallocate(this%absorber_R)
        if (allocated(this%gauge_transform))  deallocate(this%gauge_transform)
        if (allocated(this%kprop_x))          deallocate(this%kprop_x)
        if (allocated(this%kprop_R))          deallocate(this%kprop_R)
        if (allocated(this%vprop_coul_halfR)) deallocate(this%vprop_coul_halfR)
        if (allocated(this%psi_out_acc))      deallocate(this%psi_out_acc)
        if (allocated(this%psi_out_new))      deallocate(this%psi_out_new)

        print*, "Continuum 2D / Dissociation: cleaned up."

    end subroutine diss_finalize


    !===========================================================================
    ! free_continuum_prop_type methods (field-only doubly-continuum channels)
    !===========================================================================

    !> Initialize the field-only propagator.
    subroutine free_init(this, NR, Nx, gauge, name)
        class(free_continuum_prop_type), intent(inout) :: this
        integer, intent(in) :: NR, Nx
        character(*), intent(in) :: gauge
        character(*), intent(in) :: name

        this%NR = NR
        this%Nx = Nx
        this%gauge = gauge
        this%name = name

        allocate(this%kprop_R(NR), this%kprop_x(Nx))
        allocate(this%psi_out_acc(NR, Nx))
        this%psi_out_acc = (0._dp, 0._dp)

        print*, "Continuum 2D / ", trim(adjustl(name)), ": initialized."

    end subroutine free_init

    !> Enable/disable the channel.
    subroutine free_set_enabled(this, state)
        class(free_continuum_prop_type), intent(inout) :: this
        logical, intent(in) :: state
        this%enabled = state
    end subroutine free_set_enabled

    !> Add harvested flux (already in (PR,Px) representation) to the accumulator.
    subroutine free_accum(this, flux)
        class(free_continuum_prop_type), intent(inout) :: this
        complex(dp), intent(in) :: flux(:,:)
        if (.not. this%enabled) return
        this%psi_out_acc = this%psi_out_acc + flux
    end subroutine free_accum

    !> Field-only propagation (kinetic phases only, no Coulomb potential).
    subroutine free_propagate(this, A)
        use global_vars, only: dt, m_eff, m_red, Px, PR, kap, lam
        use data_au, only: im
        class(free_continuum_prop_type), intent(inout) :: this
        real(dp), intent(in) :: A
        integer :: j

        if (.not. this%enabled) return

        this%kprop_R = exp(-im * 0.5_dp * dt * (PR + lam * A)**2 / m_red)
        this%kprop_x = exp(-im * 0.5_dp * dt * (Px + kap * A)**2 / m_eff)

        do j = 1, this%Nx
            this%psi_out_acc(:, j) = this%psi_out_acc(:, j) * this%kprop_R(:) * this%kprop_x(j)
        end do

    end subroutine free_propagate

    !> Returns the channel yield (norm of accumulated wavefunction).
    subroutine free_yield(this, y)
        use global_vars, only: dR, dx
        class(free_continuum_prop_type), intent(inout) :: this
        real(dp), intent(out) :: y
        if (.not. this%enabled) then
            y = 0._dp
            return
        end if
        y = sum(abs(this%psi_out_acc)**2) * dR * dx
    end subroutine free_yield

    !> Clean up memory.
    subroutine free_finalize(this)
        class(free_continuum_prop_type), intent(inout) :: this
        if (allocated(this%kprop_R))     deallocate(this%kprop_R)
        if (allocated(this%kprop_x))     deallocate(this%kprop_x)
        if (allocated(this%psi_out_acc)) deallocate(this%psi_out_acc)
        print*, "Continuum 2D / ", trim(adjustl(this%name)), ": cleaned up."
    end subroutine free_finalize


    !===========================================================================
    ! continuum_2d_type methods (parent coordinator)
    !===========================================================================

    !> Initialize all sub-propagators.
    subroutine cont2d_init(this, ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        use global_vars, only: NR, Nx
        class(continuum_2d_type), intent(inout) :: this
        integer,  intent(in) :: ix_absorber, iR_absorber
        complex(dp), intent(in) :: abs_x(:), abs_R(:)
        character(*), intent(in) :: gauge

        call this%ion%initialize(ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        call this%diss%initialize(ix_absorber, iR_absorber, abs_x, abs_R, gauge)
        call this%diss_after_ion%initialize(NR, Nx, gauge, "dissociation-after-ionization")
        call this%ion_after_diss%initialize(NR, Nx, gauge, "ionization-after-dissociation")

        allocate(this%flux(NR, Nx))
        this%flux = (0._dp, 0._dp)

    end subroutine cont2d_init

    !> Enable/disable ionization tracking.
    subroutine cont2d_enable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%ion%set_enabled(.true.)
    end subroutine cont2d_enable

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

    !> Return the ionization yield (norm of accumulated wavefunction).
    subroutine cont2d_yield(this, ionization_yield)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(out) :: ionization_yield
        call this%ion%yield(ionization_yield)
    end subroutine cont2d_yield

    !> Enable/disable dissociation tracking.
    subroutine cont2d_diss_enable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%diss%set_enabled(.true.)
    end subroutine cont2d_diss_enable

    subroutine cont2d_diss_disable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%diss%set_enabled(.false.)
    end subroutine cont2d_diss_disable

    !> Extract newly dissociated part.
    !! Must be called BEFORE applying the R-absorber mask on psi.
    subroutine cont2d_diss_extract(this, psi, A)
        class(continuum_2d_type), intent(inout) :: this
        complex(dp), intent(in) :: psi(:,:)
        real(dp), intent(in) :: A
        call this%diss%extract(psi, A)
    end subroutine cont2d_diss_extract

    !> Return the dissociation yield.
    subroutine cont2d_diss_yield(this, dissociation_yield)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(out) :: dissociation_yield
        call this%diss%yield(dissociation_yield)
    end subroutine cont2d_diss_yield

    !> Enable/disable the dissociation-after-ionization channel.
    subroutine cont2d_diss_after_ion_enable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%diss_after_ion%set_enabled(.true.)
    end subroutine cont2d_diss_after_ion_enable

    subroutine cont2d_diss_after_ion_disable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%diss_after_ion%set_enabled(.false.)
    end subroutine cont2d_diss_after_ion_disable

    subroutine cont2d_diss_after_ion_yield(this, y)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(out) :: y
        call this%diss_after_ion%yield(y)
    end subroutine cont2d_diss_after_ion_yield

    !> Enable/disable the ionization-after-dissociation channel.
    subroutine cont2d_ion_after_diss_enable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%ion_after_diss%set_enabled(.true.)
    end subroutine cont2d_ion_after_diss_enable

    subroutine cont2d_ion_after_diss_disable(this)
        class(continuum_2d_type), intent(inout) :: this
        call this%ion_after_diss%set_enabled(.false.)
    end subroutine cont2d_ion_after_diss_disable

    subroutine cont2d_ion_after_diss_yield(this, y)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(out) :: y
        call this%ion_after_diss%yield(y)
    end subroutine cont2d_ion_after_diss_yield

    !> Propagate all accumulated continuum channels and transfer the
    !! cross-channel fluxes (diss-after-ion and ion-after-diss).
    subroutine cont2d_propagate(this, A)
        class(continuum_2d_type), intent(inout) :: this
        real(dp), intent(in) :: A

        ! 1. Ionization channel; harvest its large-R flux into diss-after-ion
        call this%ion%propagate(A, this%flux)
        call this%diss_after_ion%accumulate(this%flux)

        ! 2. Pure dissociation channel; harvest its large-|x| flux into ion-after-diss
        call this%diss%propagate(A, this%flux)
        call this%ion_after_diss%accumulate(this%flux)

        ! 3. Field-only propagation of the two doubly-continuum channels
        call this%diss_after_ion%propagate(A)
        call this%ion_after_diss%propagate(A)

    end subroutine cont2d_propagate

    !> Finalize all sub-propagators.
    subroutine cont2d_finalize(this)
        class(continuum_2d_type), intent(inout) :: this

        call this%ion%finalize()
        call this%diss%finalize()
        call this%diss_after_ion%finalize()
        call this%ion_after_diss%finalize()
        if (allocated(this%flux)) deallocate(this%flux)

    end subroutine cont2d_finalize

end module continuum_2d_mod

