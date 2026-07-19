!> Reads simulation parameters from a TOML input file.
!> Replaces the old Fortran namelist reader.  All parameters including
!> laser pulses are extracted in a single pass — pulse_gen no longer
!> opens the input file independently.
module ReadInputFile
    use global_vars
    use toml_parser
    implicit none

    type :: InputFilePath
        character(2000) :: path
    contains
        procedure :: read => read_input_file
    end type InputFilePath

contains

    subroutine read_input_file(this)
        class(InputFilePath), intent(inout) :: this
        type(toml_file) :: cfg
        type(toml_error) :: err
        integer :: i
        character(2000) :: buf

        call cfg%parse(this%path, err)
        if (err%flag) then
            print*, "ERROR: ", trim(err%message)
            stop 1
        end if

        ! ----- grid -----
        NR      = cfg%get_int   ("grid", "R_points", 512)
        Rmin    = cfg%get_real  ("grid", "R_min",    0.1_dp)
        Rmax    = cfg%get_real  ("grid", "R_max",   51.2_dp)
        Nx      = cfg%get_int   ("grid", "x_points", 1024)
        xmin    = cfg%get_real  ("grid", "x_min", -102.4_dp)
        xmax    = cfg%get_real  ("grid", "x_max",  102.4_dp)

        ! ----- system -----
        m1            = cfg%get_real  ("system", "m1",            1.0_dp)
        m2            = cfg%get_real  ("system", "m2",            1.0_dp)
        Nstates       = cfg%get_int   ("system", "Nstates",       2)
        buf           = cfg%get_string("system", "sc_kind",       "on_grid")
        sc_kind       = trim(buf)
        buf           = cfg%get_string("system", "sc_params",     "Z_interpolated_optimized_a2_parameters_with_go.txt")
        sc_params     = trim(buf)
        buf           = cfg%get_string("system", "bo_pot_kind",   "on_nuclr_grid")
        bo_pot_kind   = trim(buf)
        buf           = cfg%get_string("system", "CalcMode",      "Lab")
        CalcMode      = trim(buf)
        guess_vstates = cfg%get_int   ("system", "guess_vstates", 100)

        ! ----- time -----
        dt = cfg%get_real("time", "dt", 0.05_dp)
        Nt = cfg%get_int ("time", "Nt", 25000)

        ! ----- initial_guess (ITP) -----
        RI    = cfg%get_real("initial_guess", "RI",    0.7_dp)
        kappa = cfg%get_real("initial_guess", "kappa", -5.0_dp)

        ! ----- initial_state (TDSE) -----
        buf                  = cfg%get_string("initial_state", "distribution", "single_vibrational")
        initial_distribution = trim(buf)
        N_ini      = cfg%get_int  ("initial_state", "N_ini",      1)
        v_ini      = cfg%get_int  ("initial_state", "v_ini",      1)
        RI_tdse    = cfg%get_real ("initial_state", "RI_tdse",    1.5_dp)
        kappa_tdse = cfg%get_real ("initial_state", "kappa_tdse",-5.0_dp)

        ! ----- io -----
        buf              = cfg%get_string("io", "input_data_dir",  "input_data/")
        input_data_dir   = trim(buf)
        buf              = cfg%get_string("io", "output_data_dir", "output_data/")
        output_data_dir  = trim(buf)
        buf              = cfg%get_string("io", "adb_pot",         "H2+_BO.dat")
        adb_pot          = trim(buf)
        buf              = cfg%get_string("io", "trans_dip_prefix","")
        trans_dip_prefix = trim(buf)

        ! ----- methods -----
        buf               = cfg%get_string("methods", "absorber",        "mask")
        absorber          = trim(buf)
        buf               = cfg%get_string("methods", "propagator",      "split_operator")
        propagator_method = trim(buf)
        buf               = cfg%get_string("methods", "gauge",           "length")
        gauge_2d          = trim(buf)
        total_trans_off   = cfg%get_int   ("methods", "total_trans_off", 0)
        buf               = cfg%get_string("methods", "trans_off",       "")
        trans_off         = trim(buf)

        ! ----- parallel -----
        buf           = cfg%get_string("parallel", "prop_fftw",    "parallel")
        prop_par_FFTW = trim(buf)
        buf           = cfg%get_string("parallel", "itp_fftw",     "")
        ITP_par_FFTW  = trim(buf)
        omp_nthreads  = cfg%get_int   ("parallel", "omp_nthreads", 0)

        ! ----- laser pulses (array-of-tables) -----
        N_lasers = cfg%count_array("laser.pulses")
        if (N_lasers > 0) then
            allocate(lasers(N_lasers))
            do i = 1, N_lasers
                buf              = cfg%get_array_string("laser.pulses", i, "envelope",  "sin2")
                lasers(i)%envelope  = trim(buf)
                lasers(i)%lambda    = cfg%get_array_real  ("laser.pulses", i, "lambda",    800._dp)
                lasers(i)%tp        = cfg%get_array_real  ("laser.pulses", i, "tp",          0._dp)
                lasers(i)%t_mid     = cfg%get_array_real  ("laser.pulses", i, "t_mid",       0._dp)
                lasers(i)%alpha0    = cfg%get_array_real  ("laser.pulses", i, "alpha0",      0._dp)
                lasers(i)%phi       = cfg%get_array_real  ("laser.pulses", i, "phi",         0._dp)
                lasers(i)%rise_time = cfg%get_array_real  ("laser.pulses", i, "rise_time",   0._dp)
            end do
        end if

        ! Set default gauge if not specified
        if (gauge_2d == '') gauge_2d = 'length'

        call cfg%finalise()
    end subroutine read_input_file

end module ReadInputFile