module PrintInputVars
    use global_vars
    use pot_param
    use pulse_mod
   
    implicit none
    contains
      subroutine print_input_vars()
        print*, "_________________________"
        print*
        print*, "Input Parameters"
        print*, "_________________________"
        print*
        print'(a,i0)', "Number of R grid points: NR = ", NR
        print'(a,i0)', "Number of x grid points: Nx = ", Nx
        print'(a,f7.2,a,f7.2)', "Masses: m1 = ", m1, ", m2 = ", m2
        print'(a,f5.3,a)', "Time grid: dt = ", dt, " fs"
        print'(a,i0,a)', "Number of time steps: Nt = ", Nt, " steps"
        print'(a,i0)', "Number of electronic states: Nstates = ", Nstates
        print*, "sc parameters kind: sc_kind = ", trim(sc_kind)
        print*, "Reading from file: sc_params = ", trim(sc_params)
        print*, "potential frame: CalcMode = ", trim(CalcMode)
        print*, "Number of maximum considered vibrational states: guess_vstates = ", guess_vstates
        print*
        print*, "Guess vibrational wavefunction (Gaussian): Initial position (RI) = ", RI
        print*, "with initial width (kappa) = ", kappa
        print*
        print*
        print*, "Input and Output Directories:"
        print*, "Input data directory: ", trim(input_data_dir)
        print*, "Output data directory: ", trim(output_data_dir)
        print*, "Transition Dipole switched off: ", total_trans_off
        print*
        print*, "Absorber function: ", absorber
        print*
        print*, "TDSE Initial State:"
        print*, "Mode: ", trim(initial_distribution)
        print*, "Gauge:", trim(gauge_2d)
        print*, "electronic state(s) ", (N_ini-1)
        print*, "vibrational state(s) ", (v_ini-1)
        print*, "Gaussian Distribution TDSE: "
        print*, "centered at RI: ", RI_tdse
        print*, "standard deviation: ", kappa_tdse
        print*
        print*, "FFTW Parallelization:"
        print*, "TDSE Propagation FFTW: ", trim(prop_par_FFTW)
        print*, "ITP FFTW: ", trim(ITP_par_FFTW)
        print*
        print*, "_________________________"
        print*
        print*, "Final grid Parameters"
        print*, "_________________________"
        print*
        print*, "dt = ", SNGL(dt), "a.u."
        print*, "dR = ", SNGL(dR), "a.u."
        print*, "dPR = ", SNGL(dpR), "a.u."
        print*, "dx = ", SNGL(dx), "a.u."
        print*, "dPx = ", SNGL(dpx), "a.u."
        print*, "RI=", sngl(RI), "a.u."
        print*, "R0=", sngl(R0), "a.u.", "Rend=",sngl(Rend), "a.u."
        print*, "x0=", sngl(x0), "a.u.", "xend=",sngl(xend), "a.u."
        print*
        print*, "kap =", kap
        print*, "lam =", lam
        print*
        print*, "__________________________"
        print*
  
      end subroutine print_input_vars
  end module PrintInputVars
