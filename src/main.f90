!> Main driver for the Multi-level 1D TDSE solver.
!> This file contains:
!>  - program TDSE_main : entry point that reads input, initializes data,
!>                       computes bound vibrational states, generates pulses
!>                       and runs the real-time propagation.

program TDSE_main

  use CommandLineModule         ! parse command line options
  use ReadInputFile             ! read input.ini (namelists / parameters)
  use PrintInputVars            ! printing helper for chosen simulation variables
  use global_vars               ! shared allocatables and simulation variables
  use pulse_mod                 ! pulse generation and IO
  use data_au                   ! atomic units / unit conversion constants
  use initializer, only: setup => initializer_setup   ! subroutine to initialize grids, arrays
  use adiabatic_mod             ! compute adiabatic surfaces via ITP
  use nuclear_wavefkt           ! compute vibrational eigenstates via ITP
  use propagation1d_mod         ! 1d time propagation using split operator
  use propagation2d_mod        ! 2d time propagation
  use timeit                    ! timeit for checking elapsed time
  implicit none

  ! Local types/objects
  type(CommandLine) :: cmd_line      ! command line handler object
  type(InputFilePath) :: input_path  ! input file path holder
  type(pulse_param) :: pulse         ! pulse object with methods (read, init, gen)
  type(adiabatic_wavefkt_class) :: awf_obj ! adiabatic wavefunction
  type(nuclear_wavefkt_class) :: nwf_obj ! nuclear wavefunction object
  type(time_prop) :: prop
  type(time_prop_2d) :: prop_2d
  type(timer) :: time_this
  
  real(dp) :: Emax                 ! placeholder for energy ranges
  
  ! Start overall timing
  call time_this%start()

  ! Read command line and input file
  call cmd_line%read()
  print*, "reading input:"
  print*, "General Inputs from ", trim(cmd_line%input)
  input_path%path = trim(cmd_line%input)
  call input_path%read()

  ! Read pulse parameters from input file
  call pulse%read(input_path%path)
  print*, "Done reading input"
  print*, "_________________________"

  ! Prepare working directories, grids and arrays
  call setup()

  ! Print and initialize pulse parameters
  print*, "Printing input variables"
  call print_input_vars()
  call pulse%initialize()
  call pulse%param_print()

  call awf_obj%adiabatic_wf_calc()

  ! Allocate arrays for computed vibrational states
  ! Vstates: eigenvalues/energies, chi0: wavefunctions (R, state index, vib index)
  allocate(Vstates(Nstates), chi0(NR,guess_vstates,Nstates))
  chi0 = 0._dp
 
  ! Compute vibrational eigenstates via Imaginary Time Propagation (ITP)
  call nwf_obj%nuclear_wf_calc

  ! Generate and store pulse(s) then run time propagation
  call pulse%generate()

  print*
  print*, "Time dependent calcuations"

!  call prop%propagation_1D(pulse%El, pulse%Al, propagator_method)

  call prop_2d%propagation_2D(pulse%El, pulse%Al, pulse%alpha_t, propagator_method)
!
  ! Clean up allocated arrays and pulse internals
  deallocate (adb, mu_all, chi0)
  call pulse%deallocate_all()

  ! Final timing and report
  call time_this%stop()
  call time_this%elapsed()
  
end program


