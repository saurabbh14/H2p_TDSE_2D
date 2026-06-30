!> This module reads the variables from the input file
module ReadInputFile
    use global_vars
    use pulse_mod
    implicit none
    type :: InputFilePath
      character(2000) :: path
    contains
      procedure :: read => read_input_file
    end type InputFilePath
  contains
    subroutine read_input_file(this)
      class(InputFilePath), intent(inout) :: this
      integer :: input_tk

      namelist /R_grid/NR,Rmin,Rmax
      namelist /x_grid/Nx,xmin,xmax
      namelist /nucl_masses/m1,m2
      namelist /time_grid/dt,Nt
      namelist /elec_states/Nstates, sc_kind
      namelist /softcore_params/sc_params, CalcMode
      namelist /vib_states/guess_vstates
      namelist /ini_guess_wf/Ri, kappa
      namelist /input_files/input_data_dir, trans_dip_prefix
      namelist /output_files/output_data_dir
      namelist /trans_dip_off/total_trans_off, trans_off
      namelist /absorber_choice/absorber
      namelist /ini_state/v_ini,N_ini,initial_distribution,temperature,kappa_tdse, RI_tdse
      namelist /parallelization/prop_par_FFTW,ITP_par_FFTW
      namelist /openmp_threads/omp_nthreads
      namelist /propagation_method/propagator_method
     
      open(newunit=input_tk, file=adjustl(trim(this%path)), status='old')
      read(input_tk,nml=R_grid)
      read(input_tk,nml=x_grid)
      read(input_tk,nml=nucl_masses)
      read(input_tk,nml=time_grid)
      read(input_tk,nml=elec_states)
      read(input_tk,nml=softcore_params)
      read(input_tk,nml=vib_states)
      read(input_tk,nml=ini_guess_wf)
      read(input_tk,nml=input_files)
      read(input_tk,nml=output_files) 
      read(input_tk,nml=trans_dip_off)
      read(input_tk,nml=absorber_choice)
      read(input_tk,nml=ini_state)
      read(input_tk,nml=parallelization)
      read(input_tk,nml=openmp_threads)
      read(input_tk,nml=propagation_method)

      close(input_tk)

    end subroutine read_input_file
  end module ReadInputFile
