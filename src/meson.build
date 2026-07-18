blas = dependency('blas')
lapack = dependency('lapack')
omp = dependency('openmp')
fftw3 = dependency('fftw3')
fftw3_omp = declare_dependency(
  link_args: ['-lfftw3_omp', 
    '-fexternal-blas']
)

subdir('IO_modules')
subdir('libs')
subdir('processes')
prog_sources = files('main.f08',
           'nuclear_wv.f08',
           'adiabatic.f08',
           'propagation.f08',
           'propagation_2d.f08')

incdir = include_directories('IO_modules', 
                    'libs',
                    'processes')
executable('TDSE-2D',
            IO_mod_sources
            + libs_sources
            + processes_sources
            + prog_sources,
           dependencies: [ blas, lapack, omp, fftw3, fftw3_omp ],
           install : true,
           include_directories : incdir
           )
