project('TDSE-2D', 'fortran',
  version : '0.1',
  default_options : ['warning_level=2',
    'optimization=3'
    ]
  )
add_project_arguments([ '-flto=auto' ], language: 'fortran')
#add_project_arguments([ '-fcheck=all', '-g', '-Wall', '-fbacktrace' ], language: 'fortran')
  
subdir('src')
