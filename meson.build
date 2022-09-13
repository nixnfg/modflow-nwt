project(
  'MODFLOW-NWT',
  'fortran',
 version: '1.3.0',
  meson_version: '>= 0.59.0',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'optimization=2',
    'fortran_std=f2008',
])

if get_option('optimization') == '3'
  error('Only optimization levels <= 2 are supported')
endif

if get_option('optimization') == '2'
  profile = 'release'
else
  profile = 'develop'
endif

message('The used profile is:', profile)


fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
compile_args = []
link_args = []

# Command line options for gfortran
if fc_id == 'gcc'
  # General options
  compile_args += ['-pedantic',
                   '-Wcharacter-truncation',
                   '-Wno-unused-dummy-argument', # This makes problems with OOP
                   '-Wno-intrinsic-shadow',      # We shadow intrinsics with methods, which should be fine
                   '-Wno-maybe-uninitialized',   # "Uninitialized" flags produce false positives with allocatables
                   '-Wno-uninitialized',
                   '-Wno-unused-variable',
                   '-Wno-unused-label',
                   '-Wno-conversion',
                   '-Wno-Obsolescent',
                   ]

  # Options specific to profile
  if profile == 'release'
    compile_args += ['-ffpe-summary=overflow', '-ffpe-trap=overflow,zero,invalid']
  elif profile == 'develop'
    compile_args += ['-fcheck=all', '-ffpe-trap=overflow,zero,invalid,denormal']
  endif

  # Define OS with gfortran for OS specific code
  # These are identical to pre-defined macros available with ifort
  system = build_machine.system()
  if system == 'linux'
    compile_args += '-D__linux__'
  elif system == 'darwin'
    compile_args += '-D__APPLE__'
  elif system == 'windows'
    compile_args += '-D_WIN32'
  endif
endif

# Command line options for ifort
if fc_id == 'intel-cl'
  # windows
  compile_args += ['/fpe:0',              # Activate all floating point exceptions
                   '/fpp',                # Activate preprocessing
                   #'/Qdiag-disable:7416', # f2008 warning
                   #'/Qdiag-disable:7025', # f2008 warning
                   #'/Qdiag-disable:5268', # Line too long
                  ]
  link_args += ['/ignore:4217', # access through ddlimport might be inefficient
                '/ignore:4286'  # same as 4217, but more general
               ]
elif fc_id == 'intel'
  # linux and macOS
  compile_args += ['-fpe0',              # Activate all floating point exceptions
                  '-diag-disable:7416',  # f2008 warning
                  '-diag-disable:7025',  # f2008 warning
                  '-diag-disable:5268',  # Line too long
                  ]
  link_args += '-static-intel'
endif

add_project_arguments(fc.get_supported_arguments(compile_args), language: 'fortran')
add_project_link_arguments(fc.get_supported_arguments(link_args), language: 'fortran')

mpi = dependency('mpi', language : 'fortran', required : false)

# modflow-nwt
sources = files(
    'src' / 'mf_de47_NWT.f90',
    'src' / 'mf_gwf2ag1_NWT_mf.f90',
    'src' / 'mf_gwf2bas7_NWT.f90',
    'src' / 'mf_gwf2bcf7.f90',
    'src' / 'mf_gwf2chd7.f90',
    'src' / 'mf_gwf2drn7_NWT.f90',
    'src' / 'mf_gwf2drt7.f90',
    'src' / 'mf_gwf2ets7.f90',
    'src' / 'mf_gwf2evt7.f90',
    'src' / 'mf_gwf2fhb7.f90',
    'src' / 'mf_gwf2gag7.f90',
    'src' / 'mf_gwf2ghb7_NWT.f90',
    'src' / 'mf_gwf2hfb7_NWT.f90',
    'src' / 'mf_gwf2huf7.f90',
    'src' / 'mf_gwf2hydmod7.f90',
    'src' / 'mf_gwf2ibs7.f90',
    'src' / 'mf_gwf2lak7_NWT.f90',
    'src' / 'mf_gwf2lpf7.f90',
    'src' / 'mf_gwf2mnw17_NWT.f90',
    'src' / 'mf_gwf2mnw27_NWT.f90',
    'src' / 'mf_gwf2mnw2i7.f90',
    'src' / 'mf_gwf2rch7.f90',
    'src' / 'mf_gwf2res7.f90',
    'src' / 'mf_gwf2riv7_NWT.f90',
    'src' / 'mf_gwf2sfr7_NWT.f90',
    'src' / 'mf_gwf2str7.f90',
    'src' / 'mf_gwf2sub7_NWT.f90',
    'src' / 'mf_gwf2swi27.F90',
    'src' / 'mf_gwf2swr7.f90',
    'src' / 'mf_gwf2swr7util.f90',
    'src' / 'mf_gwf2swt7.f90',
    'src' / 'mf_gwf2upw1.f90',
    'src' / 'mf_gwf2uzf1_NWT.f90',
    'src' / 'mf_gwf2wel7_NWT.f90',
    'src' / 'mf_gwflakmodule_NWT.f90',
    'src' / 'mf_gwfsfrmodule_NWT.f90',
    'src' / 'mf_gwfuzfmodule_NWT.f90',
    'src' / 'mf_hufutl7.f90',
    'src' / 'mf_Irestart.f90',
    'src' / 'mf_lmt8_NWT.f90',
    'src' / 'mf_mach_mod.f90',
    'src' / 'mf_mhc7.f90',
    'src' / 'mf_modules.f90',
    'src' / 'mf_nogmg.f90',
    'src' / 'mf_NWT1_gmres.f90',
    'src' / 'mf_NWT1_ilupc_mod.f90',
    'src' / 'mf_NWT1_module.f90',
    'src' / 'mf_NWT1_solver.f90',
    'src' / 'mf_NWT1_xmd.f90',
    'src' / 'mf_NWT1_xmdlib.f90',
    'src' / 'mf_obs2bas7.f90',
    'src' / 'mf_obs2chd7.f90',
    'src' / 'mf_obs2drn7.f90',
    'src' / 'mf_obs2ghb7.f90',
    'src' / 'mf_obs2riv7.f90',
    'src' / 'mf_obs2str7.f90',
    'src' / 'mf_openspec.F90',
    'src' / 'mf_parutl7.f90',
    'src' / 'mf_pcg7_NWT.f90',
    'src' / 'mf_sip7_NWT.f90',
    'src' / 'mf_utl7.f90',
)
src_main = files(
    'src' / 'mf_MF_NWT.f90',
)

ext_sources = static_library('ext_sources', sources)

executable('modflow-nwt', src_main, link_with: [ext_sources], install: true)
