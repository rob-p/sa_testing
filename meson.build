project('sa_testing', ['c', 'cpp'], default_options : ['cpp_std=c++17'])

cmake = import('cmake')

# Configure the CMake project

opt_var = cmake.subproject_options()
opt_var.add_cmake_defines({'LIBSAIS_USE_OPENMP':true})
sais_proj = cmake.subproject('libsais', options: opt_var)

fs = import('fs')
if fs.exists('subprojects/libsais.git/VERSION')
  message('found VERSION file in libsais')
  cmd = run_command(
    'rm',
    'subprojects/libsais.git/VERSION',
    check: true 
  )
  cmd = run_command(
    'touch',
    'subprojects/libsais.git/VERSION',
    check: true
  )
else
  message('did not find VERSION file in libsais')
endif

kseqpp_proj = cmake.subproject('kseqpp')

# Fetch the dependency object
sais_lib = sais_proj.dependency('libsais')
kseqpp_lib = kseqpp_proj.dependency('kseq++')

quill_lib = dependency('quill')
zdep = dependency('zlib', version : '>=1.2.8')
# the static linking for OMP seems not do actually statically 
# link it.
omp = dependency('openmp', language: 'c', static: true)
cm_deps = [sais_lib, kseqpp_lib, zdep, quill_lib, omp]

incdir = include_directories('include')

executable('driver', ['src/main.cpp'], 
  include_directories: incdir, 
  dependencies: cm_deps
)
