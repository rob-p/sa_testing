project('sa_testing', ['c', 'cpp'], default_options : ['cpp_std=c++17'])

cmake = import('cmake')

# Configure the CMake project
sais_proj = cmake.subproject('libsais')

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

zdep = dependency('zlib', version : '>=1.2.8')

cm_deps = [sais_lib, kseqpp_lib, zdep]

incdir = include_directories('include')

executable('driver', ['src/main.cpp'], include_directories: incdir, dependencies: cm_deps)
