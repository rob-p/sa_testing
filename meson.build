project('sa_testing', ['c', 'cpp'], default_options : ['cpp_std=c++17'])

cmake = import('cmake')

# Configure the CMake project
sais_proj = cmake.subproject('libsais')
kseqpp_proj = cmake.subproject('kseqpp')

# Fetch the dependency object
sais_lib = sais_proj.dependency('libsais')
kseqpp_lib = kseqpp_proj.dependency('kseq++')

cm_deps = [sais_lib, kseqpp_lib]

incdir = include_directories('include')

add_project_arguments(['-lstdc++', '-lm'], language: 'cpp')
executable('driver', ['src/main.cpp'], include_directories: incdir, dependencies: cm_deps)

#https://github.com/IlyaGrebnov/libsais/blob/master/CMakeLists.txt
