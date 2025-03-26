#CMake >= 3.18: Fortran linkage. This does not affect most mixed Fortran/Kokkos builds.
#New in version 3.8: Added CUDA support.
#In CMake 3.12 and above the find_package() 
#command now searches prefixes specified 
#by the _ROOT CMake variable and the ROOT environment variable. 
# Taishan
#cmake_minimum_required (VERSION 3.16)
cmake_minimum_required (VERSION 3.18)

#project name
set(PROJECT_NAME "LICOMKpp")
project(${PROJECT_NAME})

#find_package(Kokkos REQUIRED)
# find_package(MPI REQUIRED)

enable_language(C CXX Fortran)

# Debug or Release
# set(CMAKE_BUILD_TYPE               Debug)
set(CMAKE_BUILD_TYPE               Release)
set(EXECUTABLE_OUTPUT_PATH         ${CMAKE_SOURCE_DIR}/build/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/build/mod)

# Fortran compiler, GNU or Intel
if (${LICOM_ENABLE_FORTRAN_COMPILER} STREQUAL "LICOM_ENABLE_FORTRAN_COMPILER_INTEL")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W1 -O2 -r8 -i4 -mcmodel=medium -m64 -free -fp-model precise -convert big_endian -assume byterecl -ftz")
    # set(CMAKE_Fortran_FLAGS "-W1 -O2 -r8 -i4 -mcmodel=medium -m64 -free -convert big_endian -assume byterecl -ftz -mieee-fp")
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W1 -O2 -r8 -i4 -mcmodel=medium -m64 -free -fp-model precise -convert big_endian -assume byterecl -ftz -C")

elseif (${LICOM_ENABLE_FORTRAN_COMPILER} STREQUAL "LICOM_ENABLE_FORTRAN_COMPILER_GNU")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium -O2 -ffree-form -ffree-line-length-none -fp-model=precise -fconvert=big-endian -fdefault-real-8 -fdefault-double-8")

if(${MACHINE} STREQUAL "OCEANLIGHT")
# -ffloat-store
    set(CMAKE_Fortran_FLAGS "-O3 -fdefault-double-8 -fdefault-real-8 -fconvert='big-endian' -ffree-line-length-none -mieee-with-inexact -mftz -faddress_align=64 -ffree-form -mieee")
endif ()

endif ()

set(CMAKE_C_FLAGS			"-O2")

set(CMAKE_CXX_STANDARD 17)

if(${MACHINE} STREQUAL "419")
    set(CMAKE_CXX_FLAGS			"-O2 --gpu-max-threads-per-block=1024")
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    set(CMAKE_C_FLAGS			"-O2 -mftz -mieee -faddress_align=64")
    set(CMAKE_CXX_FLAGS			"-O2 -mftz -mieee -faddress_align=64 -std=gnu++1z")
    #set(CMAKE_C_FLAGS			"-Ofast -mftz -mieee -funroll-loops  -floop-nest-optimize -finline -faddress_align=64")
    #set(CMAKE_CXX_FLAGS			"-Ofast -mftz -mieee -funroll-loops  -floop-nest-optimize -finline -faddress_align=64 -std=gnu++1z")
else ()
    set(CMAKE_CXX_FLAGS			"-O3")
    # set(CMAKE_CXX_FLAGS			"-O2 -Wall -ftz -fp-model precise -m64")
endif ()

# set(CMAKE_EXE_LINKER_FLAGS	"${CMAKE_EXE_LINKER_FLAGS} -rdynamic -finstrument-functions")

include_directories(src/head)
include_directories(${CMAKE_SOURCE_DIR}/build/mod)

if(${MACHINE} STREQUAL "IAP_GPU")

# netcdf
    link_directories(/usr/local/netcdf/lib)
    include_directories(/usr/local/netcdf/include)
# MPI
    link_directories(/home/wjl/software/oneAPI/mpi/2021.2.0/lib)
    include_directories(/home/wjl/software/oneAPI/mpi/2021.2.0/include)

# PAPI
    include_directories(/home/wjl/software/PAPI/7.0.1/include)
    link_directories(/home/wjl/software/PAPI/7.0.1/lib)
# GPTL
    include_directories(/home/wjl/software/GPTL/8.1.1/include)
    link_directories(/home/wjl/software/GPTL/8.1.1/lib)
# CUDA
    include_directories(/usr/local/cuda/include)
    link_directories(/usr/local/cuda/lib64)

elseif(${MACHINE} STREQUAL "419")
# MPI
    include_directories(/opt/hpc/software/mpi/hpcx/v2.4.1/intel-2017.5.239/include)
    link_directories(/opt/hpc/software/mpi/hpcx/v2.4.1/intel-2017.5.239/lib)
# NetCDF
    include_directories(/public/software/mathlib/netcdf/4.4.1/intel/include)
    link_directories(/public/software/mathlib/netcdf/4.4.1/intel/lib)
# GPTL
    include_directories(/public/home/esm7/software/GPTL/8.1.1/include)
    link_directories(/public/home/esm7/software/GPTL/8.1.1/lib)

elseif(${MACHINE} STREQUAL "KUNSHAN")

elseif(${MACHINE} STREQUAL "TAISHAN")
# Fortran Compiler
    link_directories(/workspace/public/software/compilers/gcc/kunpenggcc/9.3.1/gcc-9.3.1-2021.03-aarch64-linux/lib)
    link_directories(/workspace/public/software/compilers/gcc/kunpenggcc/9.3.1/gcc-9.3.1-2021.03-aarch64-linux/lib64)

# MPI
    include_directories(/workspace/public/software/mpi/openmpi/4.0.1_kunpenggcc9.3.1/include)
    link_directories(/workspace/public/software/mpi/openmpi/4.0.1_kunpenggcc9.3.1/lib)

# NETCDF
    include_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/netcdf/include)
    link_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/netcdf/lib)

    include_directories(/workspace/public/software/libs/netcdf/netcdf_c_4.7.4_fortran_4.5.3/kunpenggcc9.3.1_hmpi1.1.1/include)
    link_directories(/workspace/public/software/libs/netcdf/netcdf_c_4.7.4_fortran_4.5.3/kunpenggcc9.3.1_hmpi1.1.1/lib)
    
# HDF5
    include_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/hdf5/include)
    link_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/hdf5/lib)

    include_directories(/workspace/public/software/libs/hdf5/1.12.0/kunpenggcc9.3.1_hmpi1.1.1/include)
    link_directories(/workspace/public/software/libs/hdf5/1.12.0/kunpenggcc9.3.1_hmpi1.1.1/lib)
# GPTL
    include_directories(/workspace/home/migration/chengqian/wjl/GPTL/8.1.1/include)
    link_directories(/workspace/home/migration/chengqian/wjl/GPTL/8.1.1/lib)
# PAPI
    include_directories(/workspace/home/migration/chengqian/wjl/PAPI/7.0.1/include)
    link_directories(/workspace/home/migration/chengqian/wjl/PAPI/7.0.1/lib)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
	
    # G++
    include_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/include)
    link_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/lib)

    include_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/include)
    link_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib)
    
    # NetCDF
    include_directories(/usr/sw/yyzlib/local-saca/include)
    link_directories(/usr/sw/yyzlib/local-saca/lib)

    # MPI
    include_directories(/usr/sw/mpi/mpi_current/include)

    # Athread
    include_directories(/usr/sw/swgcc/swgcc710-tools-SEA/shared_include)
    link_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib)

    # include_directories(/usr/sw/penv/include)
    # link_directories(/usr/sw/penv/lib)
    
endif()

file(GLOB SRC-LICOM        LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/*)
file(GLOB SRC-CPPIMPL      LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/cpp-impl/*)
file(GLOB SRC-KOKKOSIMPL   LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/kokkos-impl/*.cpp)
file(GLOB SRC-SWIMPL-HOST  LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/athread-impl/*host.c)
file(GLOB SRC-SWIMPL-SLAVE LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/athread-impl/*slave.c)
file(GLOB SRC-UTIL         LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} src/util/*)

if(${MACHINE} STREQUAL "IAP_GPU")
elseif(${MACHINE} STREQUAL "419")
elseif(${MACHINE} STREQUAL "KUNSHAN")
elseif(${MACHINE} STREQUAL "TAISHAN")
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    set_source_files_properties(${SRC-SWIMPL-HOST} 
      PROPERTIES COMPILE_FLAGS "-mftz -mieee -mhost -msimd -faddress_align=64 -fPIC -O3")
    set_source_files_properties(${SRC-SWIMPL-SLAVE} 
      PROPERTIES COMPILE_FLAGS "-mftz -mieee -mslave -msimd -faddress_align=64 -fPIC -O3")
endif()

message("================================")
message("number of processes: ${NUM_PROCE}")

if(${MACHINE} STREQUAL "IAP_GPU")
    add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos/kokkos-4.1.00-athread)
elseif(${MACHINE} STREQUAL "419")
    add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos/kokkos-4.1.00-athread)
elseif(${MACHINE} STREQUAL "KUNSHAN")
    add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos/kokkos-3.5.00)
elseif(${MACHINE} STREQUAL "TAISHAN")
    add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos/kokkos-4.1.00-athread)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos/kokkos-4.1.00-athread)
endif()

add_executable(${PROJECT_NAME} ${SRC-LICOM} ${SRC-SWIMPL-HOST} ${SRC-SWIMPL-SLAVE} ${SRC-CPPIMPL} ${SRC-KOKKOSIMPL} ${SRC-UTIL})

#set Fortran compiler as the linker
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(${PROJECT_NAME} Kokkos::kokkos netcdff netcdf)

if (${LICOM_ENABLE_GPTL} STREQUAL "define" AND ${LICOM_ENABLE_TEST_TIME} STREQUAL "define")
configure_file(${CMAKE_SOURCE_DIR}/scripts/licom_test_time.hpp.in ${CMAKE_SOURCE_DIR}/src/head/licom_test_time.hpp @ONLY)
target_link_libraries(${PROJECT_NAME} gptl papi)
endif ()

if(${MACHINE} STREQUAL "IAP_GPU")
    target_link_libraries(${PROJECT_NAME} libmpifort.so.12 nvToolsExt mpi)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    target_link_libraries(${PROJECT_NAME} hdf5_hl hdf5 dl z m_slave m -lrt)
    set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "stdc++")
    set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
    set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "-Wl,--as-needed -Wl,-Bdynamic")
endif()


#===============================
# configure file
configure_file(${CMAKE_SOURCE_DIR}/scripts/licom-${LICOM_RES}-config.h.in ${CMAKE_SOURCE_DIR}/src/head/def-undef.h @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/kokkos-config.hpp.in           ${CMAKE_SOURCE_DIR}/src/head/kokkos_config.hpp @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/ocean-${LICOM_RES}-config.in   ${CMAKE_SOURCE_DIR}/run/ocn.parm @ONLY)

if(${MACHINE} STREQUAL "IAP_GPU")
configure_file(${CMAKE_SOURCE_DIR}/scripts/run-iap-gpu.in    ${CMAKE_SOURCE_DIR}/scripts/run      @ONLY)
elseif(${MACHINE} STREQUAL "419")
configure_file(${CMAKE_SOURCE_DIR}/scripts/run-419.in        ${CMAKE_SOURCE_DIR}/scripts/run      @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/dcu-bind-419.in   ${CMAKE_SOURCE_DIR}/scripts/dcu_bind @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/submit-419.in     ${CMAKE_SOURCE_DIR}/scripts/submit   @ONLY)
elseif(${MACHINE} STREQUAL "KUNSHAN")

elseif(${MACHINE} STREQUAL "TAISHAN")
configure_file(${CMAKE_SOURCE_DIR}/scripts/run-taishan.in    ${CMAKE_SOURCE_DIR}/scripts/run    @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/submit-taishan.in ${CMAKE_SOURCE_DIR}/scripts/submit @ONLY)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
configure_file(${CMAKE_SOURCE_DIR}/scripts/run-ocean-light.in ${CMAKE_SOURCE_DIR}/scripts/run    @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/submit-ocean-light.in ${CMAKE_SOURCE_DIR}/scripts/submit @ONLY)
endif()

# configure_file(${CMAKE_SOURCE_DIR}/src/forc_mod.F90-${LICOM_RES}  ${CMAKE_SOURCE_DIR}/src/forc_mod.F90 @ONLY)
# configure_file(${CMAKE_SOURCE_DIR}/src/jra_daily.F90-${LICOM_RES} ${CMAKE_SOURCE_DIR}/src/jra_daily.F90 @ONLY)
#configure_file(${CMAKE_SOURCE_DIR}/src/ssave-cdf-instant.F90-${LICOM_RES} ${CMAKE_SOURCE_DIR}/src/ssave-cdf-instant.F90 @ONLY)
