INCLUDE(IFEMOptions)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Custom profiles
IF(${CMAKE_BUILD_TYPE} MATCHES "Nopt")
  SET(CMAKE_BUILD_TYPE Debug)
ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Nomp")
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Nopt")

# IFEM includes Fortran code
ENABLE_LANGUAGE(Fortran)
IF(CMAKE_Fortran_COMPILER MATCHES ifort)
  SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} -lifcore)
ENDIF(CMAKE_Fortran_COMPILER MATCHES ifort)

# And C++ code
ENABLE_LANGUAGE(CXX)
IF(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DUSE_MKL -mkl=sequential")
  SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DUSE_MKL -mkl=sequential")
  FIND_PACKAGE(CBLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSE(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  FIND_PACKAGE(CBLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
  SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} ${CBLAS_DEFINITIONS}")
  SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} ${CBLAS_DEFINITIONS}")
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES Intel)

if(MINGW)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -static-libgfortran -static")
endif()

INCLUDE(TestCXXAcceptsFlag)

IF(${CMAKE_BUILD_TYPE} MATCHES "Release")
  CHECK_CXX_ACCEPTS_FLAG("-mtune=native" HAVE_MTUNE)
  IF (HAVE_MTUNE)
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -mtune=native")
  ENDIF(HAVE_MTUNE)
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Release")

IF(IFEM_WHOLE_PROG_OPTIM)
  CHECK_CXX_ACCEPTS_FLAG("-flto" HAVE_LINK_OPTS)
  IF(HAVE_LINK_OPTS)
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -flto")
  ENDIF(HAVE_LINK_OPTS)
ENDIF(IFEM_WHOLE_PROG_OPTIM)

# Required dependences
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(GoTrivariate REQUIRED)
FIND_PACKAGE(ARPACK REQUIRED)
find_package(CXX11 REQUIRED)

# Mimimum GoTools version
IF(GoTools_VERSION_MAJOR LESS 3 OR NOT GoTools_VERSION_MAJOR)
  MESSAGE(FATAL_ERROR "GoTools >= 3.0.0 required. bailing")
ENDIF(GoTools_VERSION_MAJOR LESS 3 OR NOT GoTools_VERSION_MAJOR)

SET(IFEM_DEPLIBS ${IFEM_DEPLIBS}
                 ${GoTrivariate_LIBRARIES}
                 ${GoTools_LIBRARIES}
                 ${ARPACK_LIBRARIES}
                 ${LAPACK_LIBRARIES}
                 ${CBLAS_LIBRARIES})

SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES}
                     ${GoTools_INCLUDE_DIRS}
                     ${GoTrivariate_INCLUDE_DIRS})

SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CXX_STD11_FLAGS}")
SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CXX_STD11_FLAGS}")

# HDF5
IF(IFEM_USE_HDF5)
  FIND_PACKAGE(HDF5)
  IF(HDF5_LIBRARIES AND HDF5_INCLUDE_DIR)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${HDF5_LIBRARIES})
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${HDF5_INCLUDE_DIR})
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_HDF5=1")
  ENDIF(HDF5_LIBRARIES AND HDF5_INCLUDE_DIR)
ENDIF(IFEM_USE_HDF5)

# TinyXML
FIND_PACKAGE(TinyXML)
IF(TINYXML_INCLUDE_DIR AND TINYXML_LIBRARIES)
  SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${TINYXML_INCLUDE_DIR})
  SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${TINYXML_LIBRARIES})
  SET(IFEM_USE_SYSTEM_TINYXML 1)
ENDIF(TINYXML_INCLUDE_DIR AND TINYXML_LIBRARIES)

# SuperLU
IF(IFEM_USE_SUPERLU OR IFEM_USE_SUPERLU_MT)
  FIND_PACKAGE(SuperLU)
  IF(SuperLU_MT_LIBRARIES AND SuperLU_MT_INCLUDES AND IFEM_USE_SUPERLU_MT)
    FIND_PACKAGE(Threads REQUIRED)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${SuperLU_MT_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${SuperLU_MT_INCLUDES})
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_SUPERLU_MT")
    MESSAGE(STATUS "Using SuperLU-MT")
  ELSEIF(SuperLU_LIBRARIES AND SuperLU_INCLUDES)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${SuperLU_LIBRARIES})
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${SuperLU_INCLUDES})
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_SUPERLU")
    MESSAGE(STATUS "Using SuperLU")
  ENDIF(SuperLU_MT_LIBRARIES AND SuperLU_MT_INCLUDES AND IFEM_USE_SUPERLU_MT)
ENDIF(IFEM_USE_SUPERLU OR IFEM_USE_SUPERLU_MT)

# LR splines
IF(IFEM_USE_LRSPLINES)
  FIND_PACKAGE(LRSpline)
  IF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)
    IF(LRSPLINE_VERSION_MINOR LESS 5 AND LRSPLINE_VERSION_MAJOR LESS 1)
      MESSAGE(STATUS "LR-spline library is out of date.
                      Found ${LRSPLINE_VERSION_MAJOR}.${LRSPLINE_VERSION_MINOR}.${LRSPLINE_VERSION_PATCH}, need at least 0.5.0.
                      Support not enabled")
      SET(LRSpline_LIBRARIES "")
      SET(LRSpline_INCLUDE_DIRS "")
    ELSE(LRSPLINE_VERSION_MINOR LESS 5 AND LRSPLINE_VERSION_MAJOR LESS 1)
      SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${LRSpline_LIBRARIES})
      SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${LRSpline_INCLUDE_DIRS})
      SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_LRSPLINE=1 ${LRSpline_DEFINITIONS}")
      SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DHAS_LRSPLINE=1 ${LRSpline_DEFINITIONS}")
    ENDIF(LRSPLINE_VERSION_MINOR LESS 5 AND LRSPLINE_VERSION_MAJOR LESS 1)
  ENDIF(LRSpline_LIBRARIES AND LRSpline_INCLUDE_DIRS)
ENDIF(IFEM_USE_LRSPLINES)

# PETSc
IF(IFEM_USE_PETSC)
  FIND_PACKAGE(Petsc)
  IF(PETSC_FOUND)
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${PETSC_INCLUDE_DIRS})
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${PETSC_LIBRARIES})
    SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DHAS_PETSC")
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_PETSC")
    IF(IFEM_USE_PARALLEL_PETSC)
      IF(IFEM_USE_OPENMP)
        MESSAGE(FATAL_ERROR "Cannot use parallel PETSc and OpenMP in combination, bailing")
      ENDIF(IFEM_USE_OPENMP)
      SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DPARALLEL_PETSC") # Needed due to usage in apps..
      SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DPARALLEL_PETSC")
    ENDIF(IFEM_USE_PARALLEL_PETSC)
    IF(IFEM_ENABLE_SLEPC)
      FIND_PACKAGE(SLEPc)
      IF(SLEPC_LIBRARIES AND SLEPC_INCLUDES)
        SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${SLEPC_INCLUDES})
        SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${SLEPC_LIBRARIES})
        SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_SLEPC")
      ENDIF(SLEPC_LIBRARIES AND SLEPC_INCLUDES)
    ENDIF(IFEM_ENABLE_SLEPC)
  ENDIF()
ENDIF(IFEM_USE_PETSC)

# SPR
IF(IFEM_USE_SPR)
  FIND_PACKAGE(SPR)
  IF(SPR_LIBRARIES)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${SPR_LIBRARIES})
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_SPR -DUSE_F77SAM")
  ENDIF(SPR_LIBRARIES)
ENDIF(IFEM_USE_SPR)

# SAMG
IF(IFEM_USE_SAMG)
  FIND_PACKAGE(SAMG)
  IF(SAMG_LIBRARIES AND SAMG_INCLUDES)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${SAMG_LIBRARIES})
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${SAMG_INCLUDES})
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_SAMG -DSAMG_UNIX_LINUX=1 -DSAMG_LCASE_USCORE")
  ENDIF(SAMG_LIBRARIES AND SAMG_INCLUDES)
ENDIF(IFEM_USE_SAMG)

# VTF writer
IF(IFEM_USE_VTFWRITER)
  FIND_PACKAGE(VTFWriter)
  IF(VTFWRITER_LIBRARIES AND VTFWRITER_INCLUDES)
    SET(IFEM_DEPLIBS ${IFEM_DEPLIBS} ${VTFWRITER_LIBRARIES})
    SET(IFEM_DEPINCLUDES ${IFEM_DEPINCLUDES} ${VTFWRITER_INCLUDES})
    SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DHAS_VTFAPI=${VTFAPI}")
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} -DHAS_VTFAPI=${VTFAPI}")
  ENDIF(VTFWRITER_LIBRARIES AND VTFWRITER_INCLUDES)
ENDIF(IFEM_USE_VTFWRITER)

# OpenMP
IF(IFEM_USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP")
    SET(IFEM_BUILD_CXX_FLAGS "${IFEM_BUILD_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP")
  ENDIF(OPENMP_FOUND)
ENDIF(IFEM_USE_OPENMP)

# Portability issues
include(CheckFunctionExists)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_FLAGS)
check_function_exists(get_current_dir_name unistd.h HAVE_GET_CURRENT_DIR_NAME) # lacks on osx
if(HAVE_GET_CURRENT_DIR_NAME)
  set(IFEM_CXX_FLAGS "${IFEM_CXX_FLAGS} -DHAVE_GET_CURRENT_DIR_NAME=1")
endif()

# Get GTest tests as CMake tests.
# Copied from FindGTest.cmake
# Thanks to Daniel Blezek <blezek@gmail.com> for the GTEST_ADD_TESTS code
function(gtest_add_tests executable working_dir)
    if(NOT ARGN)
        message(FATAL_ERROR "Missing ARGN: Read the documentation for GTEST_ADD_TESTS")
    endif()
    if(NOT UNIT_TEST_NUMBER)
      set(UNIT_TEST_NUMBER 0 CACHE INT "" FORCE)
    endif()
    foreach(source ${ARGN})
        file(READ "${source}" contents)
        string(REGEX MATCHALL "TEST_?[F]?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
        foreach(hit ${found_tests})
            string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+), *([A-Za-z_0-9]+) *\\).*" "\\1.\\2" test_name ${hit})
            math(EXPR UNIT_TEST_NUMBER "${UNIT_TEST_NUMBER}+1")
            set(UNIT_TEST${UNIT_TEST_NUMBER} ${test_name} ${working_dir} ${executable} --gtest_filter=${test_name} CACHE STRING "" FORCE)
        endforeach()
        # Groups parametrized tests under a single ctest entry
        string(REGEX MATCHALL "INSTANTIATE_TEST_CASE_P\\(([^,]+), *([^,]+)" found_tests2 ${contents})
        foreach(hit ${found_tests2})
          string(SUBSTRING ${hit} 24 -1 test_name)
          string(REPLACE "," ";" test_name "${test_name}")
          list(GET test_name 0 filter_name)
          list(GET test_name 1 test_prefix)
          string(STRIP ${test_prefix} test_prefix)
          math(EXPR UNIT_TEST_NUMBER "${UNIT_TEST_NUMBER}+1")
          set(UNIT_TEST${UNIT_TEST_NUMBER} ${test_prefix}.${filter_name} ${working_dir} ${executable} --gtest_filter=${filter_name}* CACHE STRING "" FORCE)
        endforeach()
    endforeach()
    set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
endfunction()

if(NOT IFEM_TESTING_INCLUDED)
  include(IFEMTesting)
endif()
