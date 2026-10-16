#
# Copyright by The HDF Group.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the LICENSE file, which can be found at the root of the source code
# distribution tree, or in https://www.hdfgroup.org/licenses.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#
#-----------------------------------------------------------------------------
# HDF5 Config file for compiling against hdf5 build/install directory
# Sections that depend on an ENABLE type option should first check if the
# option is enabled before using the dependent variables.
#-----------------------------------------------------------------------------

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was hdf5-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include (CMakeFindDependencyMacro)

string (TOUPPER hdf5 HDF5_PACKAGE_NAME)

set (${HDF5_PACKAGE_NAME}_VALID_COMPONENTS
    static
    shared
    C
    CXX
    Fortran
    HL
    CXX_HL
    Fortran_HL
    Java
    Tools
)

set (${HDF5_PACKAGE_NAME}_BUILD_MODE Release)

#-----------------------------------------------------------------------------
# Version Strings
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_VERSION        2.3.0)
set (${HDF5_PACKAGE_NAME}_VERSION_STRING 2.3.0)
set (${HDF5_PACKAGE_NAME}_VERSION_MAJOR  2.3)
set (${HDF5_PACKAGE_NAME}_VERSION_MINOR  0)

#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
# Languages:
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN     ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_CPP_LIB     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA        OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JNI         OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA_COMPAT ON) # Always provide the Java compatibility API
set (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN  "SHARED")
#-----------------------------------------------------------------------------
# Features:
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_SHARED_LIBS                   NO)
set (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS                   YES)
set (${HDF5_PACKAGE_NAME}_PROVIDES_HL_LIB                        ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_THREADS                       ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_THREADSAFE                    OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL                      ON)
set (${HDF5_PACKAGE_NAME}_DEFAULT_API_VERSION                    "v200")
set (${HDF5_PACKAGE_NAME}_PROVIDES_DEPRECATED_SYMBOLS            ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_DIMENSION_SCALES_WITH_NEW_REF OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_PLUGIN_SIGNATURE_VERIFICATION OFF)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_TOOLS        ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_TOOLS ON)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_NONSTANDARD_FEATURE_FLOAT16 ON)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT  OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT  OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_ENCODING OFF)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_MAP_API        OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_DIRECT_VFD     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_MIRROR_VFD     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_ROS3_VFD       OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_HDFS_VFD       OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SUBFILING_VFD  OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_PLUGIN_SUPPORT OFF)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PACKAGE_EXTLIBS  OFF)
set (${HDF5_PACKAGE_NAME}_EXPORT_LIBRARIES hdf5-static;hdf5_tools-static;hdf5_hl-static;hdf5_f90cstub-static;hdf5_fortran-static;hdf5_hl_f90cstub-static;hdf5_hl_fortran-static)
set (${HDF5_PACKAGE_NAME}_ARCHITECTURE     "")
set (${HDF5_PACKAGE_NAME}_TOOLSET          "")

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include" "${${HDF5_PACKAGE_NAME}_MPI_C_INCLUDE_DIRS}")

set (${HDF5_PACKAGE_NAME}_SHARE_DIR "${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5")
set_and_check (${HDF5_PACKAGE_NAME}_BUILD_DIR "${PACKAGE_PREFIX_DIR}")

if (${HDF5_PACKAGE_NAME}_PROVIDES_CPP_LIB)
  set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_CPP "${PACKAGE_PREFIX_DIR}/include")
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_HL_LIB)
  set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_HL "${PACKAGE_PREFIX_DIR}/include")
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_HL_LIB AND ${HDF5_PACKAGE_NAME}_PROVIDES_CPP_LIB)
  set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_HL_CPP "${PACKAGE_PREFIX_DIR}/include")
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_TOOLS)
  set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_TOOLS "${PACKAGE_PREFIX_DIR}/include")
  set_and_check (${HDF5_PACKAGE_NAME}_TOOLS_DIR "${PACKAGE_PREFIX_DIR}/bin")
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA AND ${HDF5_PACKAGE_NAME}_PROVIDES_JNI)
  set (${HDF5_PACKAGE_NAME}_JAVA_INCLUDE_DIRS
      ${PACKAGE_PREFIX_DIR}/lib/
      ${PACKAGE_PREFIX_DIR}/lib/
      ${PACKAGE_PREFIX_DIR}/lib/
  )
else ()
  set (${HDF5_PACKAGE_NAME}_JAVA_INCLUDE_DIRS
      ${PACKAGE_PREFIX_DIR}/lib/
      ${PACKAGE_PREFIX_DIR}/lib/
      ${PACKAGE_PREFIX_DIR}/lib/
  )
  if (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA_COMPAT)
    list (APPEND ${HDF5_PACKAGE_NAME}_JAVA_INCLUDE_DIRS
        ${PACKAGE_PREFIX_DIR}/lib/
    )
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Handle default components if none specified
#-----------------------------------------------------------------------------
if (NOT ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  set (${HDF5_PACKAGE_NAME}_FIND_COMPONENTS C HL)
endif ()

#-----------------------------------------------------------------------------
# Handle default shared/static libraries if none specified; prefer static
# libraries if available
#-----------------------------------------------------------------------------
if (NOT "shared" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS AND NOT "static" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  if (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS)
    list (APPEND ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS static)
    set (${HDF5_PACKAGE_NAME}_FIND_REQUIRED_static_C true)
  else ()
    list (APPEND ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS shared)
    set (${HDF5_PACKAGE_NAME}_FIND_REQUIRED_shared_C true)
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Resolve request for shared vs. static libraries before resolving
# dependencies
#-----------------------------------------------------------------------------
unset (${HDF5_PACKAGE_NAME}_LIB_TYPE)
list (REMOVE_DUPLICATES ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
if ("shared" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS AND "static" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
  set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "The 'shared' and 'static' components are mutually exclusive; select only one.")
  return ()
endif ()
if ("shared" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  if (NOT ${HDF5_PACKAGE_NAME}_PROVIDES_SHARED_LIBS)
    set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
    set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "Shared HDF5 libraries are not available.")
    return ()
  endif ()

  set (${HDF5_PACKAGE_NAME}_LIB_TYPE "shared")
  if (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN)
    if (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN STREQUAL "SHARED")
      set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/include")
    else ()
      set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/mod/shared")
    endif ()
  endif ()
elseif ("static" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  if (NOT ${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS)
    set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
    set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "Static HDF5 libraries are not available.")
    return ()
  endif ()

  set (${HDF5_PACKAGE_NAME}_LIB_TYPE "static")

  # Static libraries on Windows are built without thread support
  set (${HDF5_PACKAGE_NAME}_PROVIDES_THREADS ON)
  if (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN)
    if (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN STREQUAL "STATIC")
      set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/include")
    else ()
      set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/mod/static")
    endif ()
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Resolve dependencies linked with PUBLIC scope
# - these should always be propagated with find_dependency() calls, unless
#   they are only used by a component that wasn't requested
#-----------------------------------------------------------------------------
if (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL)
  unset (_hdf5_mpi_components)

  get_property (_project_enabled_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  if ("C" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
    if (NOT "C" IN_LIST _project_enabled_languages)
      set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
      set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "C language was not enabled; enable C language with project() or enable_language() commands.")
      unset (_project_enabled_languages)
      return ()
    endif ()
    list (APPEND _hdf5_mpi_components "C")
  endif ()
  if ("Fortran" IN_LIST ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
    if (NOT "Fortran" IN_LIST _project_enabled_languages)
      set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
      set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "Fortran language was not enabled; enable Fortran language with project() or enable_language() commands.")
      unset (_project_enabled_languages)
      return ()
    endif ()
    list (APPEND _hdf5_mpi_components "Fortran")
  endif ()
  unset (_project_enabled_languages)

  find_dependency (MPI QUIET REQUIRED COMPONENTS ${_hdf5_mpi_components})

  if ("C" IN_LIST _hdf5_mpi_components)
    set (${HDF5_PACKAGE_NAME}_MPI_C_INCLUDE_PATH "")
    set (${HDF5_PACKAGE_NAME}_MPI_C_LIBRARIES    "")
  endif ()
  if ("Fortran" IN_LIST _hdf5_mpi_components)
    set (${HDF5_PACKAGE_NAME}_MPI_Fortran_INCLUDE_PATH "")
    set (${HDF5_PACKAGE_NAME}_MPI_Fortran_LIBRARIES    "")
  endif ()

  set (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL_FILTERED_WRITES ON)
  set (${HDF5_PACKAGE_NAME}_PROVIDES_LARGE_PARALLEL_IO        ON)

  unset (_hdf5_mpi_components)
endif ()

#-----------------------------------------------------------------------------
# Resolve dependencies linked with PRIVATE scope
# - these should be propagated with find_dependency() calls when used by
#   requested components that are static libraries 
#-----------------------------------------------------------------------------
if (${HDF5_PACKAGE_NAME}_LIB_TYPE STREQUAL "static")
  # Threads (linked even where the library is built without thread support)
  if (ON)
    set (THREADS_PREFER_PTHREAD_FLAG ON)
    find_dependency (Threads QUIET REQUIRED)
  endif ()

  # zlib / zlib-ng
  if (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT AND NOT OFF)
    if (NOT OFF)
      if (OFF)
        # Expect that the default shared library is expected with FindZLIB.cmake
        find_dependency (ZLIB MODULE)
      else ()
        find_dependency (ZLIB NAMES  CONFIG OPTIONAL_COMPONENTS )
      endif ()
    else ()
      # Locate zlib-ng using installed FindZLIBNG.cmake module. If a CMake configuration
      # file is available for zlib-ng, the module will prefer to use it.
      set (_hdf5_cmake_module_path_saved ${CMAKE_MODULE_PATH})
      set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/Modules")
      find_dependency (ZLIBNG OPTIONAL_COMPONENTS )
      set (CMAKE_MODULE_PATH ${_hdf5_cmake_module_path_saved})
      unset (_hdf5_cmake_module_path_saved)
    endif ()
  endif ()

  # libaec - if szip/libaec support was enabled and libaec was found on the system at build
  # time, locate libaec to propagate the dependency 
  if (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT AND NOT OFF)
    # Locate libaec using installed Findlibaec.cmake module. If a CMake configuration
    # file is available for libaec, the module will prefer to use it.
    set (_hdf5_cmake_module_path_saved ${CMAKE_MODULE_PATH})
    set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/Modules")
    find_dependency (libaec OPTIONAL_COMPONENTS )
    set (CMAKE_MODULE_PATH ${_hdf5_cmake_module_path_saved})
    unset (_hdf5_cmake_module_path_saved)
  endif ()

  # aws-c-s3 library for ROS3 VFD
  if (${HDF5_PACKAGE_NAME}_PROVIDES_ROS3_VFD)
    find_dependency (aws-c-s3)
  endif ()

  # OpenSSL
  if (${HDF5_PACKAGE_NAME}_PROVIDES_PLUGIN_SIGNATURE_VERIFICATION)
    find_dependency (OpenSSL)
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built hdf5 as a subproject
#-----------------------------------------------------------------------------
if (NOT TARGET "hdf5")
  # If zlib support was enabled and zlib was built from an external source,
  # include the installed CMake targets file so the targets are available
  if (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT AND OFF)
    if (OFF)
      set (_zlib_targets_file_name "${PACKAGE_PREFIX_DIR}//-targets.cmake")
    else ()
      string (TOUPPER "" ${HDF5_PACKAGE_NAME}_ZLIB_INSTALL_NAME_UPPER)
      set (_zlib_targets_file_name "${PACKAGE_PREFIX_DIR}//${${HDF5_PACKAGE_NAME}_ZLIB_INSTALL_NAME_UPPER}Config.cmake")
      unset (${HDF5_PACKAGE_NAME}_ZLIB_INSTALL_NAME_UPPER)

      if (OFF)
        set (ZLIB_FIND_COMPONENTS "static")
      else ()
        set (ZLIB_FIND_COMPONENTS "shared")
      endif ()
    endif ()
    include ("${_zlib_targets_file_name}" OPTIONAL RESULT_VARIABLE _zlib_loaded)
    unset (ZLIB_FIND_COMPONENTS)
    if (NOT _zlib_loaded)
      if (OFF)
        set (_zlib_name "zlib-ng")
      else ()
        set (_zlib_name "zlib")
      endif ()
      set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
      set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "Couldn't locate ${_zlib_name} targets file ${_zlib_targets_file_name}")
      unset (_zlib_targets_file_name)
      unset (_zlib_name)
      unset (_zlib_loaded)
      return ()
    endif ()
    unset (_zlib_targets_file_name)
    unset (_zlib_loaded)
  endif ()
  # If szip/libaec support was enabled and libaec was built from an external source,
  # include the installed CMake targets file so the targets are available
  if (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT AND OFF)
    if (OFF)
      set (_libaec_targets_file_name "${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5/libaec_static-targets.cmake")
    else ()
      set (_libaec_targets_file_name "${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5/libaec_shared-targets.cmake")
    endif ()
    include ("${_libaec_targets_file_name}" OPTIONAL RESULT_VARIABLE _libaec_loaded)
    if (NOT _libaec_loaded)
      set (${HDF5_PACKAGE_NAME}_FOUND FALSE)
      set (${HDF5_PACKAGE_NAME}_NOT_FOUND_MESSAGE "Couldn't locate libaec targets file ${_libaec_targets_file_name}")
      unset (_libaec_targets_file_name)
      unset (_libaec_loaded)
      return ()
    endif ()
    unset (_libaec_targets_file_name)
    unset (_libaec_loaded)
  endif ()
  include ("${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5/hdf5-targets.cmake")
  # OPTIONAL: some packages ship HDF5 in such a way that the static and
  # Java targets files may not be available.
  if (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS AND ${HDF5_PACKAGE_NAME}_PROVIDES_SHARED_LIBS)
    include ("${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5/hdf5_static-targets.cmake" OPTIONAL RESULT_VARIABLE _hdf5_static_loaded)
    if (NOT _hdf5_static_loaded)
      set (${CMAKE_FIND_PACKAGE_NAME}_static_FOUND FALSE)
    endif ()
    unset (_hdf5_static_loaded)
  endif ()
  if (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA)
    include ("${PACKAGE_PREFIX_DIR}/lib/cmake/hdf5/hdf5_java-targets.cmake" OPTIONAL RESULT_VARIABLE _hdf5_java_loaded)
    if (NOT _hdf5_java_loaded)
      set (${CMAKE_FIND_PACKAGE_NAME}_Java_FOUND FALSE)
    endif ()
    unset (_hdf5_java_loaded)
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Handle all other requested components
#-----------------------------------------------------------------------------
set (libtype ${${HDF5_PACKAGE_NAME}_LIB_TYPE})
foreach (comp IN LISTS ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  set (hdf5_comp2)
  if (comp STREQUAL "C")
    set (hdf5_comp "hdf5")
  elseif (comp STREQUAL "CXX")
    set (hdf5_comp "hdf5_cpp")
  elseif (comp STREQUAL "HL")
    set (hdf5_comp "hdf5_hl")
  elseif (comp STREQUAL "CXX_HL")
    set (hdf5_comp "hdf5_hl_cpp")
  elseif (comp STREQUAL "Java")
    set (hdf5_comp "hdf5_java")
  elseif (comp STREQUAL "Tools")
    set (hdf5_comp "hdf5_tools")
  elseif (comp STREQUAL "Fortran")
    set (hdf5_comp2 "hdf5_f90cstub")
    set (hdf5_comp "hdf5_fortran")
  elseif (comp STREQUAL "Fortran_HL")
    set (hdf5_comp2 "hdf5_hl_f90cstub")
    set (hdf5_comp "hdf5_hl_fortran")
  endif ()
  if (comp STREQUAL "Java")
    list (FIND ${HDF5_PACKAGE_NAME}_EXPORT_LIBRARIES "${hdf5_comp}" HAVE_COMP)
  else ()
    list (FIND ${HDF5_PACKAGE_NAME}_EXPORT_LIBRARIES "${hdf5_comp}-${libtype}" HAVE_COMP)
  endif ()
  if (${HAVE_COMP} LESS 0)
    if (comp STREQUAL "Java")
      set (${HDF5_PACKAGE_NAME}_${comp}_FOUND 0)
    else ()
      set (${HDF5_PACKAGE_NAME}_${libtype}_${comp}_FOUND 0)
    endif ()
  else ()
    if (hdf5_comp2)
      list (FIND ${HDF5_PACKAGE_NAME}_EXPORT_LIBRARIES "${hdf5_comp2}-${libtype}" HAVE_COMP2)
      if (${HAVE_COMP2} LESS 0)
        set (${HDF5_PACKAGE_NAME}_${libtype}_${comp}_FOUND 0)
      else ()
        set (${HDF5_PACKAGE_NAME}_${libtype}_${comp}_FOUND 1)
        string(TOUPPER ${HDF5_PACKAGE_NAME}_${comp}_${libtype}_LIBRARY COMP_LIBRARY)
        set (${COMP_LIBRARY} ${${COMP_LIBRARY}} ${hdf5_comp2}-${libtype} ${hdf5_comp}-${libtype})
      endif ()
    elseif (comp STREQUAL "Java")
      set (${HDF5_PACKAGE_NAME}_${comp}_FOUND 1)
      string(TOUPPER ${HDF5_PACKAGE_NAME}_${comp}_LIBRARY COMP_LIBRARY)
      set (${COMP_LIBRARY} ${${COMP_LIBRARY}} ${hdf5_comp})
    else ()
      set (${HDF5_PACKAGE_NAME}_${libtype}_${comp}_FOUND 1)
      string(TOUPPER ${HDF5_PACKAGE_NAME}_${comp}_${libtype}_LIBRARY COMP_LIBRARY)
      set (${COMP_LIBRARY} ${${COMP_LIBRARY}} ${hdf5_comp}-${libtype})
    endif ()
  endif ()
endforeach ()
unset (libtype)

# Should be last
check_required_components (${HDF5_PACKAGE_NAME}_${${HDF5_PACKAGE_NAME}_LIB_TYPE})
