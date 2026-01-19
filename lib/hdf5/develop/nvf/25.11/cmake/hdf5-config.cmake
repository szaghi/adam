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

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../" ABSOLUTE)

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

include(CMakeFindDependencyMacro)

string(TOUPPER hdf5 HDF5_PACKAGE_NAME)

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

set (${HDF5_PACKAGE_NAME}_BUILD_MODE       Release)

#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
# Languages:
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN     ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_CPP_LIB     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA        OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JNI         OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_JAVA_COMPAT      ON) # Always provide the Java compatibility API
set (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN "SHARED")
#-----------------------------------------------------------------------------
# Features:
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_HL_LIB          ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SHARED_LIBS     NO)
set (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS     YES)
set (${HDF5_PACKAGE_NAME}_PROVIDES_THREADS         ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_THREADSAFE      OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL        ON)
set (${HDF5_PACKAGE_NAME}_DEFAULT_API_VERSION     "v200")
set (${HDF5_PACKAGE_NAME}_PROVIDES_DEPRECATED_SYMBOLS ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_DIMENSION_SCALES_WITH_NEW_REF )
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_TOOLS           ON)
set (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_TOOLS    ON)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_NONSTANDARD_FEATURE_FLOAT16 ON)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT   OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT   OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_ENCODING  )
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PROVIDES_MAP_API        OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_DIRECT_VFD     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_MIRROR_VFD     OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_ROS3_VFD       OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_HDFS_VFD       OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_SUBFILING_VFD  OFF)
set (${HDF5_PACKAGE_NAME}_PROVIDES_PLUGIN_SUPPORT OFF)
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_PACKAGE_EXTLIBS       OFF)
set (${HDF5_PACKAGE_NAME}_EXPORT_LIBRARIES      hdf5-static;hdf5_tools-static;hdf5_hl-static;hdf5_f90cstub-static;hdf5_fortran-static;hdf5_hl_f90cstub-static;hdf5_hl_fortran-static)
set (${HDF5_PACKAGE_NAME}_ARCHITECTURE         "")
set (${HDF5_PACKAGE_NAME}_TOOLSET              "")

#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------
if (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL)
  set (${HDF5_PACKAGE_NAME}_MPI_C_INCLUDE_PATH "")
  set (${HDF5_PACKAGE_NAME}_MPI_C_LIBRARIES    "")
  if (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN)
    set (${HDF5_PACKAGE_NAME}_MPI_Fortran_INCLUDE_PATH "")
    set (${HDF5_PACKAGE_NAME}_MPI_Fortran_LIBRARIES    "")
  endif ()
  set (${HDF5_PACKAGE_NAME}_PROVIDES_PARALLEL_FILTERED_WRITES ON)
  set (${HDF5_PACKAGE_NAME}_PROVIDES_LARGE_PARALLEL_IO        ON)

  enable_language(C) # for MPI::MPI_C
  find_dependency(MPI QUIET REQUIRED)
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_THREADS)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_dependency(Threads QUIET REQUIRED)
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

if (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT)
  if (NOT OFF)
    if ()
      # Expect that the default shared library is expected with FindZLIB.cmake
      find_dependency (ZLIB MODULE)
    else ()
      find_dependency (ZLIB NAMES  CONFIG OPTIONAL_COMPONENTS )
    endif ()
  endif ()
endif ()

if (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT)
  if (NOT OFF)
    find_dependency (LIBAEC NAMES libaec CONFIG OPTIONAL_COMPONENTS )
  endif ()
endif ()

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include" "${${HDF5_PACKAGE_NAME}_MPI_C_INCLUDE_DIRS}")

set (${HDF5_PACKAGE_NAME}_SHARE_DIR "${PACKAGE_PREFIX_DIR}/cmake")
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

#-----------------------------------------------------------------------------
# Version Strings
#-----------------------------------------------------------------------------
set (${HDF5_PACKAGE_NAME}_VERSION 2.0.1)
set (${HDF5_PACKAGE_NAME}_VERSION_STRING 2.0.1)
set (${HDF5_PACKAGE_NAME}_VERSION_MAJOR  2.0)
set (${HDF5_PACKAGE_NAME}_VERSION_MINOR  1)

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built hdf5 as a subproject
#-----------------------------------------------------------------------------
if (NOT TARGET "hdf5")
  if (${HDF5_PACKAGE_NAME}_PROVIDES_ZLIB_SUPPORT)
    if (${HDF5_PACKAGE_NAME}_PACKAGE_EXTLIBS AND OFF)
      include (${PACKAGE_PREFIX_DIR}/cmake/zlib-targets.cmake)
    endif ()
  endif ()
  if (${HDF5_PACKAGE_NAME}_PROVIDES_SZIP_SUPPORT)
    if (${HDF5_PACKAGE_NAME}_PACKAGE_EXTLIBS AND OFF)
      include (${PACKAGE_PREFIX_DIR}/cmake/libaec-targets.cmake)
    endif ()
  endif ()
  include (${PACKAGE_PREFIX_DIR}/cmake/hdf5-targets.cmake)
endif ()

# Handle default component(static) :
if (NOT ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  if (${HDF5_PACKAGE_NAME}_PROVIDES_STATIC_LIBS)
    set (${HDF5_PACKAGE_NAME}_LIB_TYPE)
    set (${HDF5_PACKAGE_NAME}_FIND_COMPONENTS C HL static)
    set (${HDF5_PACKAGE_NAME}_FIND_REQUIRED_static_C true)
  else ()
    set (${HDF5_PACKAGE_NAME}_LIB_TYPE)
    set (${HDF5_PACKAGE_NAME}_FIND_COMPONENTS C HL shared)
    set (${HDF5_PACKAGE_NAME}_FIND_REQUIRED_shared_C true)
  endif ()
endif ()

# Handle requested components:
list (REMOVE_DUPLICATES ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
foreach (comp IN LISTS ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS)
  if (comp STREQUAL "shared")
    list (REMOVE_ITEM ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS ${comp})
    set (${HDF5_PACKAGE_NAME}_LIB_TYPE ${${HDF5_PACKAGE_NAME}_LIB_TYPE} ${comp})

    if (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN)
      if (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN STREQUAL "SHARED")
        set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/include")
      else ()
        set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/mod/shared")
      endif ()
    endif ()
  elseif (comp STREQUAL "static")
    list (REMOVE_ITEM ${HDF5_PACKAGE_NAME}_FIND_COMPONENTS ${comp})
    set (${HDF5_PACKAGE_NAME}_LIB_TYPE ${${HDF5_PACKAGE_NAME}_LIB_TYPE} ${comp})

    if (${HDF5_PACKAGE_NAME}_PROVIDES_FORTRAN)
      if (${HDF5_PACKAGE_NAME}_INSTALL_MOD_FORTRAN STREQUAL "STATIC")
        set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/include")
      else ()
        set (${HDF5_PACKAGE_NAME}_INCLUDE_DIR_FORTRAN "${PACKAGE_PREFIX_DIR}/mod/static")
      endif ()
    endif ()
  endif ()
endforeach ()
foreach (libtype IN LISTS ${HDF5_PACKAGE_NAME}_LIB_TYPE)
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
endforeach ()

foreach (libtype IN LISTS ${HDF5_PACKAGE_NAME}_LIB_TYPE)
  check_required_components(${HDF5_PACKAGE_NAME}_${libtype})
endforeach ()
