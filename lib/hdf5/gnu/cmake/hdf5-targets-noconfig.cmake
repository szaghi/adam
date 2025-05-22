#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hdf5-static" for configuration ""
set_property(TARGET hdf5-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5.a"
  )

list(APPEND _cmake_import_check_targets hdf5-static )
list(APPEND _cmake_import_check_files_for_hdf5-static "${_IMPORT_PREFIX}/lib/libhdf5.a" )

# Import target "hdf5-shared" for configuration ""
set_property(TARGET hdf5-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5.so.310.5.1"
  IMPORTED_SONAME_NOCONFIG "libhdf5.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5-shared )
list(APPEND _cmake_import_check_files_for_hdf5-shared "${_IMPORT_PREFIX}/lib/libhdf5.so.310.5.1" )

# Import target "mirror_server" for configuration ""
set_property(TARGET mirror_server APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(mirror_server PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/mirror_server"
  )

list(APPEND _cmake_import_check_targets mirror_server )
list(APPEND _cmake_import_check_files_for_mirror_server "${_IMPORT_PREFIX}/bin/mirror_server" )

# Import target "mirror_server_stop" for configuration ""
set_property(TARGET mirror_server_stop APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(mirror_server_stop PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/mirror_server_stop"
  )

list(APPEND _cmake_import_check_targets mirror_server_stop )
list(APPEND _cmake_import_check_files_for_mirror_server_stop "${_IMPORT_PREFIX}/bin/mirror_server_stop" )

# Import target "hdf5_tools-static" for configuration ""
set_property(TARGET hdf5_tools-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_tools-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_tools.a"
  )

list(APPEND _cmake_import_check_targets hdf5_tools-static )
list(APPEND _cmake_import_check_files_for_hdf5_tools-static "${_IMPORT_PREFIX}/lib/libhdf5_tools.a" )

# Import target "hdf5_tools-shared" for configuration ""
set_property(TARGET hdf5_tools-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_tools-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_tools.so.310.0.6"
  IMPORTED_SONAME_NOCONFIG "libhdf5_tools.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_tools-shared )
list(APPEND _cmake_import_check_files_for_hdf5_tools-shared "${_IMPORT_PREFIX}/lib/libhdf5_tools.so.310.0.6" )

# Import target "h5diff" for configuration ""
set_property(TARGET h5diff APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5diff PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5diff"
  )

list(APPEND _cmake_import_check_targets h5diff )
list(APPEND _cmake_import_check_files_for_h5diff "${_IMPORT_PREFIX}/bin/h5diff" )

# Import target "ph5diff" for configuration ""
set_property(TARGET ph5diff APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(ph5diff PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/ph5diff"
  )

list(APPEND _cmake_import_check_targets ph5diff )
list(APPEND _cmake_import_check_files_for_ph5diff "${_IMPORT_PREFIX}/bin/ph5diff" )

# Import target "h5ls" for configuration ""
set_property(TARGET h5ls APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5ls PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5ls"
  )

list(APPEND _cmake_import_check_targets h5ls )
list(APPEND _cmake_import_check_files_for_h5ls "${_IMPORT_PREFIX}/bin/h5ls" )

# Import target "h5debug" for configuration ""
set_property(TARGET h5debug APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5debug PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5debug"
  )

list(APPEND _cmake_import_check_targets h5debug )
list(APPEND _cmake_import_check_files_for_h5debug "${_IMPORT_PREFIX}/bin/h5debug" )

# Import target "h5repart" for configuration ""
set_property(TARGET h5repart APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5repart PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5repart"
  )

list(APPEND _cmake_import_check_targets h5repart )
list(APPEND _cmake_import_check_files_for_h5repart "${_IMPORT_PREFIX}/bin/h5repart" )

# Import target "h5mkgrp" for configuration ""
set_property(TARGET h5mkgrp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5mkgrp PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5mkgrp"
  )

list(APPEND _cmake_import_check_targets h5mkgrp )
list(APPEND _cmake_import_check_files_for_h5mkgrp "${_IMPORT_PREFIX}/bin/h5mkgrp" )

# Import target "h5clear" for configuration ""
set_property(TARGET h5clear APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5clear PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5clear"
  )

list(APPEND _cmake_import_check_targets h5clear )
list(APPEND _cmake_import_check_files_for_h5clear "${_IMPORT_PREFIX}/bin/h5clear" )

# Import target "h5delete" for configuration ""
set_property(TARGET h5delete APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5delete PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5delete"
  )

list(APPEND _cmake_import_check_targets h5delete )
list(APPEND _cmake_import_check_files_for_h5delete "${_IMPORT_PREFIX}/bin/h5delete" )

# Import target "h5import" for configuration ""
set_property(TARGET h5import APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5import PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5import"
  )

list(APPEND _cmake_import_check_targets h5import )
list(APPEND _cmake_import_check_files_for_h5import "${_IMPORT_PREFIX}/bin/h5import" )

# Import target "h5repack" for configuration ""
set_property(TARGET h5repack APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5repack PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5repack"
  )

list(APPEND _cmake_import_check_targets h5repack )
list(APPEND _cmake_import_check_files_for_h5repack "${_IMPORT_PREFIX}/bin/h5repack" )

# Import target "h5jam" for configuration ""
set_property(TARGET h5jam APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5jam PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5jam"
  )

list(APPEND _cmake_import_check_targets h5jam )
list(APPEND _cmake_import_check_files_for_h5jam "${_IMPORT_PREFIX}/bin/h5jam" )

# Import target "h5unjam" for configuration ""
set_property(TARGET h5unjam APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5unjam PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5unjam"
  )

list(APPEND _cmake_import_check_targets h5unjam )
list(APPEND _cmake_import_check_files_for_h5unjam "${_IMPORT_PREFIX}/bin/h5unjam" )

# Import target "h5copy" for configuration ""
set_property(TARGET h5copy APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5copy PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5copy"
  )

list(APPEND _cmake_import_check_targets h5copy )
list(APPEND _cmake_import_check_files_for_h5copy "${_IMPORT_PREFIX}/bin/h5copy" )

# Import target "h5stat" for configuration ""
set_property(TARGET h5stat APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5stat PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5stat"
  )

list(APPEND _cmake_import_check_targets h5stat )
list(APPEND _cmake_import_check_files_for_h5stat "${_IMPORT_PREFIX}/bin/h5stat" )

# Import target "h5dump" for configuration ""
set_property(TARGET h5dump APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5dump PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5dump"
  )

list(APPEND _cmake_import_check_targets h5dump )
list(APPEND _cmake_import_check_files_for_h5dump "${_IMPORT_PREFIX}/bin/h5dump" )

# Import target "h5format_convert" for configuration ""
set_property(TARGET h5format_convert APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5format_convert PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5format_convert"
  )

list(APPEND _cmake_import_check_targets h5format_convert )
list(APPEND _cmake_import_check_files_for_h5format_convert "${_IMPORT_PREFIX}/bin/h5format_convert" )

# Import target "h5perf_serial" for configuration ""
set_property(TARGET h5perf_serial APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5perf_serial PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5perf_serial"
  )

list(APPEND _cmake_import_check_targets h5perf_serial )
list(APPEND _cmake_import_check_files_for_h5perf_serial "${_IMPORT_PREFIX}/bin/h5perf_serial" )

# Import target "h5perf" for configuration ""
set_property(TARGET h5perf APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5perf PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5perf"
  )

list(APPEND _cmake_import_check_targets h5perf )
list(APPEND _cmake_import_check_files_for_h5perf "${_IMPORT_PREFIX}/bin/h5perf" )

# Import target "hdf5_hl-static" for configuration ""
set_property(TARGET hdf5_hl-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl.a"
  )

list(APPEND _cmake_import_check_targets hdf5_hl-static )
list(APPEND _cmake_import_check_files_for_hdf5_hl-static "${_IMPORT_PREFIX}/lib/libhdf5_hl.a" )

# Import target "hdf5_hl-shared" for configuration ""
set_property(TARGET hdf5_hl-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl.so.310.0.6"
  IMPORTED_SONAME_NOCONFIG "libhdf5_hl.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_hl-shared )
list(APPEND _cmake_import_check_files_for_hdf5_hl-shared "${_IMPORT_PREFIX}/lib/libhdf5_hl.so.310.0.6" )

# Import target "h5watch" for configuration ""
set_property(TARGET h5watch APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5watch PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/h5watch"
  )

list(APPEND _cmake_import_check_targets h5watch )
list(APPEND _cmake_import_check_files_for_h5watch "${_IMPORT_PREFIX}/bin/h5watch" )

# Import target "hdf5_f90cstub-static" for configuration ""
set_property(TARGET hdf5_f90cstub-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_f90cstub-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.a"
  )

list(APPEND _cmake_import_check_targets hdf5_f90cstub-static )
list(APPEND _cmake_import_check_files_for_hdf5_f90cstub-static "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.a" )

# Import target "hdf5_f90cstub-shared" for configuration ""
set_property(TARGET hdf5_f90cstub-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_f90cstub-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.so.310.3.2"
  IMPORTED_SONAME_NOCONFIG "libhdf5_f90cstub.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_f90cstub-shared )
list(APPEND _cmake_import_check_files_for_hdf5_f90cstub-shared "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.so.310.3.2" )

# Import target "hdf5_fortran-static" for configuration ""
set_property(TARGET hdf5_fortran-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_fortran-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_fortran.a"
  )

list(APPEND _cmake_import_check_targets hdf5_fortran-static )
list(APPEND _cmake_import_check_files_for_hdf5_fortran-static "${_IMPORT_PREFIX}/lib/libhdf5_fortran.a" )

# Import target "hdf5_fortran-shared" for configuration ""
set_property(TARGET hdf5_fortran-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_fortran-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_fortran.so.310.3.2"
  IMPORTED_SONAME_NOCONFIG "libhdf5_fortran.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_fortran-shared )
list(APPEND _cmake_import_check_files_for_hdf5_fortran-shared "${_IMPORT_PREFIX}/lib/libhdf5_fortran.so.310.3.2" )

# Import target "hdf5_hl_f90cstub-static" for configuration ""
set_property(TARGET hdf5_hl_f90cstub-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl_f90cstub-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.a"
  )

list(APPEND _cmake_import_check_targets hdf5_hl_f90cstub-static )
list(APPEND _cmake_import_check_files_for_hdf5_hl_f90cstub-static "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.a" )

# Import target "hdf5_hl_f90cstub-shared" for configuration ""
set_property(TARGET hdf5_hl_f90cstub-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl_f90cstub-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.so.310.0.6"
  IMPORTED_SONAME_NOCONFIG "libhdf5_hl_f90cstub.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_hl_f90cstub-shared )
list(APPEND _cmake_import_check_files_for_hdf5_hl_f90cstub-shared "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.so.310.0.6" )

# Import target "hdf5_hl_fortran-static" for configuration ""
set_property(TARGET hdf5_hl_fortran-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl_fortran-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.a"
  )

list(APPEND _cmake_import_check_targets hdf5_hl_fortran-static )
list(APPEND _cmake_import_check_files_for_hdf5_hl_fortran-static "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.a" )

# Import target "hdf5_hl_fortran-shared" for configuration ""
set_property(TARGET hdf5_hl_fortran-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(hdf5_hl_fortran-shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.so.310.0.6"
  IMPORTED_SONAME_NOCONFIG "libhdf5_hl_fortran.so.310"
  )

list(APPEND _cmake_import_check_targets hdf5_hl_fortran-shared )
list(APPEND _cmake_import_check_files_for_hdf5_hl_fortran-shared "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.so.310.0.6" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
