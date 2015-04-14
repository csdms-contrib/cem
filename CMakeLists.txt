cmake_minimum_required(VERSION 2.6)

project (deltas)
include (CTest)

set (CMAKE_C_FLAGS "-ansi -Wall -Wundef -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -O2")
set (BUILD_SHARED_LIBS ON)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/waves.pc.cmake ${CMAKE_CURRENT_SOURCE_DIR}/waves.pc )
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/deltas.pc.cmake ${CMAKE_CURRENT_SOURCE_DIR}/deltas.pc )

add_subdirectory( data )

########### Look for glib-2.0 ###############

if (DEFINED CMAKE_GLIB_DIR)
  set (ENV{PKG_CONFIG_PATH}
       "$ENV{PKG_CONFIG_PATH}:${CMAKE_GLIB_DIR}/lib/pkgconfig")
endif (DEFINED CMAKE_GLIB_DIR)

include (FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
pkg_check_modules (GTHREAD2 gthread-2.0)
include_directories (${GLIB2_INCLUDE_DIRS} ${GTHREAD2_INCLUDE_DIRS})
link_directories (${GLIB2_LIBRARY_DIRS} ${GTHREAD2_LIBRARY_DIRS})

########### enable/disable graphics  ###############

option( WITH_X11 "Enable X11 graphics" OFF)
if (WITH_X11)
  include(FindPkgConfig)
  pkg_check_modules( X11 x11 )
  find_library( OPENGL GL )
  find_library( NCURSES ncurses )

  if (X11_FOUND AND OPENGL AND NCURSES)
    include_directories( ${X11_INCLUDE_DIRS} )
    link_directories( ${X11_LIBRARY_DIRS} )
    set(OPENGL_LIBRARIES ncurses GL ${X11_LIBRARIES})
  else (X11_FOUND AND OPENGL AND NCURSES)
    set(WITH_X11 OFF)
    message( STATUS "Disabling X11" )
  endif (X11_FOUND AND OPENGL AND NCURSES)

endif (WITH_X11)

if (WITH_X11)
  add_definitions( -DWITH_OPENGL )
else (WITH_X11)
  add_definitions( -UWITH_OPENGL )
endif (WITH_X11)

########### libdeltas ###############

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set( deltas_lib_SRCS
  bmi_cem.c
  deltas_cli.c
  ndelta4.c
  deltas_api.c)
set_source_files_properties (${deltas_lib_SRCS} PROPERTIES LANGUAGE C)

add_library(bmicem ${deltas_lib_SRCS})
add_library(bmicem-static STATIC ${deltas_lib_SRCS})
target_link_libraries (bmicem m ${OPENGL_LIBRARIES})
target_link_libraries(bmicem-static m ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

install(TARGETS bmicem DESTINATION lib COMPONENT deltas)

add_library(bmi_cem ${deltas_lib_SRCS})
add_library(bmi_cem-static STATIC ${deltas_lib_SRCS})
install(TARGETS bmi_cem DESTINATION lib)
install(TARGETS bmi_cem-static DESTINATION lib)
install(FILES bmi.h bmi_cem.h DESTINATION include/cem)

########### libwaves ###############

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

set (waves_lib_SRCS
  bmi_waves.c
  waves_model.c
  waves.c)
set_source_files_properties (${waves_lib_SRCS} PROPERTIES LANGUAGE C)

add_library(bmi_waves ${waves_lib_SRCS})
add_library(bmi_waves-static STATIC ${waves_lib_SRCS})

target_link_libraries(bmi_waves ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

install(TARGETS bmi_waves DESTINATION lib)
install(TARGETS bmi_waves-static DESTINATION lib)
install(FILES bmi.h bmi_waves.h DESTINATION include/waves)

########### deltas main program ###############

set( deltas_SRCS deltas_main.c)
set_source_files_properties (${deltas_SRCS} PROPERTIES LANGUAGE C)

add_executable(run_deltas ${deltas_SRCS})
target_link_libraries(run_deltas bmicem-static)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_deltas
        DESTINATION bin
        RENAME deltas
        COMPONENT deltas)

########### waves main program ###############

set (waves_SRCS waves_main.c waves_cli.c)
set_source_files_properties (${waves_SRCS} PROPERTIES LANGUAGE C)

add_executable(run_waves ${waves_SRCS})
target_link_libraries(run_waves bmi_waves-static ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_waves
        DESTINATION bin
        RENAME waves
        COMPONENT waves)

########### Python bindings ###############

option (PYTHON_BINDINGS "Use SWIG to build Python bindings" OFF)

if (PYTHON_BINDINGS )

# You may need to use -DCMAKE_INCLUDE_PATH and -DCMAKE_LIBRARY_PATH
# to find the correct libraries/include paths
find_package (SWIG REQUIRED)
include (${SWIG_USE_FILE})

find_package (PythonLibs 2.6 EXACT)
include_directories (${PYTHON_INCLUDE_PATH})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})

set (CMAKE_SWIG_FLAGS "")

#set_source_files_properties (deltas_mod.i PROPERTIES CPLUSPLUS ON) 
#set_source_files_properties (deltas_mod.i PROPERTIES SWIG_FLAGS "-includeall")
swig_add_module (deltas_mod python deltas_mod.i ndelta4.c deltas_api.c)
swig_link_libraries (deltas_mod ${PYTHON_LIBRARIES})

#add_library( _deltas_mod deltas_mod_wrap.c ndelta4.c deltas_api.c )
#target_link_libraries( deltas_mod _deltas_mod)

message ("Python include path " ${PYTHON_INCLUDE_PATH})
message ("Python libs " ${PYTHON_LIBRARIES})
message ("Python lib path " ${PYTHON_LIBRARY_PATH})

install (FILES ${CMAKE_BINARY_DIR}/deltas_mod.py
               ${CMAKE_BINARY_DIR}/_deltas_mod.so
         DESTINATION lib/deltas-python COMPONENT deltas)
endif (PYTHON_BINDINGS )

########### Install files ###############

install(FILES deltas_api.h deltas_cli.h
        DESTINATION include
        COMPONENT deltas)

install (FILES deltas_api.h RENAME bmi_cem.h
        DESTINATION include
        COMPONENT deltas)

install(FILES deltas.pc DESTINATION lib/pkgconfig  COMPONENT deltas)
install(FILES waves.pc DESTINATION lib/pkgconfig  COMPONENT deltas)


