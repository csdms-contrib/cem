cmake_minimum_required(VERSION 2.6)

project (deltas)
include (CTest)
include(CheckIncludeFile)

set(CMAKE_MACOSX_RPATH 1)

#set (CMAKE_C_FLAGS "-Wall -Wundef -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -O2")
set (CMAKE_C_FLAGS "-Wall -O2")
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


########### libdeltas ###############
check_include_file(unistd.h WITH_UNISTD)

if (WITH_UNISTD)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWITH_UNISTD")
endif (WITH_UNISTD)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set( deltas_lib_SRCS
  bmi_cem.c
  cem_model.c
  ndelta4.c)
set_source_files_properties (${deltas_lib_SRCS} PROPERTIES LANGUAGE C)

find_library(LIB_M m)

if (LIB_M)
  set(MATH_LIB "m")
elseif (LIB_M)
  set(MATH_LIB "")
endif (LIB_M)

add_library(bmi_cem ${deltas_lib_SRCS})
add_library(bmi_cem-static STATIC ${deltas_lib_SRCS})
set_target_properties(bmi_cem-static PROPERTIES OUTPUT_NAME bmi_cem)

target_link_libraries(bmi_cem ${MATH_LIB} ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

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
set_target_properties(bmi_waves-static PROPERTIES OUTPUT_NAME bmi_waves)

target_link_libraries(bmi_waves ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

install(TARGETS bmi_waves DESTINATION lib)
install(TARGETS bmi_waves-static DESTINATION lib)
install(FILES bmi.h bmi_waves.h DESTINATION include/waves)

########### deltas main program ###############

set( deltas_SRCS deltas_main.c)
set_source_files_properties (${deltas_SRCS} PROPERTIES LANGUAGE C)

add_executable(run_deltas ${deltas_SRCS})
target_link_libraries(run_deltas bmi_cem-static ${MATH_LIB})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_deltas${CMAKE_EXECUTABLE_SUFFIX}
        DESTINATION bin
        RENAME deltas${CMAKE_EXECUTABLE_SUFFIX}
        COMPONENT deltas)

########### waves main program ###############

set (waves_SRCS waves_main.c waves_cli.c)
set_source_files_properties (${waves_SRCS} PROPERTIES LANGUAGE C)

add_executable(run_waves ${waves_SRCS})
target_link_libraries(run_waves bmi_waves-static ${GLIB2_LIBRARIES}
  ${GTHREAD2_LIBRARIES} ${MATH_LIB})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_waves${CMAKE_EXECUTABLE_SUFFIX}
        DESTINATION bin
        RENAME waves${CMAKE_EXECUTABLE_SUFFIX}
        COMPONENT waves)

########### cem main program ###############

set (cem_SRCS cem_main.c)
set_source_files_properties (${cem_SRCS} PROPERTIES LANGUAGE C)

add_executable(run_cem ${cem_SRCS})
target_link_libraries(run_cem bmi_waves-static bmi_cem-static
  ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES} ${MATH_LIB})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_cem${CMAKE_EXECUTABLE_SUFFIX}
        DESTINATION bin
        RENAME cem${CMAKE_EXECUTABLE_SUFFIX}
        COMPONENT cem)

########### Install files ###############

install(FILES deltas_api.h deltas_cli.h
        DESTINATION include
        COMPONENT deltas)

install (FILES deltas_api.h RENAME bmi_cem.h
        DESTINATION include
        COMPONENT deltas)

install(FILES deltas.pc DESTINATION lib/pkgconfig  COMPONENT deltas)
install(FILES waves.pc DESTINATION lib/pkgconfig  COMPONENT deltas)


