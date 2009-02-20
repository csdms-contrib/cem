cmake_minimum_required(VERSION 2.6)

project( deltas )

add_subdirectory( data )

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
  deltas_cli.c
  ndelta4.c
  deltas_api.c)

add_library( deltas ${deltas_lib_SRCS} )
target_link_libraries( deltas m ${OPENGL_LIBRARIES} )

install(TARGETS deltas DESTINATION lib COMPONENT deltas)


########### deltas main program ###############

set( deltas_SRCS deltas_main.c)

add_executable(run_deltas ${deltas_SRCS})
target_link_libraries(run_deltas deltas)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_deltas
        DESTINATION bin
        RENAME deltas
        COMPONENT deltas)

########### Install files ###############

install(FILES deltas_api.h deltas_cli.h
        DESTINATION include
        COMPONENT deltas)

