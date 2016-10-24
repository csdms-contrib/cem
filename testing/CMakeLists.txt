cmake_minimum_required(VERSION 2.6)

include(FindPkgConfig)
pkg_check_modules(GLIB2 glib-2.0)
pkg_check_modules(GTHREAD2 gthread-2.0)
include_directories(${GLIB2_INCLUDE_DIRS} ${GTHREAD2_INCLUDE_DIRS})
link_directories(${GLIB2_LIBRARY_DIRS} ${GTHREAD2_LIBRARY_DIRS})


set(test_irf_srcs test_irf.c)
add_executable(test-irf ${test_irf_srcs})
target_link_libraries(test-irf
  bmi_cem ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

set(test_meta_srcs test_meta.c)
add_executable(test-meta ${test_meta_srcs})
target_link_libraries(test-meta
  bmi_cem ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})

set(test_var_srcs test_var.c)
add_executable(test-var ${test_var_srcs})
target_link_libraries(test-var
  bmi_cem ${GLIB2_LIBRARIES} ${GTHREAD2_LIBRARIES})
