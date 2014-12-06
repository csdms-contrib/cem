cmake_minimum_required(VERSION 2.6)

project (cem)

set (CEM_VERSION 0.1)

set (CMAKE_C_FLAGS
  "-ansi -Wall -Wundef -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -O2")

########### CEM BMI library ###############

set (BUILD_SHARED_LIBS ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set (bmicem_sources cem.c bmi.c)

add_library (bmicem ${bmicem_sources} )

install(TARGETS bmicem DESTINATION lib COMPONENT cem)

########### CEM main ###############

set (cem_sources bmi_main.c)
add_executable (cem ${cem_sources})
target_link_libraries(cem bmicem)

########### Install files ###############

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cem DESTINATION bin COMPONENT cem)
install(FILES bmi.h bmi_const.h
        DESTINATION include/cem
        COMPONENT cem)
