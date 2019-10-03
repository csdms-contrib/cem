cmake_minimum_required(VERSION 2.6)

project (cem)

set (CEM_VERSION 0.1)

# set (CMAKE_C_FLAGS "-Wall -Wundef -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -O2")
set (CMAKE_C_FLAGS "-Wall -O2")

########### CEM BMI library ###############

set (BUILD_SHARED_LIBS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(bmicem_sources cem/cem.c cem/rocks.c cem/conditions.c cem/utils.c bmi/bmi_cem.c)

add_library(bmi_cem ${bmicem_sources})
add_library(bmi_cem-static STATIC ${bmicem_sources})
set_target_properties(bmi_cem-static PROPERTIES OUTPUT_NAME bmi_cem)

# target_link_libraries(bmi_cem m)

########### CEM main ###############

set (cem_sources main.c)
add_executable (cem ${cem_sources})
# target_link_libraries(cem bmi_cem-static m)
target_link_libraries(cem bmi_cem-static)

########### Install files ###############
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cem${CMAKE_EXECUTABLE_SUFFIX} DESTINATION bin)
#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cem DESTINATION bin)
install(FILES bmi/bmi.h bmi/bmi_cem.h DESTINATION include/cem)
install(TARGETS bmi_cem DESTINATION lib)
install(TARGETS bmi_cem-static DESTINATION lib)

add_subdirectory(testing)
