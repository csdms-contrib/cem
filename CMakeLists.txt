cmake_minimum_required(VERSION 2.6)

project (cem)

set (CEM_VERSION 0.1)

set (BUILD_SHARED_LIBS ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set (bmicem_sources bmi.c)

add_library (bmicem ${bmicem_sources} )

install(TARGETS bmicem DESTINATION lib COMPONENT cem)

########### Install files ###############

install(FILES bmi.h bmi_const.h
        DESTINATION include/cem
        COMPONENT cem)
