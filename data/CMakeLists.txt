
set( DELTAS_EXE ${CMAKE_CURRENT_BINARY_DIR}/../run_deltas )

########### Testing ###############

ENABLE_TESTING()

add_test(DELTAS_VERSION ${DELTAS_EXE} --version )
add_test(DELTAS_HELP ${DELTAS_EXE} --help )
add_test(DELTAS_TEST ${DELTAS_EXE} --stop-time=10 --out-prefix=output )
add_test(DELTAS_DIFF diff output.50 ${CMAKE_CURRENT_SOURCE_DIR}/output/output.50 )

install(DIRECTORY input
        DESTINATION share/deltas
        COMPONENT deltas
        PATTERN ".svn" EXCLUDE)
install(DIRECTORY output
        DESTINATION share/deltas
        COMPONENT deltas
        PATTERN "output/output*"
        PATTERN ".svn" EXCLUDE)

