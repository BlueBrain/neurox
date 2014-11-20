set_source_files_properties(ring/ring.cpp ring_IClamp/ring_IClamp.cpp PROPERTIES COMPILE_FLAGS "-Wno-error")

# Includes the install_files function below
include(InstallFiles)

# TEST_LIBRARIES variable is used by the CommonCTest.cmake script to link against the given libraries
set(TEST_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} corebluron)

# CommonCTest, in the current folder recursively compiles targets for *.cpp files using TEST_LIBRARIES
#include(CommonCTest)


make_mpi_exeprefix(1 MPI_PREFIX)
string(REPLACE ";" " " SRUN_PREFIX "${MPI_PREFIX}")

#set(ring_ring_TEST_ARGS "--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring --celsius=6.3 -mpi")
#set(ring_IClamp_ring_IClamp_TEST_ARGS "--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring_IClamp --celsius=6.3 -mpi")

# List of tests with arguments
set(TEST_CASE_ARGS
  "ring!--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring --outpath=${CMAKE_CURRENT_BINARY_DIR}/ring --celsius=6.3 --tstop=100. -mpi"
  "ring_IClamp!--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring_IClamp --outpath=${CMAKE_CURRENT_BINARY_DIR}/ring_IClamp --celsius=6.3 --tstop=100. -mpi")

#Configure test scripts
foreach (args_line ${TEST_CASE_ARGS})
  string(REPLACE "!" ";" string_line ${args_line})
  list(GET string_line 0 TEST_NAME)
  list(GET string_line 1 TEST_ARGS)
  set(SIM_NAME ${TEST_NAME})
  configure_file(integration_test.sh.in ${TEST_NAME}/integration_test.sh @ONLY)
  add_test(NAME ${TEST_NAME}_TEST
           COMMAND "/bin/sh" ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/integration_test.sh
           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
endforeach()

# This test is for corebluron havign all the mechanisms integrated into the
# library. Runs when the corresponding condition is utilized plus the
# executable is built.
# Set of 10Cells tests
if(COREBLURON_NEURODAMUS_MECHS AND COREBLURON_MAIN)

  make_mpi_exeprefix(2 MPI_PREFIX)
  string(REPLACE ";" " " SRUN_PREFIX "${MPI_PREFIX}")

  function(exec_test TEST_NAME SUBTEST_NAME)
    set(SIM_NAME "${TEST_NAME}/${SUBTEST_NAME}")
    if(${SUBTEST_NAME} MATCHES "generic")
      set(TEST_ARGS "--datpath=${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME} --outpath=${CMAKE_CURRENT_BINARY_DIR}/${SIM_NAME} --tstop=100. -mpi")
    elseif(${SUBTEST_NAME} MATCHES "forwardskip")
      set(TEST_ARGS "--datpath=${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME} --outpath=${CMAKE_CURRENT_BINARY_DIR}/${SIM_NAME} --forwardskip=5000. --tstop=100. -mpi")
    elseif(${SUBTEST_NAME} MATCHES "spikereplay")
      set(TEST_ARGS "--datpath=${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME} --outpath=${CMAKE_CURRENT_BINARY_DIR}/${SIM_NAME} --pattern=${CMAKE_CURRENT_SOURCE_DIR}/${SIM_NAME}/out.std --tstop=100. -mpi")
    endif()

    configure_file(integration_test.sh.in ${SIM_NAME}/integration_test.sh @ONLY)

    add_test(NAME ${TEST_NAME}_${SUBTEST_NAME}_TEST
             COMMAND "/bin/sh" ${CMAKE_CURRENT_BINARY_DIR}/${SIM_NAME}/integration_test.sh
             WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SIM_NAME}")
  endfunction(exec_test)
  
  set(TEST_NAME "10Cells")
  set(SUBTEST_NAME "generic" "forwardskip" "spikereplay")
  foreach(subtest ${SUBTEST_NAME})
    exec_test(${TEST_NAME} ${subtest})
  endforeach()

  # Installs the compiled tests ${TEST_FILES} into the examples directory
#  set(TEST_FILES
#  10Cells/14397_1.dat
#  10Cells/14397_2.dat
#  10Cells/14457_1.dat
#  10Cells/14457_2.dat
#  10Cells/bbcore_mech.dat
#  10Cells/files.dat
#  )

endif()

#install_files(share/CoreBluron/tests FILES ${TEST_FILES} COMPONENT examples)
