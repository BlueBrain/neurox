

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
