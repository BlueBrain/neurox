set_source_files_properties(ring/ring.cpp ring_IClamp/ring_IClamp.cpp PROPERTIES COMPILE_FLAGS "-Wno-error")

# Includes the install_files function below
include(InstallFiles)

# TEST_LIBRARIES variable is used by the CommonCTest.cmake script to link against the given libraries
set(TEST_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} corebluron)

make_mpi_exeprefix(1 TEST_PREFIX)

# CommonCTest, in the current folder recursively compiles targets for *.cpp files using TEST_LIBRARIES
include(CommonCTest)

# This test is for corebluron havign all the mechanisms integrated into the
# library. Runs when the corresponding condition is utilized plus the
# executable is built.
if(COREBLURON_NEURODAMUS_MECHS AND COREBLURON_MAIN)

  make_mpi_exeprefix(2 MPI_PREFIX)
  string(REPLACE ";" " " SRUN_PREFIX "${MPI_PREFIX}")

  configure_file(10Cells/10cells_test.sh.in ${CMAKE_SOURCE_DIR}/tests/integration/x86/10Cells/10cells_test.sh @ONLY)

  # Installs the compiled tests ${TEST_FILES} into the examples directory
  set(TEST_FILES
  10Cells/10cells_test.sh
  10Cells/14397_1.dat
  10Cells/14397_2.dat
  10Cells/14457_1.dat
  10Cells/14457_2.dat
  10Cells/bbcore_mech.dat
  10Cells/files.dat
  10Cells/inputs.dat
  10Cells/out.dat.ref
  )

  add_test(NAME 10CELL_TEST
            COMMAND "/bin/sh" ${CMAKE_CURRENT_SOURCE_DIR}/10Cells/10cells_test.sh
            WORKING_DIRECTORY "${${CMAKE_PROJECT_NAME}_BINARY_DIR}") 
endif()

install_files(share/CoreBluron/tests FILES ${TEST_FILES} COMPONENT examples)
