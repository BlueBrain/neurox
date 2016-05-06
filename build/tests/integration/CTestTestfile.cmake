# CMake generated Testfile for 
# Source directory: /home/bmagalha/Workspace/neurox/tests/integration
# Build directory: /home/bmagalha/Workspace/neurox/build/tests/integration
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(ring_TEST "/bin/sh" "/home/bmagalha/Workspace/neurox/build/tests/integration/ring/integration_test.sh")
SET_TESTS_PROPERTIES(ring_TEST PROPERTIES  WORKING_DIRECTORY "/home/bmagalha/Workspace/neurox/build/tests/integration/ring")
ADD_TEST(ring_IClamp_TEST "/bin/sh" "/home/bmagalha/Workspace/neurox/build/tests/integration/ring_IClamp/integration_test.sh")
SET_TESTS_PROPERTIES(ring_IClamp_TEST PROPERTIES  WORKING_DIRECTORY "/home/bmagalha/Workspace/neurox/build/tests/integration/ring_IClamp")
