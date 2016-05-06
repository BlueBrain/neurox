# CMake generated Testfile for 
# Source directory: /home/bmagalha/Workspace/neurox/tests/unit/mechbuild
# Build directory: /home/bmagalha/Workspace/neurox/build/tests/unit/mechbuild
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(bacm-extra-no-nd "bash" "/home/bmagalha/Workspace/neurox/build/tests/unit/mechbuild/bacm-extra-no-nd.sh")
SET_TESTS_PROPERTIES(bacm-extra-no-nd PROPERTIES  LABELS "unit;modinclude" WORKING_DIRECTORY "/home/bmagalha/Workspace/neurox/build/tests/unit/mechbuild")
ADD_TEST(bacm-nd-extra2only "bash" "/home/bmagalha/Workspace/neurox/build/tests/unit/mechbuild/bacm-nd-extra2only.sh")
SET_TESTS_PROPERTIES(bacm-nd-extra2only PROPERTIES  LABELS "unit;modinclude" WORKING_DIRECTORY "/home/bmagalha/Workspace/neurox/build/tests/unit/mechbuild")
