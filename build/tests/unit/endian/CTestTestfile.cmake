# CMake generated Testfile for 
# Source directory: /home/bmagalha/Workspace/neurox/tests/unit/endian
# Build directory: /home/bmagalha/Workspace/neurox/build/tests/unit/endian
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(endianness_test.cpp_test "/usr/bin/srun" "/home/bmagalha/Workspace/neurox/build/tests/unit/endian/endianness_test.cpp_test_bin")
ADD_TEST(swap_endian_noasm.cpp_test "/usr/bin/srun" "/home/bmagalha/Workspace/neurox/build/tests/unit/endian/swap_endian_noasm.cpp_test_bin")
ADD_TEST(swap_endian_default.cpp_test "/usr/bin/srun" "/home/bmagalha/Workspace/neurox/build/tests/unit/endian/swap_endian_default.cpp_test_bin")
ADD_TEST(swap_endian_nounroll.cpp_test "/usr/bin/srun" "/home/bmagalha/Workspace/neurox/build/tests/unit/endian/swap_endian_nounroll.cpp_test_bin")
ADD_TEST(swap_endian_oddunroll.cpp_test "/usr/bin/srun" "/home/bmagalha/Workspace/neurox/build/tests/unit/endian/swap_endian_oddunroll.cpp_test_bin")
