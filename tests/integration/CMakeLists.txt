# Copyright (c) 2016, Blue Brain Project
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.


# List of tests with arguments
set(TEST_CASE_ARGS
  "ring!--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring --outpath=${CMAKE_CURRENT_BINARY_DIR}/ring --celsius=6.3 --tstop=100. -mpi"
  "ring_IClamp!--datpath=${CMAKE_CURRENT_SOURCE_DIR}/ring_IClamp --outpath=${CMAKE_CURRENT_BINARY_DIR}/ring_IClamp --celsius=6.3 --tstop=100. -mpi")

#Configure test scripts
foreach (args_line ${TEST_CASE_ARGS})
  string(REPLACE "!" ";" string_line ${args_line})
  string(REPLACE ";" " " SRUN_PREFIX  "${TEST_MPI_EXEC_BIN};-n;1")
  list(GET string_line 0 TEST_NAME)
  list(GET string_line 1 TEST_ARGS)
  set(SIM_NAME ${TEST_NAME})
  configure_file(integration_test.sh.in ${TEST_NAME}/integration_test.sh @ONLY)
  add_test(NAME ${TEST_NAME}_TEST
           COMMAND "/bin/sh" ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/integration_test.sh
           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
endforeach()
