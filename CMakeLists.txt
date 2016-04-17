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


# Ecole Polytechnique Federale de Lausanne
# Brain Mind Institute,
# Blue Brain Project
# (c) 2006-2017. All rights reserved.
#
# Author: Bruno Magalhaes
# NeuroX


# Initial Setup

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(neurox)

set(NEUROX_DESCRIPTION "BBP NeuroX simulator")

set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "1")
set(VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
add_definitions("-DVERSION_STRING=\"${VERSION_STRING}\"")

## global options
option(ENABLE_SOA "Enable SoA Memory Layout" ON)
option(EXPORT_MECHS_FUNCTIONS "Enable exporting of mod init, jacob, state, cur and BeforeAfter functions of mod-based mechanisms" ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
      ${PROJECT_SOURCE_DIR}/CMake/portability
      ${PROJECT_SOURCE_DIR}/CMake/packages
      ${PROJECT_SOURCE_DIR}/CMake/config)

include(BlueGenePortability)
include(DefineInstallationPaths)
include(ReleaseDebugAutoFlags)
include(SetRpath)

include(FindPkgConfig)
include(CTest)

#Memory layout transformation: 0 for SoA and 1 for AoS
if (ENABLE_SOA)
    add_definitions("-DLAYOUT=0")
else()
    add_definitions("-DLAYOUT=1")
endif()

find_package(ClangFormat)

if(CLANG_FORMAT_FOUND)
    message("clang-format : ${CLANG_FORMAT_EXECUTABLE} : ${CLANG_FORMAT_VERSION}")
else()
    message("clang-format executable not found")
endif()

if(CLANG_FORMAT_FOUND)
    file(COPY ${PROJECT_SOURCE_DIR}/.clang-format DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    #to indent files using clang-format
    file(GLOB_RECURSE SRC_FILES_FOR_CLANG_FORMAT ${PROJECT_SOURCE_DIR}/neurox/*.c*
        ${PROJECT_SOURCE_DIR}/neurox/*.h* ${PROJECT_SOURCE_DIR}/neurox/*.ipp*)

    add_custom_target(formatsource COMMAND ${CMAKE_COMMAND}
        -DSOURCE_FILES:STRING="${SRC_FILES_FOR_CLANG_FORMAT}"
        -DCLANG_FORMAT_EXECUTABLE=${CLANG_FORMAT_EXECUTABLE}
        -P "${PROJECT_SOURCE_DIR}/CMake/ClangFormatUtils.cmake"
        )
endif()

#Turn off inline assembly for PGI
if(${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
    add_definitions(-DSWAP_ENDIAN_DISABLE_ASM)
endif()

if(EXPORT_MECHS_FUNCTIONS)
    add_definitions("-DEXPORT_MECHS_FUNCTIONS")
endif()

# Core library and application
include_directories(${CMAKE_INSTALL_PREFIX}/include) #coreneuron includes
link_directories(${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/lib64) #coreneuron libs
add_subdirectory(neurox)

# Documentation
#include(DoxygenRule)
