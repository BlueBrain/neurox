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
# (c) 2006-2016. All rights reserved.
#
# Author: Aleksandr Ovcharenko
# Core Neuron


# Initial Setup

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(neurox)

set(CORENEURON_DESCRIPTION "BBP NeuroX simulator")

set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

## global options
option(ENABLE_MPI "Enable MPI-based execution" ON)
option(ENABLE_SOA "Enable SoA Memory Layout" OFF)
option(NEUROX_TIME_STEPPING_VERBOSE "Enable indication of initialization and execution progress" OFF)

## mech / mod options
set(ADDITIONAL_MECHPATH "" CACHE PATH "Search path for optional additional mechanism MOD files")
set(ADDITIONAL_MECHS "" CACHE FILEPATH "File containing list of additional mechanism MOD files")
option(USE_INCLUDED_MECHS "Use mechs in /misc/mechs instead of the ones specified in ADDITIONAL_MECHPATH ADDITIONAL_MECHS" ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
      ${PROJECT_SOURCE_DIR}/CMake/portability
      ${PROJECT_SOURCE_DIR}/CMake/packages
      ${PROJECT_SOURCE_DIR}/CMake/config)

set(Boost_NO_BOOST_CMAKE TRUE)

include(BlueGenePortability)
include(DefineInstallationPaths)
include(ReleaseDebugAutoFlags)
include(SetRpath)

include(FindPkgConfig)
include(CTest)

## Automatically add C++11 flags (eg -libc++ or -lstdc++)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

#Include tclap
FIND_PATH(TCLAP_INCLUDE_DIR tclap/CmdLine.h
        /usr/local/include
        /usr/include
        /opt/local/include
)

IF(TCLAP_INCLUDE_DIR)
        message(STATUS "Found TCLAP: ${TCLAP_INCLUDE_DIR}")
ELSE()
        message(FATAL_ERROR "Could not find TCLAP")
ENDIF()
include_directories(${TCLAP_INCLUDE_DIR})

# Include HPX
#As a requirement for pkg_check_modules we need to provide the path of
#PKG_CONFIG_EXECUTABLE (in this case we pass only the binary as it's in the path)
set(PKG_CONFIG_EXECUTABLE pkg-config)
IF (NOT PKG_CONFIG_FOUND)
  message (FATAL_ERROR "pkg-config not found")
ELSE()
  #use pkg-config to find hpx
  message (STATUS "Found pkg-config ${PKG_CONFIG_VERSION_STRING}: ${PKG_CONFIG_EXECUTABLE}")
  pkg_check_modules(HPX REQUIRED hpx )
  IF (NOT HPX_FOUND)
    message (FATAL_ERROR "pkg-config did not find HPX")
  ELSE()
    #message (STATUS "pkg-config found hpx ${HPX_VERSION}:")
    message (STATUS "  HPX_INCLUDE_DIRS: ${HPX_INCLUDE_DIRS}")
    #message (STATUS "  HPX_LIBRARIES: ${HPX_LIBRARIES}")
    message (STATUS "  HPX_LIBRARY_DIRS: ${HPX_LIBRARY_DIRS}")
    #message (STATUS "  HPX_LDFLAGS: ${HPX_LDFLAGS}")
    #message (STATUS "  ldflags_others: ${HPX_LDFLAGS_OTHERS}")
    #message (STATUS "  HPX_CFLAGS: ${HPX_CFLAGS}")
    #message (STATUS "  cflags_other: ${HPX_CFLAGS_OTHER}") 
    include_directories(${HPX_INCLUDE_DIRS})
    link_directories(${HPX_LIBRARY_DIRS})
  ENDIF()
ENDIF()

#Find all dependencies
if (ENABLE_MPI)
    #### find_package(MPI REQUIRED)  #Not needed, it's part of HPX installation (or hpx_init will fail)
    add_definitions("-DNRNMPI=1")
else()
    message(STATUS "Use of MPI disabled by user-provided option")
    add_definitions("-DNRNMPI=0")
endif()

#Memory layout transformation: 0 for SoA and 1 for AoS
if (ENABLE_SOA)
    add_definitions("-DLAYOUT=0")
else()
    add_definitions("-DLAYOUT=1")
endif()

find_package(mod2c 2.0.0  REQUIRED)
find_package(Boost 1.41.0 QUIET COMPONENTS filesystem system atomic unit_test_framework)

# Some mechanisms use NEURON specific code which is skipped
# using this macro for CoreNeuron build
add_definitions(-DCORENEURON_BUILD)

include(CheckIncludeFiles)

#Turn off inline assembly for PGI
if(${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
    add_definitions(-DSWAP_ENDIAN_DISABLE_ASM)
endif()

#OpenACC pragmas are not guarded, disable all unknown pragm warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${IGNORE_UNKNOWN_PRAGMA_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IGNORE_UNKNOWN_PRAGMA_FLAGS}")

#allow extra function pointers
add_definitions("-DEXPORT_MECHS_FUNCTIONS")

# The mechanism registration functon will be defined depending on whether we
# are using neurodamus list of mechanisms
if(ADDITIONAL_MECHS)
    add_definitions(-DADDITIONAL_MECHS)
endif()

#Exception: add lc++ flag instead of lstdc++ on Apple machines
IF (APPLE)
  message (STATUS "Found MacOS: Added gcc LLVM flag lc++")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lc++")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lc++")
ENDIF()

# Core library
if (NEUROX_TIME_STEPPING_VERBOSE)
    add_definitions("-DNEUROX_TIME_STEPPING_VERBOSE")
endif()

add_subdirectory(neurox)
add_subdirectory(apps)
