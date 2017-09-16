# Tries to find Sundials CVODE.
#
# This module will define the following variables:
#  CVODES_INCLUDE_DIRS - Location of the CVODE includes
#  CVODES_FOUND - true if CVODE was found on the system
#  CVODES_LIBRARIES - Required libraries for all requested components
#
# This module accepts the following environment or CMake vars
#  CVODES_ROOT - Install location to search for

include(FindPackageHandleStandardArgs)

if(NOT "$ENV{CVODES_ROOT}" STREQUAL "" OR NOT "${CVODES_ROOT}" STREQUAL "")
    list(APPEND CMAKE_INCLUDE_PATH "$ENV{CVODES_ROOT}" "${CVODES_ROOT}")
    list(APPEND CMAKE_LIBRARY_PATH "$ENV{CVODES_ROOT}" "${CVODES_ROOT}")
endif()

find_path(CVODES_INCLUDE_DIRS cvodes/cvodes.h
    ENV CVODES_ROOT
    PATH_SUFFIXES include include/cvodes
)

find_library(CVODES_LIBRARY
    NAMES sundials_cvodes
    ENV CVODES_ROOT
    PATH_SUFFIXES lib Lib
)

find_library(CVODES_NVECSERIAL
    NAMES sundials_nvecserial
    ENV CVODES_ROOT
    PATH_SUFFIXES lib Lib
)

find_library(CVODES_KLU klu)

find_package_handle_standard_args(CVODES DEFAULT_MSG
    CVODES_LIBRARY
    CVODES_NVECSERIAL
    CVODES_INCLUDE_DIRS
)

if(CVODES_FOUND)
    set(CVODES_LIBRARIES ${CVODES_LIBRARY} ${CVODES_NVECSERIAL} CACHE INTERNAL "")
    if(CVODES_KLU)
        set(CVODES_LIBRARIES ${CVODES_LIBRARIES} ${CVODES_KLU} CACHE INTERNAL "")
    endif()
endif()

mark_as_advanced(CVODES_INCLUDE_DIRS CVODES_LIBRARY CVODES_NVECSERIAL)

