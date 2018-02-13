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

# Include HPX via hpx-config
# As a requirement for pkg_check_modules we need to provide the path of
# PKG_CONFIG_EXECUTABLE (in this case we pass only the binary as it's in the $PATH)
find_package(PkgConfig REQUIRED)
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

# Include CVODES (Cmake/FindCVODES.cmake)
find_package(CVODES 3.1.0 REQUIRED)
include_directories(${CVODES_INCLUDE_DIRS})
link_directories(${CVODES_LIBRARIES})
#message(STATUS "  CVODES_INCLUDE_DIRS: ${CVODES_INCLUDE_DIRS}")
#message(STATUS "  CVODES_LIBRARIES: ${CVODES_LIBRARIES}")

#binary (copy from apps/CMakeLists.txt)
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/coreneuron
    ${PROJECT_SOURCE_DIR}/neurox
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}/coreneuron
    ${PROJECT_BINARY_DIR}/neurox
    )

FILE(GLOB_RECURSE neurox_all_headers   "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB_RECURSE neurox_all_templates "${CMAKE_CURRENT_SOURCE_DIR}/*.ipp")
FILE(GLOB_RECURSE neurox_all_c_files   "${CMAKE_CURRENT_SOURCE_DIR}/*.c*")

#Fix UINT64_C undefined error in bbpviz1.cscs.ch
add_compile_options("-D__STDC_CONSTANT_MACROS")

add_library(neurox ${COMPILE_LIBRARY_TYPE} ${neurox_all_headers} ${neurox_all_templates} ${neurox_all_c_files} )

target_link_libraries(neurox coreneuron ${HPX_LIBRARIES} ${CVODES_LIBRARIES})

set_target_properties(neurox PROPERTIES
                             VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                             SOVERSION ${VERSION_MAJOR}
                             CLEAN_DIRECT_OUTPUT 1)

install(TARGETS neurox
             LIBRARY DESTINATION ${LIB_INSTALL_DIR} 
             ARCHIVE DESTINATION ${LIB_INSTALL_DIR} )

add_executable(neurox_exec "${CMAKE_SOURCE_DIR}/apps/neurox.cc")

set_target_properties(neurox_exec PROPERTIES
                      OUTPUT_NAME "neurox_exec"
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_link_libraries(neurox_exec neurox coreneuron ${HPX_LIBRARIES} )

install(TARGETS neurox_exec  DESTINATION ${BIN_INSTALL_DIR}/)

