## Automatically add C++11 flags (eg -libc++ or -lstdc++)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

#Incldude tclap
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

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB_RECURSE neurox_all_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB_RECURSE neurox_all_templates "${CMAKE_CURRENT_SOURCE_DIR}/*.ipp")
FILE(GLOB_RECURSE neurox_all_c_files "${CMAKE_CURRENT_SOURCE_DIR}/*.c*")

add_library(neurox ${COMPILE_LIBRARY_TYPE} ${neurox_all_headers} ${neurox_all_templates} ${neurox_all_c_files} )

target_link_libraries(neurox coreneuron ${HPX_LIBRARIES})

set_target_properties(neurox PROPERTIES
                            VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                            SOVERSION ${VERSION_MAJOR}
                            CLEAN_DIRECT_OUTPUT 1)

install(TARGETS neurox
             LIBRARY DESTINATION ${LIB_INSTALL_DIR} 
             ARCHIVE DESTINATION ${LIB_INSTALL_DIR} )


#binary (copy from apps/CMakeLists.txt)
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/coreneuron
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}/coreneuron
    ${PROJECT_BINARY_DIR}/neurox
    ${PROJECT_SOURCE_DIR}/neurox
    )

# The mechanism registration function will be defined depending on whether we
# are using neurodamus list of mechanisms
if(ADDITIONAL_MECHS)
    add_definitions(-DADDITIONAL_MECHS)
endif()

FILE(GLOB neurox_exec_src "*.c*")
message(INFO "---------------------- ${neurox_exec_src}")

add_executable(neurox_exec ${neurox_exec_src})

set_target_properties(neurox_exec PROPERTIES
                            OUTPUT_NAME "neurox_exec"
                            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_link_libraries(neurox_exec neurox coreneuron ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${HPX_LIBRARIES} ${reportinglib_LIBRARY})

install(TARGETS neurox_exec
         DESTINATION ${BIN_INSTALL_DIR}/)

