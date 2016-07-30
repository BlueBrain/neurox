include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB_RECURSE neurox_all_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB_RECURSE neurox_all_templates "${CMAKE_CURRENT_SOURCE_DIR}/*.ipp")
FILE(GLOB_RECURSE neurox_all_c_files "${CMAKE_CURRENT_SOURCE_DIR}/*.c*")

set(neurox_all_c_files ${neurox_all_c_files})
add_library(neurox ${COMPILE_LIBRARY_TYPE} ${neurox_all_headers} ${neurox_all_templates} ${neurox_all_c_files})

target_link_libraries(neurox coreneuron ${HPX_LIBRARIES})

set_target_properties(neurox PROPERTIES
                            VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                            SOVERSION ${VERSION_MAJOR}
                            CLEAN_DIRECT_OUTPUT 1)

install(TARGETS neurox
             LIBRARY DESTINATION ${LIB_INSTALL_DIR} 
             ARCHIVE DESTINATION ${LIB_INSTALL_DIR} )


#binary (copy from apps/CMakeLists.txt)
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/coreneuron ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/coreneuron ${PROJECT_BINARY_DIR}/neurox ${PROJECT_SOURCE_DIR}/neurox)

# The mechanism registration function will be defined depending on whether we
# are using neurodamus list of mechanisms
if(ADDITIONAL_MECHS)
    add_definitions(-DADDITIONAL_MECHS)
endif()

FILE(GLOB neurox_exec_src "*.c*")

add_executable(neurox_exec ${neurox_exec_src})

set_target_properties(neurox_exec PROPERTIES
                            OUTPUT_NAME "neurox_exec"
                            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_link_libraries(neurox_exec coreneuron neurox ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${HPX_LIBRARIES} ${reportinglib_LIBRARY})

install(TARGETS neurox_exec
         DESTINATION ${BIN_INSTALL_DIR}/)

