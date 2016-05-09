include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

#message (INFO "-------------- ${CMAKE_SOURCE_DIR}")
#message (INFO "-------------- ${CMAKE_CURRENT_BINARY_DIR}")
#message (INFO "-------------- ${CMAKE_CURRENT_SOURCE_DIR}")

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
