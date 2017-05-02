include_directories(${PROJECT_SOURCE_DIR})    #include neurox/neurox.h
link_directories(${CMAKE_INSTALL_PREFIX}/lib) #-lcoreneuron

add_executable(neurox_exec "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
set_target_properties(neurox_exec PROPERTIES
                            OUTPUT_NAME "neurox"
                            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(neurox_exec neurox )
install(TARGETS neurox_exec DESTINATION ${BIN_INSTALL_DIR}/)
