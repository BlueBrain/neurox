include_directories(${PROJECT_SOURCE_DIR}             #include neurox/neurox.h
                    ${PROJECT_SOURCE_DIR}/neurox      #include neurox.h
                    ${CMAKE_INSTALL_PREFIX}/include   #coreneuron/coreneuron.h
                    )
link_directories(${CMAKE_INSTALL_PREFIX}/lib)         #-lcoreneuron

FILE(GLOB neurox_exec_src "main.c*")
add_executable(neurox_exec ${neurox_exec_src})

set_target_properties(neurox_exec PROPERTIES
                            OUTPUT_NAME "neurox"
                            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_link_libraries(neurox_exec neurox )

install(TARGETS neurox_exec 
        DESTINATION ${BIN_INSTALL_DIR}/)
