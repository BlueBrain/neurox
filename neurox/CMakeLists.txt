FILE(GLOB_RECURSE neurox_all_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB_RECURSE neurox_all_templates "${CMAKE_CURRENT_SOURCE_DIR}/*.ipp")
FILE(GLOB_RECURSE neurox_all_c_files "${CMAKE_CURRENT_SOURCE_DIR}/*.c*")

include_directories( ${CMAKE_CURRENT_BINARY_DIR}
                     ${CMAKE_SOURCE_DIR}              #include neurox/neurox.h
                     ${CMAKE_CURRENT_SOURCE_DIR}      #include neurox.h TODO remove this
                     ${CMAKE_INSTALL_PREFIX}/include 
                     ${CMAKE_INSTALL_PREFIX}/include/coreneuron/utils/randoms 
                   )
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

#################### Generate mechanisms mod files

# Compile and include MOD files
find_program(MOD2C mod2c_core)
if(MOD2C STREQUAL "MOD2C-NOTFOUND")
message(FATAL_ERROR "Unable to locate mod2c_core executable")
endif()

# mod2c_core wants MODLUNIT set to point to its units file,
# but CMake does not support setting an environment for custom_commands.

get_filename_component(mod2c_install_prefix ${MOD2C} DIRECTORY)
find_file(MOD2C_UNITS nrnunits.lib HINT "${mod2c_install_prefix}" PATH_SUFFIXES "share")
if(NOT MOD2C_UNITS STREQUAL "MOD2C_UNITS-NOTFOUND")
set(MOD2C env "MODLUNIT=${MOD2C_UNITS}" ${MOD2C})
endif()

# Macro sets up build rule for .c files from .mod files.
# Parameters:
#    name     An arbitrary name to keep track of output .c files
#    input    The path to the mod file
#
# Because mod2c_core wants to write the .c file in the same directory
# as the mod file, we copy the mod file to the binary directory first
#
# The macro appends the names of the output files to MOD2C_${name}_OUTPUTS
# and the names of the mod files (without directories) to MOD2C_${name}_MODS

macro(mod2c_target name input)
    get_filename_component(mod2c_source_ ${input} ABSOLUTE)
    get_filename_component(mod2c_modname_ ${input} NAME)
    string(REGEX REPLACE "\\.mod$" ".c" mod2c_cname_ "${mod2c_modname_}")

    set(mod2c_output_ "${CMAKE_CURRENT_BINARY_DIR}/${mod2c_cname_}")
    list(APPEND MOD2C_${name}_OUTPUTS "${mod2c_output_}")
    list(APPEND MOD2C_${name}_MODS "${mod2c_modname_}")

    add_custom_command(OUTPUT "${mod2c_output_}"
      DEPENDS ${input} "${mod2c_install_prefix}/mod2c_core"
      COMMAND ${CMAKE_COMMAND} -E copy "${mod2c_source_}" "${CMAKE_CURRENT_BINARY_DIR}"
      COMMAND ${MOD2C} "${mod2c_modname_}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()

# Macro mod2c_from_file attempts to add modules specified in a file.
# Parameters:
#     name        Key for mod2c_target
#     modlist     Path to the file containing the list of mod files.
#     searchpath  List of directories to search for mod files.
#
# This has to be a macro in order to properly accumulate the MOD2C_xxx_OUTPUTS
# variables. Go CMake! You rock!

macro(mod2c_from_file name modlist searchpath)
    unset(tmp_path_ CACHE)
    if(NOT EXISTS "${modlist}")
        string(REGEX REPLACE ":" ";" spath "${searchpath}")
        foreach(path_ ${spath})
            file(GLOB_RECURSE mod_files "${path_}/*.mod")
            list(APPEND mods_ "${mod_files};")
        endforeach()
    else()
        file(READ "${modlist}" mods_)
    endif()
    string(REGEX REPLACE "[ \t\n]+" ";" mods_ "${mods_}")
    list(REMOVE_ITEM mods_ "")

    foreach(mod_ ${mods_})
        unset(modpath_)
        if(NOT IS_ABSOLUTE "${mod_}")
            find_path(tmp_path_ "${mod_}" PATH ${searchpath} NO_DEFAULT_PATH)
            if(tmp_path_)
                set(modpath_ "${tmp_path_}/${mod_}")
            endif()
            unset(tmp_path_ CACHE)
        else()
            set(modpath_ "${mod_}")
        endif()
        if((NOT EXISTS "${modpath_}") OR (IS_DIRECTORY "${modpath_}"))
            unset(modpath_)
        endif()

        if(modpath_)
            mod2c_target(${name} "${modpath_}")
            list(APPEND MOD_PATHS ${modpath_})
            message("-- extra mechanism: ${modpath_}")
        else()
            message(WARNING "Unable to locate MOD file: ${mod_}")
        endif()
    endforeach()
endmacro()

# Standard MOD files are always built

set(STDMODLIBPATH "${CMAKE_INSTALL_PREFIX}/share/mod")
set(MECH_STDMODS expsyn.mod  hh.mod  netstim.mod  passive.mod  pattern.mod stim.mod)

foreach(MOD ${MECH_STDMODS})
    message("-- basic mechanism: ${STDMODLIBPATH}/${MOD}")
    list(APPEND MOD_PATHS "${STDMODLIBPATH}/${MOD}")
    mod2c_target(STDMECH "${STDMODLIBPATH}/${MOD}")
endforeach()

set(MECH_FILE_LISTS)

# Additional optional mechanisms specified in ADDITIONAL_MECHS
# and ADDITIONAL_MECHPATH (set as a PATH cache variable, but
# can in fact be a list of paths.)
# If ADDITIONAL_MECHS is not provided, the entire ADDITIONAL_MECHPATH
# will be checked for the mod files.


if(EXISTS ${ADDITIONAL_MECHS})
  if(IS_DIRECTORY "${modlist}")
      message(WARNING "Unable to read MOD list file: ${modlist}")
  else()
      list(APPEND MECH_FILE_LISTS "${ADDITIONAL_MECHS}")
  endif()
  mod2c_from_file(OPTMECH "${ADDITIONAL_MECHS}" "${ADDITIONAL_MECHPATH}")
else()
  set(ADDITIONAL_MECHPATH   "${CMAKE_SOURCE_DIR}/misc/mechs/coreneuron_modlist.txt")
  set(ADDITIONAL_MECHS      "${CMAKE_SOURCE_DIR}/misc/mechs/mod")
  FILE(GLOB_RECURSE MOD2C_STDMECH_OUTPUTS "${CMAKE_SOURCE_DIR}/misc/mechs/c/*.c*")
  message("-- Using pre-existing mechanisms from ${CMAKE_SOURCE_DIR}/misc/mechs")
endif()

set(MOD_FUNC_PTRS_C    "${CMAKE_CURRENT_BINARY_DIR}/mod_func_ptrs.c")
set(MOD_FUNC_PTRS_C_PL "${CMAKE_SOURCE_DIR}/misc/mod_func_ptrs.c.pl")
message("-- mechanisms script: ${MOD_FUNC_PTRS_C_PL}")
message("-- mechanisms output: ${MOD_FUNC_PTRS_C}")

# ... pass as argument to the perl script the list of mods full paths;
add_custom_command(OUTPUT "${MOD_FUNC_PTRS_C}"
    DEPENDS ${MECH_FILE_LISTS} "${MOD_FUNC_PTRS_C_PL}"
    COMMAND perl "${MOD_FUNC_PTRS_C_PL}" ${MOD_PATHS} > "${MOD_FUNC_PTRS_C}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add generated sources from MOD files
set(GENERATED_MECH_C_FILES ${MOD_FUNC_PTRS_C} ${MOD2C_STDMECH_OUTPUTS} ${MOD2C_OPTMECH_OUTPUTS})
set(neurox_all_c_files ${neurox_all_c_files} ${GENERATED_MECH_C_FILES})

##################### library

add_library(neurox ${COMPILE_LIBRARY_TYPE} ${neurox_all_headers} ${neurox_all_templates} ${neurox_all_c_files} )

target_link_libraries(neurox coreneuron ${HPX_LIBRARIES})

set_target_properties(neurox PROPERTIES
                            VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
                            SOVERSION ${VERSION_MAJOR}
                            CLEAN_DIRECT_OUTPUT 1)

install(TARGETS neurox
        LIBRARY DESTINATION ${LIB_INSTALL_DIR} 
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR} )
