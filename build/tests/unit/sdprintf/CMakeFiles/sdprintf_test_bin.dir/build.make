# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bmagalha/Workspace/neurox

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bmagalha/Workspace/neurox/build

# Include any dependencies generated for this target.
include tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/depend.make

# Include the progress variables for this target.
include tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/progress.make

# Include the compile flags for this target's objects.
include tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/flags.make

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/flags.make
tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o: ../tests/unit/sdprintf/test_sdprintf.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bmagalha/Workspace/neurox/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o"
	cd /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf && /usr/bin/mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o -c /home/bmagalha/Workspace/neurox/tests/unit/sdprintf/test_sdprintf.cpp

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.i"
	cd /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bmagalha/Workspace/neurox/tests/unit/sdprintf/test_sdprintf.cpp > CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.i

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.s"
	cd /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bmagalha/Workspace/neurox/tests/unit/sdprintf/test_sdprintf.cpp -o CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.s

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.requires:
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.requires

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.provides: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.requires
	$(MAKE) -f tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/build.make tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.provides.build
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.provides

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.provides.build: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o

# Object files for target sdprintf_test_bin
sdprintf_test_bin_OBJECTS = \
"CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o"

# External object files for target sdprintf_test_bin
sdprintf_test_bin_EXTERNAL_OBJECTS =

tests/unit/sdprintf/sdprintf_test_bin: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o
tests/unit/sdprintf/sdprintf_test_bin: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/build.make
tests/unit/sdprintf/sdprintf_test_bin: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
tests/unit/sdprintf/sdprintf_test_bin: coreneuron/libcoreneuron.so.0.8.1
tests/unit/sdprintf/sdprintf_test_bin: /usr/lib/libmpi.so
tests/unit/sdprintf/sdprintf_test_bin: /usr/lib/x86_64-linux-gnu/libdl.so
tests/unit/sdprintf/sdprintf_test_bin: /usr/lib/x86_64-linux-gnu/libhwloc.so
tests/unit/sdprintf/sdprintf_test_bin: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable sdprintf_test_bin"
	cd /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sdprintf_test_bin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/build: tests/unit/sdprintf/sdprintf_test_bin
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/build

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/requires: tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/test_sdprintf.cpp.o.requires
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/requires

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/clean:
	cd /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf && $(CMAKE_COMMAND) -P CMakeFiles/sdprintf_test_bin.dir/cmake_clean.cmake
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/clean

tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/depend:
	cd /home/bmagalha/Workspace/neurox/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bmagalha/Workspace/neurox /home/bmagalha/Workspace/neurox/tests/unit/sdprintf /home/bmagalha/Workspace/neurox/build /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf /home/bmagalha/Workspace/neurox/build/tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/unit/sdprintf/CMakeFiles/sdprintf_test_bin.dir/depend

