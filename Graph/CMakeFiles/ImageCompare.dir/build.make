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
CMAKE_COMMAND = /usr/sci/projects/neuro/software/cmake-2.8.7-bin/bin/cmake

# The command to remove a file.
RM = /usr/sci/projects/neuro/software/cmake-2.8.7-bin/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/sci/projects/neuro/software/cmake-2.8.7-bin/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph

# Include any dependencies generated for this target.
include CMakeFiles/ImageCompare.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ImageCompare.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ImageCompare.dir/flags.make

CMakeFiles/ImageCompare.dir/ImageCompare.o: CMakeFiles/ImageCompare.dir/flags.make
CMakeFiles/ImageCompare.dir/ImageCompare.o: ImageCompare.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/ImageCompare.dir/ImageCompare.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/ImageCompare.dir/ImageCompare.o -c /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/ImageCompare.cxx

CMakeFiles/ImageCompare.dir/ImageCompare.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImageCompare.dir/ImageCompare.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/ImageCompare.cxx > CMakeFiles/ImageCompare.dir/ImageCompare.i

CMakeFiles/ImageCompare.dir/ImageCompare.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImageCompare.dir/ImageCompare.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/ImageCompare.cxx -o CMakeFiles/ImageCompare.dir/ImageCompare.s

CMakeFiles/ImageCompare.dir/ImageCompare.o.requires:
.PHONY : CMakeFiles/ImageCompare.dir/ImageCompare.o.requires

CMakeFiles/ImageCompare.dir/ImageCompare.o.provides: CMakeFiles/ImageCompare.dir/ImageCompare.o.requires
	$(MAKE) -f CMakeFiles/ImageCompare.dir/build.make CMakeFiles/ImageCompare.dir/ImageCompare.o.provides.build
.PHONY : CMakeFiles/ImageCompare.dir/ImageCompare.o.provides

CMakeFiles/ImageCompare.dir/ImageCompare.o.provides.build: CMakeFiles/ImageCompare.dir/ImageCompare.o

# Object files for target ImageCompare
ImageCompare_OBJECTS = \
"CMakeFiles/ImageCompare.dir/ImageCompare.o"

# External object files for target ImageCompare
ImageCompare_EXTERNAL_OBJECTS =

ImageCompare: CMakeFiles/ImageCompare.dir/ImageCompare.o
ImageCompare: /usr/local/lib/libITKCommon-4.4.so.1
ImageCompare: /usr/local/lib/libitksys-4.4.so.1
ImageCompare: /usr/local/lib/libITKVNLInstantiation-4.4.so.1
ImageCompare: /usr/local/lib/libitkvnl_algo-4.4.so.1
ImageCompare: /usr/local/lib/libitkv3p_lsqr-4.4.so.1
ImageCompare: /usr/local/lib/libitkvnl-4.4.so.1
ImageCompare: /usr/local/lib/libitkvcl-4.4.so.1
ImageCompare: /usr/local/lib/libitkv3p_netlib-4.4.so.1
ImageCompare: /usr/local/lib/libitkdouble-conversion-4.4.so.1
ImageCompare: CMakeFiles/ImageCompare.dir/build.make
ImageCompare: CMakeFiles/ImageCompare.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ImageCompare"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ImageCompare.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ImageCompare.dir/build: ImageCompare
.PHONY : CMakeFiles/ImageCompare.dir/build

CMakeFiles/ImageCompare.dir/requires: CMakeFiles/ImageCompare.dir/ImageCompare.o.requires
.PHONY : CMakeFiles/ImageCompare.dir/requires

CMakeFiles/ImageCompare.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ImageCompare.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ImageCompare.dir/clean

CMakeFiles/ImageCompare.dir/depend:
	cd /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/CMakeFiles/ImageCompare.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ImageCompare.dir/depend

