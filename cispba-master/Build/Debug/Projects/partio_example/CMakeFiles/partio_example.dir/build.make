# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug

# Include any dependencies generated for this target.
include Projects/partio_example/CMakeFiles/partio_example.dir/depend.make

# Include the progress variables for this target.
include Projects/partio_example/CMakeFiles/partio_example.dir/progress.make

# Include the compile flags for this target's objects.
include Projects/partio_example/CMakeFiles/partio_example.dir/flags.make

Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.o: Projects/partio_example/CMakeFiles/partio_example.dir/flags.make
Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.o: ../../Projects/partio_example/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.o"
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/partio_example.dir/main.cpp.o -c /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Projects/partio_example/main.cpp

Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/partio_example.dir/main.cpp.i"
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Projects/partio_example/main.cpp > CMakeFiles/partio_example.dir/main.cpp.i

Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/partio_example.dir/main.cpp.s"
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Projects/partio_example/main.cpp -o CMakeFiles/partio_example.dir/main.cpp.s

# Object files for target partio_example
partio_example_OBJECTS = \
"CMakeFiles/partio_example.dir/main.cpp.o"

# External object files for target partio_example
partio_example_EXTERNAL_OBJECTS =

../../Projects/partio_example/partio_example_debug: Projects/partio_example/CMakeFiles/partio_example.dir/main.cpp.o
../../Projects/partio_example/partio_example_debug: Projects/partio_example/CMakeFiles/partio_example.dir/build.make
../../Projects/partio_example/partio_example_debug: partio-build/lib/libpartio.a
../../Projects/partio_example/partio_example_debug: Projects/partio_example/CMakeFiles/partio_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../Projects/partio_example/partio_example_debug"
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/partio_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Projects/partio_example/CMakeFiles/partio_example.dir/build: ../../Projects/partio_example/partio_example_debug

.PHONY : Projects/partio_example/CMakeFiles/partio_example.dir/build

Projects/partio_example/CMakeFiles/partio_example.dir/clean:
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example && $(CMAKE_COMMAND) -P CMakeFiles/partio_example.dir/cmake_clean.cmake
.PHONY : Projects/partio_example/CMakeFiles/partio_example.dir/clean

Projects/partio_example/CMakeFiles/partio_example.dir/depend:
	cd /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Projects/partio_example /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example /Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Debug/Projects/partio_example/CMakeFiles/partio_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Projects/partio_example/CMakeFiles/partio_example.dir/depend

