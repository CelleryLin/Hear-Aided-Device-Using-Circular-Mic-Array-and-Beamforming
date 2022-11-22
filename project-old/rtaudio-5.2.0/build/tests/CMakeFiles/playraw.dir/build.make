# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pi/project/rtaudio-5.2.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pi/project/rtaudio-5.2.0/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/playraw.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/playraw.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/playraw.dir/flags.make

tests/CMakeFiles/playraw.dir/playraw.cpp.o: tests/CMakeFiles/playraw.dir/flags.make
tests/CMakeFiles/playraw.dir/playraw.cpp.o: ../tests/playraw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/project/rtaudio-5.2.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/playraw.dir/playraw.cpp.o"
	cd /home/pi/project/rtaudio-5.2.0/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/playraw.dir/playraw.cpp.o -c /home/pi/project/rtaudio-5.2.0/tests/playraw.cpp

tests/CMakeFiles/playraw.dir/playraw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/playraw.dir/playraw.cpp.i"
	cd /home/pi/project/rtaudio-5.2.0/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/project/rtaudio-5.2.0/tests/playraw.cpp > CMakeFiles/playraw.dir/playraw.cpp.i

tests/CMakeFiles/playraw.dir/playraw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/playraw.dir/playraw.cpp.s"
	cd /home/pi/project/rtaudio-5.2.0/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/project/rtaudio-5.2.0/tests/playraw.cpp -o CMakeFiles/playraw.dir/playraw.cpp.s

# Object files for target playraw
playraw_OBJECTS = \
"CMakeFiles/playraw.dir/playraw.cpp.o"

# External object files for target playraw
playraw_EXTERNAL_OBJECTS =

tests/playraw: tests/CMakeFiles/playraw.dir/playraw.cpp.o
tests/playraw: tests/CMakeFiles/playraw.dir/build.make
tests/playraw: librtaudio.so.6.0.2
tests/playraw: /usr/lib/arm-linux-gnueabihf/libasound.so
tests/playraw: /usr/lib/arm-linux-gnueabihf/libpulse.so
tests/playraw: /usr/lib/arm-linux-gnueabihf/libpulse-simple.so
tests/playraw: tests/CMakeFiles/playraw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/project/rtaudio-5.2.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable playraw"
	cd /home/pi/project/rtaudio-5.2.0/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/playraw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/playraw.dir/build: tests/playraw

.PHONY : tests/CMakeFiles/playraw.dir/build

tests/CMakeFiles/playraw.dir/clean:
	cd /home/pi/project/rtaudio-5.2.0/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/playraw.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/playraw.dir/clean

tests/CMakeFiles/playraw.dir/depend:
	cd /home/pi/project/rtaudio-5.2.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/project/rtaudio-5.2.0 /home/pi/project/rtaudio-5.2.0/tests /home/pi/project/rtaudio-5.2.0/build /home/pi/project/rtaudio-5.2.0/build/tests /home/pi/project/rtaudio-5.2.0/build/tests/CMakeFiles/playraw.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/playraw.dir/depend

