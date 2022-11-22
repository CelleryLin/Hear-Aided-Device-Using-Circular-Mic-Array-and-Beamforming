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
include CMakeFiles/rtaudio.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/rtaudio.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rtaudio.dir/flags.make

CMakeFiles/rtaudio.dir/RtAudio.cpp.o: CMakeFiles/rtaudio.dir/flags.make
CMakeFiles/rtaudio.dir/RtAudio.cpp.o: ../RtAudio.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/project/rtaudio-5.2.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rtaudio.dir/RtAudio.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rtaudio.dir/RtAudio.cpp.o -c /home/pi/project/rtaudio-5.2.0/RtAudio.cpp

CMakeFiles/rtaudio.dir/RtAudio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rtaudio.dir/RtAudio.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/project/rtaudio-5.2.0/RtAudio.cpp > CMakeFiles/rtaudio.dir/RtAudio.cpp.i

CMakeFiles/rtaudio.dir/RtAudio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rtaudio.dir/RtAudio.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/project/rtaudio-5.2.0/RtAudio.cpp -o CMakeFiles/rtaudio.dir/RtAudio.cpp.s

CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o: CMakeFiles/rtaudio.dir/flags.make
CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o: ../rtaudio_c.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/project/rtaudio-5.2.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o -c /home/pi/project/rtaudio-5.2.0/rtaudio_c.cpp

CMakeFiles/rtaudio.dir/rtaudio_c.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rtaudio.dir/rtaudio_c.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/project/rtaudio-5.2.0/rtaudio_c.cpp > CMakeFiles/rtaudio.dir/rtaudio_c.cpp.i

CMakeFiles/rtaudio.dir/rtaudio_c.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rtaudio.dir/rtaudio_c.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/project/rtaudio-5.2.0/rtaudio_c.cpp -o CMakeFiles/rtaudio.dir/rtaudio_c.cpp.s

# Object files for target rtaudio
rtaudio_OBJECTS = \
"CMakeFiles/rtaudio.dir/RtAudio.cpp.o" \
"CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o"

# External object files for target rtaudio
rtaudio_EXTERNAL_OBJECTS =

librtaudio.so.6.0.2: CMakeFiles/rtaudio.dir/RtAudio.cpp.o
librtaudio.so.6.0.2: CMakeFiles/rtaudio.dir/rtaudio_c.cpp.o
librtaudio.so.6.0.2: CMakeFiles/rtaudio.dir/build.make
librtaudio.so.6.0.2: /usr/lib/arm-linux-gnueabihf/libasound.so
librtaudio.so.6.0.2: /usr/lib/arm-linux-gnueabihf/libpulse.so
librtaudio.so.6.0.2: /usr/lib/arm-linux-gnueabihf/libpulse-simple.so
librtaudio.so.6.0.2: CMakeFiles/rtaudio.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/project/rtaudio-5.2.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library librtaudio.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rtaudio.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library librtaudio.so.6.0.2 librtaudio.so.6 librtaudio.so

librtaudio.so.6: librtaudio.so.6.0.2
	@$(CMAKE_COMMAND) -E touch_nocreate librtaudio.so.6

librtaudio.so: librtaudio.so.6.0.2
	@$(CMAKE_COMMAND) -E touch_nocreate librtaudio.so

# Rule to build all files generated by this target.
CMakeFiles/rtaudio.dir/build: librtaudio.so

.PHONY : CMakeFiles/rtaudio.dir/build

CMakeFiles/rtaudio.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rtaudio.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rtaudio.dir/clean

CMakeFiles/rtaudio.dir/depend:
	cd /home/pi/project/rtaudio-5.2.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/project/rtaudio-5.2.0 /home/pi/project/rtaudio-5.2.0 /home/pi/project/rtaudio-5.2.0/build /home/pi/project/rtaudio-5.2.0/build /home/pi/project/rtaudio-5.2.0/build/CMakeFiles/rtaudio.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rtaudio.dir/depend

