# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/xl/quadprog_test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xl/quadprog_test/build

# Include any dependencies generated for this target.
include CMakeFiles/demo.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/demo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/demo.dir/flags.make

CMakeFiles/demo.dir/src/rtGetInf.c.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/src/rtGetInf.c.o: ../src/rtGetInf.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/demo.dir/src/rtGetInf.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/demo.dir/src/rtGetInf.c.o   -c /home/xl/quadprog_test/src/rtGetInf.c

CMakeFiles/demo.dir/src/rtGetInf.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/demo.dir/src/rtGetInf.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xl/quadprog_test/src/rtGetInf.c > CMakeFiles/demo.dir/src/rtGetInf.c.i

CMakeFiles/demo.dir/src/rtGetInf.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/demo.dir/src/rtGetInf.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xl/quadprog_test/src/rtGetInf.c -o CMakeFiles/demo.dir/src/rtGetInf.c.s

CMakeFiles/demo.dir/src/rtGetInf.c.o.requires:

.PHONY : CMakeFiles/demo.dir/src/rtGetInf.c.o.requires

CMakeFiles/demo.dir/src/rtGetInf.c.o.provides: CMakeFiles/demo.dir/src/rtGetInf.c.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/src/rtGetInf.c.o.provides.build
.PHONY : CMakeFiles/demo.dir/src/rtGetInf.c.o.provides

CMakeFiles/demo.dir/src/rtGetInf.c.o.provides.build: CMakeFiles/demo.dir/src/rtGetInf.c.o


CMakeFiles/demo.dir/src/rtGetNaN.c.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/src/rtGetNaN.c.o: ../src/rtGetNaN.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/demo.dir/src/rtGetNaN.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/demo.dir/src/rtGetNaN.c.o   -c /home/xl/quadprog_test/src/rtGetNaN.c

CMakeFiles/demo.dir/src/rtGetNaN.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/demo.dir/src/rtGetNaN.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xl/quadprog_test/src/rtGetNaN.c > CMakeFiles/demo.dir/src/rtGetNaN.c.i

CMakeFiles/demo.dir/src/rtGetNaN.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/demo.dir/src/rtGetNaN.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xl/quadprog_test/src/rtGetNaN.c -o CMakeFiles/demo.dir/src/rtGetNaN.c.s

CMakeFiles/demo.dir/src/rtGetNaN.c.o.requires:

.PHONY : CMakeFiles/demo.dir/src/rtGetNaN.c.o.requires

CMakeFiles/demo.dir/src/rtGetNaN.c.o.provides: CMakeFiles/demo.dir/src/rtGetNaN.c.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/src/rtGetNaN.c.o.provides.build
.PHONY : CMakeFiles/demo.dir/src/rtGetNaN.c.o.provides

CMakeFiles/demo.dir/src/rtGetNaN.c.o.provides.build: CMakeFiles/demo.dir/src/rtGetNaN.c.o


CMakeFiles/demo.dir/src/rt_nonfinite.c.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/src/rt_nonfinite.c.o: ../src/rt_nonfinite.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/demo.dir/src/rt_nonfinite.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/demo.dir/src/rt_nonfinite.c.o   -c /home/xl/quadprog_test/src/rt_nonfinite.c

CMakeFiles/demo.dir/src/rt_nonfinite.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/demo.dir/src/rt_nonfinite.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xl/quadprog_test/src/rt_nonfinite.c > CMakeFiles/demo.dir/src/rt_nonfinite.c.i

CMakeFiles/demo.dir/src/rt_nonfinite.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/demo.dir/src/rt_nonfinite.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xl/quadprog_test/src/rt_nonfinite.c -o CMakeFiles/demo.dir/src/rt_nonfinite.c.s

CMakeFiles/demo.dir/src/rt_nonfinite.c.o.requires:

.PHONY : CMakeFiles/demo.dir/src/rt_nonfinite.c.o.requires

CMakeFiles/demo.dir/src/rt_nonfinite.c.o.provides: CMakeFiles/demo.dir/src/rt_nonfinite.c.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/src/rt_nonfinite.c.o.provides.build
.PHONY : CMakeFiles/demo.dir/src/rt_nonfinite.c.o.provides

CMakeFiles/demo.dir/src/rt_nonfinite.c.o.provides.build: CMakeFiles/demo.dir/src/rt_nonfinite.c.o


CMakeFiles/demo.dir/src/test_quadp.c.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/src/test_quadp.c.o: ../src/test_quadp.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/demo.dir/src/test_quadp.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/demo.dir/src/test_quadp.c.o   -c /home/xl/quadprog_test/src/test_quadp.c

CMakeFiles/demo.dir/src/test_quadp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/demo.dir/src/test_quadp.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xl/quadprog_test/src/test_quadp.c > CMakeFiles/demo.dir/src/test_quadp.c.i

CMakeFiles/demo.dir/src/test_quadp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/demo.dir/src/test_quadp.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xl/quadprog_test/src/test_quadp.c -o CMakeFiles/demo.dir/src/test_quadp.c.s

CMakeFiles/demo.dir/src/test_quadp.c.o.requires:

.PHONY : CMakeFiles/demo.dir/src/test_quadp.c.o.requires

CMakeFiles/demo.dir/src/test_quadp.c.o.provides: CMakeFiles/demo.dir/src/test_quadp.c.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/src/test_quadp.c.o.provides.build
.PHONY : CMakeFiles/demo.dir/src/test_quadp.c.o.provides

CMakeFiles/demo.dir/src/test_quadp.c.o.provides.build: CMakeFiles/demo.dir/src/test_quadp.c.o


CMakeFiles/demo.dir/examples/main.c.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/examples/main.c.o: ../examples/main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/demo.dir/examples/main.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/demo.dir/examples/main.c.o   -c /home/xl/quadprog_test/examples/main.c

CMakeFiles/demo.dir/examples/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/demo.dir/examples/main.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xl/quadprog_test/examples/main.c > CMakeFiles/demo.dir/examples/main.c.i

CMakeFiles/demo.dir/examples/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/demo.dir/examples/main.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xl/quadprog_test/examples/main.c -o CMakeFiles/demo.dir/examples/main.c.s

CMakeFiles/demo.dir/examples/main.c.o.requires:

.PHONY : CMakeFiles/demo.dir/examples/main.c.o.requires

CMakeFiles/demo.dir/examples/main.c.o.provides: CMakeFiles/demo.dir/examples/main.c.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/examples/main.c.o.provides.build
.PHONY : CMakeFiles/demo.dir/examples/main.c.o.provides

CMakeFiles/demo.dir/examples/main.c.o.provides.build: CMakeFiles/demo.dir/examples/main.c.o


# Object files for target demo
demo_OBJECTS = \
"CMakeFiles/demo.dir/src/rtGetInf.c.o" \
"CMakeFiles/demo.dir/src/rtGetNaN.c.o" \
"CMakeFiles/demo.dir/src/rt_nonfinite.c.o" \
"CMakeFiles/demo.dir/src/test_quadp.c.o" \
"CMakeFiles/demo.dir/examples/main.c.o"

# External object files for target demo
demo_EXTERNAL_OBJECTS =

demo: CMakeFiles/demo.dir/src/rtGetInf.c.o
demo: CMakeFiles/demo.dir/src/rtGetNaN.c.o
demo: CMakeFiles/demo.dir/src/rt_nonfinite.c.o
demo: CMakeFiles/demo.dir/src/test_quadp.c.o
demo: CMakeFiles/demo.dir/examples/main.c.o
demo: CMakeFiles/demo.dir/build.make
demo: CMakeFiles/demo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xl/quadprog_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C executable demo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/demo.dir/build: demo

.PHONY : CMakeFiles/demo.dir/build

CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/src/rtGetInf.c.o.requires
CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/src/rtGetNaN.c.o.requires
CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/src/rt_nonfinite.c.o.requires
CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/src/test_quadp.c.o.requires
CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/examples/main.c.o.requires

.PHONY : CMakeFiles/demo.dir/requires

CMakeFiles/demo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/demo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/demo.dir/clean

CMakeFiles/demo.dir/depend:
	cd /home/xl/quadprog_test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xl/quadprog_test /home/xl/quadprog_test /home/xl/quadprog_test/build /home/xl/quadprog_test/build /home/xl/quadprog_test/build/CMakeFiles/demo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/demo.dir/depend

