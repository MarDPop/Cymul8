# CMakeList.txt : CMake project for Full6DOF, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (Full6DOF "Full6DOF.cpp" "Full6DOF.h" "src/util/Vector.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Full6DOF PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
