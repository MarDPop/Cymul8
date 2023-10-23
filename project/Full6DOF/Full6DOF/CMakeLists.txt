# CMakeList.txt : CMake project for Full6DOF, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (Full6DOF "Full6DOF.cpp" "Full6DOF.h" "src/util/Vector.h" "src/ode/ODE_standalone.h" "src/util/fast_math.h" "src/simulation/Loader.h" "src/simulation/Simulation.cpp" "src/vehicle/6DOF/Aerodynamics.h" "src/vehicle/6DOF/Vehicle_6DOF.h" "src/vehicle/3DOF/Vehicle_3DOF.h" "src/vehicle/Component/Vehicle_Component.h" )

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Full6DOF PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
