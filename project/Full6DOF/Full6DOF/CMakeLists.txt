# CMakeList.txt : CMake project for Full6DOF, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (Full6DOF "Full6DOF.cpp" "Full6DOF.h" "src/util/Vector.h" "src/ode/ODE_standalone.h" "src/util/fast_math.h" "src/simulation/Loader.h" "src/simulation/Simulation.cpp" "src/vehicle/aerodynamics/Aerodynamics.h" "src/vehicle/Component/Vehicle_Component.h"  "src/physics/Ephemeris.cpp" "src/util/functions.h" "src/physics/SolarSystem.cpp" "src/vehicle/Environment.cpp" "src/physics/Body.cpp"  "src/vehicle/propulsion/Propulsion.h" "src/vehicle/propulsion/Thruster.h" "src/physics/Constants.h" "src/physics/Action.cpp" "src/physics/Geometry.cpp" "src/physics/Atmosphere.cpp" "src/physics/Time.cpp" "src/vehicle/propulsion/Tank.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Full6DOF PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
