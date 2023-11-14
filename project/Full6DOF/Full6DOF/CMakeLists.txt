# CMakeList.txt : CMake project for Full6DOF, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (Full6DOF "Full6DOF.cpp" "Full6DOF.h")

add_subdirectory (src/ode)
add_subdirectory (src/physics)
add_subdirectory (src/simulation)
add_subdirectory (src/util)
add_subdirectory (src/vehicle)
add_subdirectory (test)

target_include_directories(Full6DOF PUBLIC lib/Eigen lib/fpm lib/tinyxml)

target_link_libraries(Full6DOF odeLib physicsLib simulationLib utilsLib vehiclesLib)

set_property(TARGET Full6DOF PROPERTY CXX_STANDARD 20)

# TODO: Add tests and install targets if needed.
