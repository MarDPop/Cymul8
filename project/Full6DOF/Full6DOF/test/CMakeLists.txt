add_executable (test6DOF "main.cpp")

target_include_directories(test6DOF PUBLIC ../lib/Eigen)

target_link_libraries(test6DOF odeLib physicsLib simulationLib utilsLib vehiclesLib)

target_compile_features(test6DOF PRIVATE cxx_std_20)