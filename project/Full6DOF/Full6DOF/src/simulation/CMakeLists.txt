add_library(simulationLib STATIC
    Simulation.cpp
)

set_target_properties(simulationLib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(simulationLib PUBLIC .)

target_compile_features(simulationLib PRIVATE cxx_std_20)

install(TARGETS simulationLib
    LIBRARY DESTINATION lib)
