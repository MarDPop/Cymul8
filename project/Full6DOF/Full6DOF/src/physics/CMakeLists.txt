cmake_minimum_required(VERSION 3.12)

add_library(physicsLib STATIC
    Action.cpp
    Atmosphere.cpp
    Body.cpp
    Ephemeris.cpp
    Geometry.cpp
    SolarSystem.cpp
    "EpochTime.cpp"
)

set_target_properties(physicsLib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(physicsLib PUBLIC .)

target_compile_features(physicsLib PRIVATE cxx_std_20)

install(TARGETS physicsLib
    LIBRARY DESTINATION lib)

