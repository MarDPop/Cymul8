add_library(vehiclesLib STATIC
    Environment.cpp
    Vehicle.cpp
)

set_target_properties(vehiclesLib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(vehiclesLib PUBLIC .
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/aerodynamics>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/component>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/propulsion>)


target_compile_features(vehiclesLib PRIVATE cxx_std_20)

install(TARGETS vehiclesLib
    LIBRARY DESTINATION lib)
