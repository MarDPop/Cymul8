add_library(utilsLib STATIC
    functions.cpp
)

set_target_properties(utilsLib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(utilsLib PUBLIC .)

target_compile_features(utilsLib PRIVATE cxx_std_20)

install(TARGETS utilsLib
    LIBRARY DESTINATION lib)
