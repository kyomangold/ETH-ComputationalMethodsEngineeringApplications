# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

add_library(base_template
    "boundaries.cpp"
    "coefficient_functions.cpp"
    "initial_data.cpp"
    "run_boundaries.cpp"
    "run_stability.cpp")
target_link_libraries(base_template Eigen::headers)
target_include_directories(base_template PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
