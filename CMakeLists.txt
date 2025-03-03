cmake_minimum_required(VERSION 3.27)

# Set the target architecture.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the project name
project(gitBooleans)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)
add_subdirectory(arrangements/external/oneTBB)

# Add the executable
add_executable(${PROJECT_NAME} main.cpp
        code/intersect_custom.h code/intersect_custom.cpp
        code/intersect_point_rationals.h code/intersect_point_rationals.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
        ./
        code/
        mesh_booleans/
        arrangements/code/
)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_CGAL_GMP_MPFR ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/cinolib)
find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)
target_link_libraries(${PROJECT_NAME} tbb)
target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Add the executable for tests
add_executable(${PROJECT_NAME}_tests main_tests.cpp
        code/intersect_custom.h code/intersect_custom.cpp
        code/intersect_point_rationals.h code/intersect_point_rationals.cpp)

target_include_directories(${PROJECT_NAME}_tests PUBLIC
        ./
        code/
        mesh_booleans/
        arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_tests cinolib)
target_link_libraries(${PROJECT_NAME}_tests tbb)
target_compile_definitions(${PROJECT_NAME}_tests PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_tests PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_tests PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Add the executable for parse_diff
add_executable(${PROJECT_NAME}_parse_diff main_parse_diff.cpp
        code/intersect_custom.h code/intersect_custom.cpp
        code/intersect_point_rationals.h code/intersect_point_rationals.cpp)

target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC
        ./
        code/
        mesh_booleans/
        arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_parse_diff cinolib)
target_link_libraries(${PROJECT_NAME}_parse_diff tbb)
target_compile_definitions(${PROJECT_NAME}_parse_diff PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
    target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
    if(ENABLE_AVX2)
        target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
    elseif(ENABLE_SSE2)
        target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
    endif()
    target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
    target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
else()
    target_compile_options(${PROJECT_NAME} PUBLIC -O2)
    target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
    target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
    if(ENABLE_AVX2)
        target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
    elseif(ENABLE_SSE2)
        target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
    endif()
endif()