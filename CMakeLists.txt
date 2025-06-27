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

# Path to external mesh boolean library
set(EXT_BOOLS_PATH "/Users/michele/Documents/GitHub/InteractiveAndRobustMeshBooleans")


# Add the executable
add_executable(${PROJECT_NAME} main.cpp
        ${EXT_BOOLS_PATH}/rationals_code/intersect_custom.h  ${EXT_BOOLS_PATH}/rationals_code/intersect_custom.cpp
        ${EXT_BOOLS_PATH}/rationals_code/intersect_point_rationals.h  ${EXT_BOOLS_PATH}/rationals_code/intersect_point_rationals.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
        ./
        #code/
        #mesh_booleans/
        #arrangements/code/
        ${EXT_BOOLS_PATH}/code/
        ${EXT_BOOLS_PATH}/rationals_code/
        ${EXT_BOOLS_PATH}/arrangements/code/
        ${EXT_BOOLS_PATH}/arrangements/external/NFG/include/
        ${EXT_BOOLS_PATH}/arrangements/external/Indirect_Predicates/include/
)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_CGAL_GMP_MPFR ON)
set(CINOLIB_USES_INDIRECT_PREDICATES OFF)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/cinolib)
find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)
target_link_libraries(${PROJECT_NAME} tbb)
target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Add the executable for tests
add_executable(${PROJECT_NAME}_tests main_tests.cpp
        #code/intersect_custom.h code/intersect_custom.cpp
        #code/intersect_point_rationals.h code/intersect_point_rationals.cpp
)

target_include_directories(${PROJECT_NAME}_tests PUBLIC
        ./
        #code/
        #mesh_booleans/
        #arrangements/code/
        ${EXT_BOOLS_PATH}/code/
        ${EXT_BOOLS_PATH}/rationals_code/
        ${EXT_BOOLS_PATH}/arrangements/code/
        ${EXT_BOOLS_PATH}/arrangements/external/
)

target_link_libraries(${PROJECT_NAME}_tests cinolib)
target_link_libraries(${PROJECT_NAME}_tests tbb)
target_compile_definitions(${PROJECT_NAME}_tests PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_tests PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_tests PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Add the executable for parse_diff
add_executable(${PROJECT_NAME}_parse_diff main_parse_diff.cpp
        ${EXT_BOOLS_PATH}/rationals_code/intersect_custom.h  ${EXT_BOOLS_PATH}/rationals_code/intersect_custom.cpp
        ${EXT_BOOLS_PATH}/rationals_code/intersect_point_rationals.h  ${EXT_BOOLS_PATH}/rationals_code/intersect_point_rationals.cpp
)

target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC
        ./
        #code/
        #mesh_booleans/
        #arrangements/code/
        ${EXT_BOOLS_PATH}/code/
        ${EXT_BOOLS_PATH}/rationals_code/
        ${EXT_BOOLS_PATH}/arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_parse_diff cinolib)
target_link_libraries(${PROJECT_NAME}_parse_diff tbb)
target_compile_definitions(${PROJECT_NAME}_parse_diff PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_parse_diff PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

# Compiler-specific options
#if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
 #   target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
  #  target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
   # if(ENABLE_AVX2)
    #    target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
    # elseif(ENABLE_SSE2)
#   target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
#    endif()
#    target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
#    target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
#else()
#    target_compile_options(${PROJECT_NAME} PUBLIC -O2)
#    target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)

    # Solo su architetture x86_64, aggiungi -mavx2 o -msse2
#   if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
#       if(ENABLE_AVX2)
#           target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
#       elseif(ENABLE_SSE2)
#           target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
#       endif()
#   endif()

    # Opzionale: linker option corretta per macOS
    # Opzionale: linker option corretta per macOS
#    target_link_options(${PROJECT_NAME} PUBLIC "-Wl,-stack_size,0x800000") # 8MB = 0x800000
    # 8MB allineato

#endif()


if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # --- Windows / MSVC ---
    target_compile_options(${PROJECT_NAME} PUBLIC
            "/fp:strict"
            "/Oi"
            "/D _CRT_SECURE_NO_WARNINGS"
    )

    if(ENABLE_AVX2)
        target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
    elseif(ENABLE_SSE2)
        target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
    endif()

    # Stack size for Windows (8MB)
    target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")

else()
    # --- macOS & Linux (Clang or GCC) ---
    target_compile_options(${PROJECT_NAME} PUBLIC -O2 -frounding-math)

    # SIMD only on x86_64
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        if(ENABLE_AVX2)
            target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
        elseif(ENABLE_SSE2)
            target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
        endif()
    endif()

    # Stack size linker option per sistema
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        # macOS linker (ld64)
        target_link_options(${PROJECT_NAME} PUBLIC "-Wl,-stack_size,0x800000")  # 8MB
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        # Linux linker (ld)
        target_link_options(${PROJECT_NAME} PUBLIC "-Wl,-z,stack-size=8388608")  # 8MB
    endif()
endif()
