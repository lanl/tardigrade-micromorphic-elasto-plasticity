set(TEST_NAME "test_${PROJECT_NAME}")
add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_target_properties(${TEST_NAME} PROPERTIES CXX_STANDARD 11)
target_compile_options(${TEST_NAME} PRIVATE "-lrt")
target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${TEST_NAME} PUBLIC
                               ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                               ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${tardigrade_micromorphic_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${tardigrade_micromorphic_linear_elasticity_SOURCE_DIR}/${CPP_SRC_PATH})
endif()
