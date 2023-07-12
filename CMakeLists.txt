#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.14)
project(micromorphic_elasto_plasticity)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")

# Add flag for whether to build the python bindings
set(MICROMORPHIC_ELASTO_PLASTICITY_BUILD_PYTHON_BINDINGS ON CACHE BOOL "Flag for whether to build the python bindings")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
endif()

# Get version number from Git
if(${not_conda_test} STREQUAL "true")
    # TODO: On osx-arm64 some c++ projects struggle to find $PREFIX/bin/python during conda-builds and others don't.
    # Trace the source and remove the hint when CMake configures with the correct python version during conda-build.
    # FIXME: troubleshoot sstbigbird/RHEL8 python interpretter mismatch
    # https://re-git.lanl.gov/aea/material-models/micromorphic_tools/-/issues/9
    if (DEFINED ENV{PREFIX})
        set(Python_ROOT_DIR "$ENV{PREFIX}/bin")
    elseif (DEFINED ENV{CONDA_PREFIX})
        set(Python_ROOT_DIR "$ENV{CONDA_PREFIX}/bin")
    endif()
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(COMMAND python -m setuptools_scm
		            OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${${PROJECT_NAME}_VERSION_STRING_FULL} STREQUAL "")
        set(${PROJECT_NAME}_VERSION 0.0.0)
    else()
         string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    endif()
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "error_tools" "vector_tools" "constitutive_tools" "solver_tools" "micromorphic_tools" "micromorphic_linear_elasticity")
if(${not_conda_test} STREQUAL "false")
    # During conda-build testing, we must find the installed project files as if they were an external project
    set(upstream_packages ${upstream_packages} ${PROJECT_NAME})
endif()
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS ${MICROMORPHIC_ELASTO_PLASTICITY_BUILD_PYTHON_BINDINGS})
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # TODO: Decide how to handle main/dev build differences. bash env variables for GIT TAG?
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY ssh://git@re-git.lanl.gov:10022/aea/material-models/${package}.git
            GIT_TAG        origin/dev
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
if(${not_conda_test} STREQUAL "true")
    include_directories(${CPP_SRC_PATH})
    add_subdirectory(${CPP_SRC_PATH})
endif()

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost and add tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    add_subdirectory(${CPP_TEST_PATH})
    # Add docs
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
if(${not_conda_test} STREQUAL "true")
    foreach(package ${PROJECT_NAME} ${UMAT})
        include(CMakePackageConfigHelpers)
        write_basic_package_version_file("${package}ConfigVersion.cmake"
                                         VERSION ${PROJECT_VERSION}
                                         COMPATIBILITY SameMajorVersion)
        configure_package_config_file(
          "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
          "${PROJECT_BINARY_DIR}/${package}Config.cmake"
          INSTALL_DESTINATION
          ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

        # CMake won't build the targets for local builds of upstream projects
        if(cmake_build_type_lower STREQUAL release)
            install(EXPORT ${package}_Targets
                    FILE ${package}Targets.cmake
                    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
        endif()

        install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                      "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
    endforeach(package)
endif()
