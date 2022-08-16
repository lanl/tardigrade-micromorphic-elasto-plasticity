# Project
add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11
                                                 PUBLIC_HEADER ${PROJECT_NAME}.h)
target_link_libraries(${PROJECT_NAME} error_tools solver_tools micromorphic_tools micromorphic_linear_elasticity)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${PROJECT_NAME} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${solver_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${micromorphic_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${micromorphic_linear_elasticity_SOURCE_DIR}/${CPP_SRC_PATH}")
endif()

foreach(package ${PROJECT_NAME} ${UMAT})
    install(TARGETS ${package}
            EXPORT ${package}_Targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endforeach(package)
