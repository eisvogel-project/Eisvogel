find_package(Python3 COMPONENTS Interpreter)

find_package(MPI REQUIRED)

set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(PY_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/")

configure_file(${SETUP_PY_IN} ${SETUP_PY})

add_custom_command(OUTPUT "${PY_OUTPUT}" 
                   COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build_ext --inplace
)

add_custom_target(pyeisvogel ALL DEPENDS ${PY_OUTPUT})
add_dependencies(pyeisvogel eisvogel)
