set(SETUP_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in")
set(SETUP_OUT "${PROJECT_BINARY_DIR}/setup.sh")

configure_file(${SETUP_IN} ${SETUP_OUT})
