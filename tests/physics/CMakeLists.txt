set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/tests/physics/data)

add_executable(testCoulomb testCoulomb.cxx)
target_link_libraries(testCoulomb eisvogel TestUtils)
add_test(NAME testCoulomb COMMAND testCoulomb ${TEST_DATA_DIR}/testCoulomb)

add_executable(testAskaryan testAskaryan.cxx)
target_link_libraries(testAskaryan eisvogel TestUtils)
add_test(NAME testAskaryan COMMAND testAskaryan ${TEST_DATA_DIR}/testAskaryan)

if(BUILD_MEEP)
    add_executable(testCherenkov testCherenkov.cxx)
    target_link_libraries(testCherenkov eisvogel-meep TestUtils)
    # add_test(NAME testCherenkov COMMAND testCherenkov)
endif(BUILD_MEEP)
