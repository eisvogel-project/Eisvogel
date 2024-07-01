add_executable(testCoulomb testCoulomb.cxx)
target_link_libraries(testCoulomb eisvogel TestUtils)

add_executable(testAskaryan testAskaryan.cxx)
target_link_libraries(testAskaryan eisvogel TestUtils)

add_executable(testCherenkov testCherenkov.cxx)
target_link_libraries(testCherenkov eisvogel-meep TestUtils)
