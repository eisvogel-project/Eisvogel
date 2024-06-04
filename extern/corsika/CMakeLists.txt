add_library(eisvogel-corsika SHARED
            src/C8SignalCalculator.cxx)
set_target_properties(eisvogel-corsika PROPERTIES PUBLIC_HEADER include/C8SignalCalculator.hh)

target_include_directories(eisvogel-corsika PUBLIC include/)
target_link_libraries(eisvogel-corsika eisvogel)

configure_file(eisvogel-corsika.pc.in eisvogel-corsika.pc @ONLY)

message(CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

install(FILES ${CMAKE_BINARY_DIR}/extern/corsika/eisvogel-corsika.pc DESTINATION ./${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
install(TARGETS eisvogel-corsika
        LIBRARY DESTINATION ./${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ./${CMAKE_INSTALL_INCLUDEDIR}/Eisvogel)

add_subdirectory(tests)
