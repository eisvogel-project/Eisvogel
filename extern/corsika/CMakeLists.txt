add_library(eisvogel-corsika SHARED src/C8.cxx)

if(BUILD_MEEP)

  target_compile_definitions(eisvogel-corsika PRIVATE WITHMEEP=1)  
  set_target_properties(eisvogel-corsika PROPERTIES
    PUBLIC_HEADER include/C8SignalCalculator.hh
    PUBLIC_HEADER include/C8GreensFunctionCalculator.hh)
  target_link_libraries(eisvogel-corsika eisvogel-meep)
  
else()

  set_target_properties(eisvogel-corsika PROPERTIES
    PUBLIC_HEADER include/C8SignalCalculator.hh)
  target_link_libraries(eisvogel-corsika eisvogel)
  
endif(BUILD_MEEP)

target_include_directories(eisvogel-corsika
  PUBLIC include/ impl/)

configure_file(eisvogel-corsika.pc.in eisvogel-corsika.pc @ONLY)

message(CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

install(FILES ${CMAKE_BINARY_DIR}/extern/corsika/eisvogel-corsika.pc DESTINATION ./${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
install(TARGETS eisvogel-corsika
        LIBRARY DESTINATION ./${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ./${CMAKE_INSTALL_INCLUDEDIR}/Eisvogel)
