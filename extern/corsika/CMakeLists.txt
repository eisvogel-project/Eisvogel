add_library(eisvogel-corsika SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/C8.cxx)
target_include_directories(eisvogel-corsika
  PUBLIC
  ./include/
  ./impl/
)

if(BUILD_MEEP)

  target_compile_definitions(eisvogel-corsika PRIVATE WITHMEEP=1)
  target_sources(
    eisvogel-corsika
    PUBLIC FILE_SET HEADERS
    FILES
    ./include/Eisvogel/C8SignalCalculator.hh
    ./include/Eisvogel/C8GreensFunctionCalculator.hh
  )
  target_link_libraries(eisvogel-corsika eisvogel-meep)
  
else()

  target_sources(
    eisvogel-corsika
    PUBLIC FILE_SET HEADERS
    FILES
    ./include/Eisvogel/C8SignalCalculator.hh
  )  
  target_link_libraries(eisvogel-corsika eisvogel)
  
endif(BUILD_MEEP)

message("Eisvogel will be installed at '${CMAKE_INSTALL_PREFIX}'")

install(TARGETS eisvogel-corsika
  LIBRARY DESTINATION ./${CMAKE_INSTALL_LIBDIR}
  FILE_SET HEADERS
  DESTINATION ./${CMAKE_INSTALL_INCLUDEDIR})

configure_file(eisvogel-corsika.pc.in eisvogel-corsika.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/extern/corsika/eisvogel-corsika.pc DESTINATION ./${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
