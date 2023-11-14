if(BUILD_WFC)
    add_executable(testHDF5 testHDF5.cxx)
    target_link_libraries(testHDF5 eisvogel_wfc)
endif(BUILD_WFC)
