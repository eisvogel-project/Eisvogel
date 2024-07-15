if(BUILD_EXTERN)
    add_subdirectory(shower_profile)
    add_subdirectory(corsika)
endif(BUILD_EXTERN)

if(BUILD_CORSIKA)
    add_subdirectory(corsika)
endif(BUILD_CORSIKA)
