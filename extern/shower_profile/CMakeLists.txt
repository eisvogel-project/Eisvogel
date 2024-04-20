set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_library(ShowerCreator SHARED
            shower_creator/shower_creator.cpp
	    shower_1D/shower_1D.cpp
	    environment/ice_profile.cpp
	    )

target_include_directories(ShowerCreator PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_1D"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/environment"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/charge_excess_profile"
                          "${PROJECT_SOURCE_DIR}/extern/utilities"
			  "${PROJECT_SOURCE_DIR}/src/"
                          )

target_link_libraries(ShowerCreator eisvogel)
