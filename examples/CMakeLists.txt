set(CMAKE_CXX_FLAGS "-O3 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -fno-trapping-math -ftree-vectorizer-verbose=2 -march=native -save-temps -fverbose-asm -g")

if(BUILD_WFC)
    add_subdirectory(dipole_free_space)
    add_subdirectory(dipole_ice)
endif(BUILD_WFC)

add_subdirectory(coulomb)
add_subdirectory(point_charge_coulomb)
add_subdirectory(point_charge_cherenkov)
add_subdirectory(point_charge_askaryan)
add_subdirectory(shower_profile)
add_subdirectory(impacting_shower)
add_subdirectory(plot_greens_function)
