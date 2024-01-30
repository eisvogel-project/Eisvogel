if(BUILD_WFC)
    add_subdirectory(dipole_free_space)
    add_subdirectory(dipole_ice)
endif(BUILD_WFC)

add_subdirectory(shower_profile)
add_subdirectory(electric_dipole_antenna)

add_subdirectory(coulomb)
add_subdirectory(point_charge_coulomb)
add_subdirectory(point_charge_cherenkov)
add_subdirectory(point_charge_askaryan)
