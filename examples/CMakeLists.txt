add_subdirectory(dipole_free_space)
add_subdirectory(point_charge_coulomb)
add_subdirectory(point_charge_askaryan)
add_subdirectory(plot_greens_function)

if(BUILD_MEEP)
    add_subdirectory(dipole_ice)
endif(BUILD_MEEP)

if(BUILD_EXTERN)
    add_subdirectory(shower_profile)
endif(BUILD_EXTERN)

