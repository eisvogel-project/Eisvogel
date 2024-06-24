[![CI](https://github.com/philippwindischhofer/Eisvogel/actions/workflows/build-ci.yml/badge.svg)](https://github.com/philippwindischhofer/Eisvogel/actions/workflows/build-ci.yml)

# Eisvogel

**Eisvogel** is a tool to calculate the time-domain antenna signal in a radio neutrino experiment.
It uses an electrodynamic Green's function to efficiently propagate radiation through complex environments.

## Key features

- All electrodynamic effects are included: no ray-tracing or geometric optics are assumed
- Supports arbitrary linear, inhomogeneous, anisotropic materials
- Interfaces to CORSIKA 8 to simulate signals from microscopic showers

## Contact us

Eisvogel is still under development. If you want to help, or just learn more about the project, feel free to contact us.

## Installation

```
mkdir build && cd build
cmake -DBUILD_TESTS=ON -DBUILD_EXAMPLES=ON ..
make -j5
```

To also build the components for the calculation of the Green's function do (this requires MEEP to be installed and available)

```
cmake -DBUILD_MEEP=ON ..
```

## More information

- P. Windischhofer, C. Welling, C. Deaconu, "Eisvogel: Exact and efficient calculations of radio emissions from in-ice neutrino showers", PoS(ICRC2023)
- P. Windischhofer, W. Riegler, "Electrical signals induced in detectors by cosmic rays: a reciprocal look at electrodynamics", PoS(ICRC2021), 184 (2021)
- W. Riegler, P. Windischhofer, "Signals induced on electrodes by moving charges, a general theorem for Maxwell’s equations based on Lorentz-reciprocity", Nucl. Instrum. Meth. A, 980, 164471 (2020)
