# constellation\_sim v0.1

Basic modeling and simulation environment for a Walker-Delta constellation about a planetary body.



# Current Status

All necessary classes are in place for spacecraft propagation around a planetary body with Keplerian + J2 force modeling. Classes are also in place that define and initialize a Walker-Delta constellation. Code currently establishes a 12/3/1 Walker-Delta constellation at Earth with a 56 degree inclination and 25,000 km semi-major axis.



Next steps:

1. Data product generation and result plotting. Data will include ground station access across the constellation as well as inter-satellite geometries.
2. Station keeping.



# Dependencies

* astrokit: A header-only library with basic astrodynamics functions. Comes in the include/ directory in this repo so there's no need to download it separately. If needed, however, it can be found here: https://github.com/alec-mudek/astrokit/
* cspice: Spice utility that must be downloaded by the user. It can be found here: https://naif.jpl.nasa.gov/naif/toolkit\_C.html
* Eigen: A header-only linear algebra library that must also be downloaded by the user. It can be found here: https://libeigen.gitlab.io/



Before building constellation\_sim, make sure to fill in the correct paths for EIGEN\_INCLUDE\_DIR, CSPICE\_INCLUDE\_DIR, CSPICE\_LIB\_DIR, and ASTROKIT\_INCLUDE\_DIR in CMakeLists\_template.txt.

