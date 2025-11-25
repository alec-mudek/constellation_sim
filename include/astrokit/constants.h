#pragma once

#include <numbers>

// this file last generated on 2025-10-11 16:43
// data pulled from C:/astrokit/data/planet_database.json

namespace astrokit 
{

    constexpr double MU_SUN_km3_s2 = 132712000000.0;
    constexpr double AU_km         = 149598000.0;
    constexpr double PI            = std::numbers::pi;
    constexpr double DEG2RAD       = PI / 180.0;
	constexpr double RAD2DEG       = 180.0 / PI;

    struct Planet {
        double MU_km3_s2;
        double R_MEAN_km;
        double R_EQUATOR_km;
        double R_POLE_km;
        double SMA_km;
        double ECC_deg;
        double INC_deg;
        double T_days;
        double J2;
    };

    // Following planetary parameters from the GSFC Planetary Fact Sheets
    inline constexpr Planet MERCURY =
    {
        .MU_km3_s2    = 22032.0,
        .R_MEAN_km    = 2439.7,
        .R_EQUATOR_km = 2440.5,
        .R_POLE_km    = 2438.3,
        .SMA_km       = 57909000.0,
        .ECC_deg      = 0.2056,
        .INC_deg      = 7.004,
        .T_days       = 87.9691,
        .J2           = 5.03e-5
    };

    inline constexpr Planet VENUS =
    {
        .MU_km3_s2    = 324860.0,
        .R_MEAN_km    = 6051.8,
        .R_EQUATOR_km = 6051.8,
        .R_POLE_km    = 6051.8,
        .SMA_km       = 108210000.0,
        .ECC_deg      = 0.0068,
        .INC_deg      = 3.395,
        .T_days       = 224.701,
        .J2           = 4.458e-6
    };                

    inline constexpr Planet EARTH =
    {
        .MU_km3_s2    = 398600.4418,
        .R_MEAN_km    = 6371.0,
        .R_EQUATOR_km = 6378.1,
        .R_POLE_km    = 6356.8,
        .SMA_km       = 149598000.0,
        .ECC_deg      = 0.0167,
        .INC_deg      = 0.0,
        .T_days       = 365.256,
        .J2           = 1.08263e-3
    };
    
    inline constexpr Planet MARS =
    {
        .MU_km3_s2    = 42828.0,
        .R_MEAN_km    = 3389.5,
        .R_EQUATOR_km = 3396.2,
        .R_POLE_km    = 3376.2,
        .SMA_km       = 227956000.0,
        .ECC_deg      = 0.0935,
        .INC_deg      = 1.848,
        .T_days       = 686.98,
        .J2           = 1.96045e-3
    };
    
    inline constexpr Planet JUPITER =
    {
        .MU_km3_s2    = 126687000.0,
        .R_MEAN_km    = 69911.0,
        .R_EQUATOR_km = 71492.0,
        .R_POLE_km    = 66854.0,
        .SMA_km       = 778479000.0,
        .ECC_deg      = 0.0487,
        .INC_deg      = 1.304,
        .T_days       = 4332.589,
        .J2           = 0.014736
    };
    
    inline constexpr Planet SATURN =
    {
        .MU_km3_s2    = 37931000.0,
        .R_MEAN_km    = 58232.0,
        .R_EQUATOR_km = 60268.0,
        .R_POLE_km    = 54364.0,
        .SMA_km       = 1432041000.0,
        .ECC_deg      = 0.052,
        .INC_deg      = 2.486,
        .T_days       = 10755.699,
        .J2           = 0.016298
    };
    
    inline constexpr Planet URANUS =
    {
        .MU_km3_s2    = 5794000.0,
        .R_MEAN_km    = 25362.0,
        .R_EQUATOR_km = 25559.0,
        .R_POLE_km    = 24973.0,
        .SMA_km       = 2867043000.0,
        .ECC_deg      = 0.0469,
        .INC_deg      = 0.77,
        .T_days       = 30685.4,
        .J2           = 0.00334343
    };
    
    inline constexpr Planet NEPTUNE =
    {
        .MU_km3_s2    = 6835100.0,
        .R_MEAN_km    = 24622.0,
        .R_EQUATOR_km = 24764.0,
        .R_POLE_km    = 24341.0,
        .SMA_km       = 4514953000.0,
        .ECC_deg      = 0.0097,
        .INC_deg      = 1.77,
        .T_days       = 60189.018,
        .J2           = 0.003411
    };
    
} // namespace astrokit