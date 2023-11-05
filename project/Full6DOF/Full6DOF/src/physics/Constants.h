#pragma once

namespace CONSTANTS
{
    // SI Units
    constexpr double BOLTZMAN = 1.380649e-23; // Joule/Kelvin , SI exact by definition

    constexpr double AVOGADRO = 6.02214076e23; // # particles , SI exact by definition

    constexpr double E_CHARGE = 1.602176634e-19; // Amp sec , SI exact by definition

    constexpr double PLANCK = 6.62607015e-34; // Joule sec , SI exact by definition

    constexpr double C = 299792458; // meters / sec , SI exact by definition
    constexpr int C_ = 299792458; 

    constexpr double CS_FREQ = 9192631770; // Hertz , SI exact by definition
    constexpr long CS_FREQ_ = 9192631770; 

    constexpr double FREQ_CANDELA = 540e12; // Hertz, monochromatic radiation of frequency for candela defintion
    constexpr long FREQ_CANDELA = 540'000'000'000'000L;
    constexpr double LUM_EFFICIENCY_CANDELA = 686; // lumens  / Watt
    constexpr int LUM_EFFICIENCY_CANDELA = 686;

    // Other Exact Numbers
    constexpr double PI = 3.1415926535897932384626433832795;
    constexpr double TWOPI = 6.283185307179586476925286766559;
    constexpr double FOURPI = 12.566370614359172953850573533118;
    constexpr double HALFPI = 1.5707963267948966192313216916398;
    constexpr double QUARTERPI = 0.78539816339744830961566084581988;
    constexpr double LN2 = 0.69314718055994530941723212145818;
    constexpr double LOG2 = 0.30102999566398119521373889472449;
    constexpr double E = 2.7182818284590452353602874713527;
    constexpr double SQRT2 = 1.4142135623730950488016887242097;
    constexpr double SQRT3 = 1.7320508075688772935274463415059;

    constexpr double C_SQ = 89875517873681764;
    constexpr long C_SQ_ = 89875517873681764L;

    // Uncertain Numbers
    constexpr double VACUUM_PERMEABILITY = 1.2566370621219e-6; // 4*pi*1e-7 Newton / Ampere^2
    constexpr double VACUUM_PERMITTIVITY = 8.854187812813e-12; // Farad / meter
    constexpr double COULOMB = 8.987551792314e9; // Newton meters^2 / coulomb^2

    constexpr double GRAVITATIONAL = 6.6743015e-11; // Newton meters^2 / kilogram
    constexpr double GRAVITATIONAL_RELATIVE_UNCERTAINTY = 2.2e-5;

    constexpr double FINE_STRUCTURE_CONSTANT = 0.007297352569311;

    // same for VACUUM_PERMEABILITY and VACUUM_PERMITTIVITY and COULOMB
    constexpr double FINE_STRUCTURE_CONSTANT_RELATIVE_UNCERTAINTY = 1.5e-10; 

    constexpr double ABSOLUTE_ZERO_KELVIN = 273.16001;
    constexpr double ABSOLUTE_ZERO_KELVIN_UNCERTAINTY = 1e-4;

}