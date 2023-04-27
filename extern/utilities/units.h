#ifndef __UNITS
#define __UNITS

namespace units {
    // base units
    const double second = 1e+9;
    const double meter = 1;
    const double electronvolt = 1;
    const double e_charge = 1;
    
    // prefix
    const double exa = 1e+18;
    const double peta = 1e+15;
    const double tera = 1e+12;
    const double giga = 1e+9;
    const double mega = 1e+6;
    const double kilo = 1e+3;
    const double deci = 1e-1;
    const double centi = 1e-2;
    const double milli = 1e-3;
    const double micro = 1e-6;
    const double nano = 1e-9;
    const double pico = 1e-12;
    const double femto = 1e-15;
    
    // time units
    
    const double nanosecond = nano * second;
    const double ns = nanosecond;
    const double microsecond = second * micro;
    
    const double hertz = 1 / second;
    const double Hz = hertz;
    const double kilohertz = kilo * hertz;
    const double kHz = kilohertz;
    const double megahertz = mega * hertz;
    const double MHz = megahertz;
    const double gigahertz = giga * hertz;
    const double GHz = gigahertz;
    
    //distance units
    
    const double kilometer = kilo * meter;
    const double km = kilometer;
    const double decimeter = deci * meter;
    const double dm = decimeter;
    const double centimeter = centi * meter;
    const double cm = centimeter;
    const double millimeter = milli * meter;
    const double mm = millimeter;
    const double micrometer = micro * meter;
    const double nanometer = nano * meter;
    const double nm = nanometer;
    
    // volume
    
    const double square_meter = meter * meter;
    const double cubic_meter = meter * meter * meter;
    const double liter = dm * dm * dm;
    const double L = liter;
    const double milliliter = milli * liter;
    const double mL = milliliter;
    
    
    //energy & mass units
    const double eV = electronvolt;
    const double keV = kilo * eV;
    const double MeV = mega * eV;
    const double GeV = giga * eV;
    const double TeV = tera * eV;
    const double PeV = peta * eV;
    const double EeV = exa * eV;
    const double joule = 0.624150974e19 * electronvolt;
    const double J = joule;
    const double kilogram = joule * second * second / meter / meter;
    const double kg = kilogram;
    const double gram = 1e-3 * kilogram;
    const double g = gram;
    const double ton = 1e3 * kg;
    
    // charge & electrical units
    const double coulomb = 0.624150974 * e_charge;
    const double C = coulomb;
    const double ampere = coulomb / second;
    const double A = ampere;
    const double volt = joule / coulomb;
    const double V = volt;
}

#endif