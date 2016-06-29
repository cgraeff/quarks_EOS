//
//  FiniteTemperatureEOS.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef FiniteTemperatureEOS_h
#define FiniteTemperatureEOS_h

typedef struct _multi_dim_gap_eq_param {
    double barionic_density;
} multi_dim_gap_eq_param;

typedef struct _fermi_dirac_distrib_integrand{
    double mass;
    double chemical_potential;
    double temperature;
} fermi_dirac_distrib_integrand;

typedef struct _therm_pot_free_gas_contrib_params{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} therm_pot_free_gas_contrib_params;

double FermiDiracDistributionFromDensityIntegral(double mass,
                                                 double renormalized_chemical_potential);
double FermiDiracDistributionIntegralFromGapEquation(double mass,
                                                     double renormalized_chemical_potential);
void CalculateMassAndRenormalizedChemicalPotentialSimultaneously(double barionic_density,
                                                                 double * return_mass,
                                                                 double * return_renormalized_chemical_potential);

double ZeroedGapEquationForFiniteTemperature(double mass, double renormalized_chemical_potential);

double ZeroedBarionicDensityEquationForFiniteDensity(double mass,
                                                     double renormalized_chemical_potential,
                                                     double barionic_density);

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);
double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);

double ThermodynamicPotentialForFiniteTemperature(double mass,
                                                  double chemical_potential,
                                                  double renormalized_chemical_potential);

double ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(double mass,
                                                                     double renormalized_chemical_potential);

double PressureForFiniteTemperature(double regularized_thermodynamic_potential);
double EnergyForFiniteTemperature(double regularized_thermodynamic_potential,
                                  double chemical_potential,
                                  double barionic_density,
                                  double temperature,
                                  double entropy);

double Entropy(double mass, double temperature, double renormalized_chemical_potential);
double EntropyIntegrand(double momentum, void * parameters);


#endif /* FiniteTemperatureEOS_h */
