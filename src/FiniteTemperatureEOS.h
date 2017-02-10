//
//  FiniteTemperatureEOS.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef FiniteTemperatureEOS_h
#define FiniteTemperatureEOS_h

typedef struct _EntropyParameters{
    double lower_limit;
    double upper_limit;
    double abs_error;
    double rel_error;
    int max_sub_interval;
    int integration_key;
} EntropyParameters;

// Parameters for simultaneous solution of mass and renormalized
// chemical potential at finite temperature
typedef struct _FiniteTemperatureParameters{
   	double temperature; // (MeV)

   	int max_iter;
   	double mass_guess; // (MeV)
   	double renor_chem_pot_guess; //(MeV)
   	double abs_error;
   	double rel_error;
   	double renor_chem_pot_lower_bound;
   	double renor_chem_pot_upper_bound;

	EntropyParameters entropy;

} FiniteTemperatureParameters;

typedef struct _RenomalizedChemicalPotentialIntegrationParameters{
   	int max_iter;
   	double lower_bound;
   	double upper_bound;
   	double abs_error;
   	double rel_error;
} RenormalizedChemicalPotentialIntegrationParameters;

typedef struct _multi_dim_gap_eq_param {
    double barionic_density;
} multi_dim_gap_eq_param;

typedef struct _therm_pot_free_gas_contrib_params{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} therm_pot_free_gas_contrib_params;

typedef struct _entropy_integrand_parameters{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} entropy_integrand_parameters;

void SolveMultiRoots(double  barionic_density,
                     double *return_mass,
                     double *return_renorm_chem_pot);

double ZeroedGapEquationForFiniteTemperature(double mass,
                                             double renormalized_chemical_potential,
                                             void * params);

double ZeroedBarionicDensityEquationForFiniteTemperature(double mass,
                                                         double renormalized_chemical_potential,
                                                         void * barionic_density);

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

// FIXME: Remove once I'm sure the other version works fine
double EntropyIntegrandFromDerivative(double momentum, void * parameters);
double EntropyIntegrandArt(double momentum, void * parameters);



#endif /* FiniteTemperatureEOS_h */
