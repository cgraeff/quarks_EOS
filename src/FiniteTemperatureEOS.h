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

typedef struct _entropy_integrand_parameters{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} entropy_integrand_parameters;

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

// FIXME: Remove once I'm sure the other version works fine
//double EntropyIntegrandFromDerivative(double momentum, void * parameters);
//double EntropyIntegrandArt(double momentum, void * parameters);

// The following function takes another function and verifies whether its absolute return value is
// less than a given tolerance, that is, whether |f(x,y)| <= tolerance is true
// INPUT:
//      double (*f)(double x, double y): a pointer to a function which takes two doubles and return a double
//      min_*, max_*, *_num_pts: minimum and maximum values of each variable, as well as the number of points
//                               in which the function will be probed in a given axis
//      tolerance: the maximum value the function may return (in absolute numbers) for the point to be
//                 registered
//      output_*: Return the values of x and y for each point that satisfies the condition |f(x,y)| <= tolerance
//                Those vectors must be pre-allocated and its sizes must be x_num_pts * y_num_pts (the vectors
//                must be able to hold all points). The i-th ordinate pair (x, y) is given by
//                (x, y)_i = (gsl_vector_get(output_x, i), gsl_vector_get(output_y, i))
//      num_points: returns the number of points that satisfies |f(x,y)| <= tolerance
//      show_progress: The execution of the map procedure may take a while if the map is too large
//                     (too many points), or if the provided function takes too long to calculate.
//                     In these ocasions, percentual may be shown to show progress.
void MapFunction(double (*f)(double, double),
                 double min_x,
                 double max_x,
                 int x_num_pts,
                 double min_y,
                 double max_y,
                 int y_num_pts,
                 double tolerance,
                 gsl_vector * output_x,
                 gsl_vector * output_y,
                 int * num_points,
                 bool show_progress);

#endif /* FiniteTemperatureEOS_h */
