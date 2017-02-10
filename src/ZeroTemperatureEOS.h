//
//  ZeroTemperatureEOS.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef ZeroTemperatureEOS_h
#define ZeroTemperatureEOS_h


/* Parameters for root finding in the solution of the Gap Equation, Vacuum Mass,
   and Renormalized Chemical Potential:
       - The root is assumed to be in the region bracketed by *_lower_bound and *_upper_bound;
 		- Two kinds of error are defined:
            - *_const_abs_error: An absolute error in MeV for the width of the region
                                 which contais the root
            - *_const_rel_error: Error as a fraction of the value of MIN(x_lower, x_upper),
                                 where x_lower and x_upper are the values that bracket the
                                 root after each iteration of the root findind algorithm
       - The root is accepted when
            |x_lower - x_upper| < *_const_abs_error
                                  + *_const_rel_error * MIN(|x_upper|, |x_lower|)
*/
typedef struct _ZeroTemperatureParameters{
	int max_iterations;
    double lower_bound;
    double upper_bound;
    double abs_error;
    double rel_error;
} ZeroTemperatureParameters;

typedef struct _VacuumMassDeterminationParameters{
   	int max_iterations;
   	double lower_bound;
   	double upper_bound;
   	double abs_error;
   	double rel_error;
} VacuumMassDeterminationParameters;

typedef struct _gap_equation_input{
    double fermi_momentum;
} gap_equation_input;

typedef struct _renorm_chem_pot_equation_input{
    double chemical_potential;
    double mass;
} renorm_chem_pot_equation_input;

double GapEquationSolver();
double ZeroedGapEquation(double mass, void * input);

double ScalarDensity(double mass, double fermi_momentum);

double VacuumMassDetermination();
double VacuumMassEquation(double mass, void * input);

double ZeroedRenormalizedChemicalPotentialEquation(double renor_chem_pot,
                                                   void * input);

double ThermodynamicPotential(double mass,
                              double fermi_momentum,
                              double chemical_potential,
                              double renormalized_chemical_potential);

double EnergyDensity(double regularized_thermodynamic_potential,
                     double chemical_potential,
                     double barionic_density);

double Pressure(double regularized_thermodynamic_potential);

#endif /* ZeroTemperatureEOS_h */
