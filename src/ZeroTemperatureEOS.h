//
//  ZeroTemperatureEOS.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef ZeroTemperatureEOS_h
#define ZeroTemperatureEOS_h

typedef struct _gap_equation_input{
    double fermi_momentum;
} gap_equation_input;

typedef struct _renorm_chem_pot_equation_input{
    double chemical_potential;
    double mass;
} renorm_chem_pot_equation_input;

double UnidimensionalRootFinder(gsl_function * F,
								double lower_bound,
								double upper_bound,
								double abs_error,
								double rel_error,
								int max_iter);
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
