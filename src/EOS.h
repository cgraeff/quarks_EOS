//
//  EOS.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef EOS_h
#define EOS_h

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

double ScalarDensity(double mass, double fermi_momentum);

int WriteZeroedGapEquation(char * filename, double minimum_mass, double maximum_mass, int points_number, double fermi_momentum);

double VacuumMassDetermination();
int WriteVacuumMassEquation(char * filename, double minimum_mass, double maximum_mass, int points_number);

double ThermodynamicPotential(double mass,
                              double fermi_momentum,
                              double chemical_potential,
                              double renormalized_chemical_potential);
/*
double VacuumThermodynamicPotential2(double vacuum_mass, double barionic_density, double chemical_potential);
double ThermodynamicPotential2(double mass,
                              double barionic_density,
                              double fermi_momentum,
                              double scalar_density,
                              double chemical_potential,
                              double vacuum_thermodynamic_potential);
*/
double EnergyDensity(double regularized_thermodynamic_potential, double chemical_potential, double barionic_density);
double Pressure(double regularized_thermodynamic_potential);

double F0(double mass, double momentum);


#endif /* EOS_h */
