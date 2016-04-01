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
	
	double scalar_density;
    double fermi_momentum;
	
} gap_equation_input;

double GapEquationSolver();

double ScalarDensity(double mass, double fermi_momentum);

int WriteZeroedGapEquation(char * filename, double minimum_mass, double maximum_mass, int points_number, double fermi_momentum);

double VacuumMassDetermination();
int WriteVacuumMassEquation(char * filename, double minimum_mass, double maximum_mass, int points_number);

double VacuumThermodynamicPotential(double vacuum_mass, double fermi_momentum);

double ThermodynamicPotential(double mass, double fermi_momentum, double barionic_density, double scalar_density, double vacuum_termodynamic_potential);

double EnergyDensity(double thermodynamic_potential, double chemical_potential, double barionic_density);
double Pressure(double thermodynamic_potential);


#endif /* EOS_h */
