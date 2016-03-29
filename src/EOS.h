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

double TermodynamicPotential();
double EnergyDensity();
double Pressure();


#endif /* EOS_h */
