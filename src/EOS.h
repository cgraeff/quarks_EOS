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
double ZeroedGapEquation(double mass, void * input);

double ScalarDensity(double mass, double fermi_momentum);

double F0(double mass, double momentum);

#endif /* EOS_h */
