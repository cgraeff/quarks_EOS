//
//  EOS.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"

double GapEquationSolver(double fermi_momentum)
{
    // Prepare input for ZeroedGapEquation
	gap_equation_input input;
    input.fermi_momentum = fermi_momentum;
	
    // Prepare function to be passed to the root finding algorithm
	gsl_function F;
	F.function = &ZeroedGapEquation;
	F.params = &input;

    // Setup root finding solver
    const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection; // Maybe this would be better: gsl_root_fsolver_brent
    
    gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, parameters.gap_eq_solver_lower_bound, parameters.gap_eq_solver_upper_bound);
	
    // Iterate the algorithm until
    // |x_lower - x_upper| < const_abs_error_gap_eq_solving + const_rel_error_gap_eq_solving * MIN(|x_upper|, |x_lower|)
    // or i = gap_eq_solver_max_iterations
	int i = 0;
	double x_lower;
	double x_upper;
	do{
		i++;
		
		int status = gsl_root_fsolver_iterate(s);
		
		if (status != GSL_SUCCESS){
			printf("ERROR: No solution to the gap equation was found!\n");
			exit(EXIT_FAILURE);
		}
		
		x_lower = gsl_root_fsolver_x_lower(s);
		x_upper = gsl_root_fsolver_x_upper(s);
	} while(GSL_CONTINUE == gsl_root_test_interval(x_lower, x_upper, parameters.gap_eq_solver_abs_error, parameters.gap_eq_solver_rel_error)
			&& i <= parameters.gap_eq_solver_max_iterations);

	double result = gsl_root_fsolver_root(s);

	void gsl_root_fsolver_free(gsl_root_fsolver * S);
		
	return result;
}

double ZeroedGapEquation(double mass, void * input)
{
	gap_equation_input * param = (gap_equation_input *)input;
	
	double scalar_density = ScalarDensity(mass, param->fermi_momentum);
	
	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.G_S * scalar_density;

	return mass + gap_1st_term - parameters.bare_mass;
}

double ScalarDensity(double mass, double fermi_momentum)
{
	return NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (mass / pow(M_PI, 2.0))
           * (F0(mass, fermi_momentum) - F0(mass, parameters.cutoff));
}

double F0(double mass, double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
	
	return (1.0 / 2.0) * (momentum * E - pow(mass, 2.0) * log((momentum + E) / mass));
}
