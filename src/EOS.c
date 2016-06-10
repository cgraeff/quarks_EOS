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

double F0(double mass, double momentum);
double F2(double mass, double momentum);
double F_E(double mass, double momentum);


double UnidimensionalRootFinder(gsl_function * F,
								double lower_bound,
								double upper_bound,
								double abs_error,
								double rel_error,
								int max_iter)
{
    // Setup root finding solver. Maybe this would be better: gsl_root_fsolver_brent
    const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection;

    gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, F, lower_bound, upper_bound);
	
    // Iterate the algorithm until
    // |x_lower - x_upper| < const_abs_error_gap_eq_solving
	// 						 + const_rel_error_gap_eq_solving * MIN(|x_upper|, |x_lower|)
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
	} while(GSL_CONTINUE == gsl_root_test_interval(x_lower,
												   x_upper,
												   abs_error,
												   rel_error)
			&& i <=  max_iter);

	double result = gsl_root_fsolver_root(s);

	void gsl_root_fsolver_free(gsl_root_fsolver * S);
		
	return result;
}

double GapEquationSolver(double fermi_momentum)
{
    // Prepare input for ZeroedGapEquation
    gap_equation_input input;
    input.fermi_momentum = fermi_momentum;
    
    // Prepare function to be passed to the root finding algorithm
    gsl_function F;
    F.function = &ZeroedGapEquation;
    F.params = &input;
    
    double root = UnidimensionalRootFinder(&F,
                                           parameters.gap_eq_solver_lower_bound,
                                           parameters.gap_eq_solver_upper_bound,
                                           parameters.gap_eq_solver_abs_error,
                                           parameters.gap_eq_solver_rel_error,
                                           parameters.gap_eq_solver_max_iterations);
    
    return root;
}

double ZeroedGapEquation(double mass, void * input)
{
	gap_equation_input * param = (gap_equation_input *)input;
	
	double scalar_density = ScalarDensity(mass, param->fermi_momentum);

	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.G_S * scalar_density;

	return mass + gap_1st_term - parameters.bare_mass;
}

double VacuumMassDetermination()
{
    // Prepare function to be passed to the root finding algorithm.
    // No parameters are needed.
    gsl_function F;
    F.function = &VacuumMassEquation;
    
    double root = UnidimensionalRootFinder(&F,
                                           parameters.vac_mass_det_lower_bound,
                                           parameters.vac_mass_det_upper_bound,
                                           parameters.vac_mass_det_abs_error,
                                           parameters.vac_mass_det_rel_error,
                                           parameters.vac_mass_det_max_iterations);
    
    return root;
}

double VacuumMassEquation(double mass, void * input)
{
    double F_diff = F0(mass, parameters.cutoff) - F0(mass, 0.0);
    double term = 2.0 * NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -2.0)
  				  * parameters.G_S * mass * F_diff
  				  / pow(M_PI, 2.0);
    
    return mass - parameters.bare_mass - term;
}

double ZeroedRenormalizedChemicalPotentialEquation(double renor_chem_pot,
												   void * input)
{
  	renorm_chem_pot_equation_input * param = (renorm_chem_pot_equation_input *) input;

  	double c = 2.0 * parameters.G_V * NUM_COLORS * NUM_FLAVORS
  				   / (3.0 * pow(M_PI * CONST_HBAR_C, 2.0));
  	double arg = pow(renor_chem_pot, 2.0) - pow(param->mass, 2.0);

  	return renor_chem_pot - param->chemical_potential + c * pow(arg, 3.0 / 2.0);
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


double F2(double mass, double momentum)
{
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
    
    return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum + 2.0 * pow(momentum, 3.0)) * E
    + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}

double ThermodynamicPotential(double mass,
                              double fermi_momentum,
                              double chemical_potential,
                              double renormalized_chemical_potential)

{
    double F_diff = F_E(mass, parameters.cutoff) - F_E(mass, fermi_momentum);
    
    double first_term = - NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0)
                        * (F_diff + chemical_potential * pow(fermi_momentum, 3.0) / 3.0)
                        / pow(M_PI, 2.0);
    double second_term = pow(mass - parameters.bare_mass, 2.0)
  						 / (4.0 * parameters.G_S * CONST_HBAR_C);
    
    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0)
                     / (4.0 * parameters.G_V * CONST_HBAR_C);
    
    return first_term + second_term - third_term;
}

double F_E(double mass, double momentum)
{
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
    
    return (momentum * pow(E, 3.0)
  			- 0.5 * pow(mass, 2.0) * momentum * E
			- 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass))
  		   / 4.0;
}

double EnergyDensity(double regularized_thermodynamic_potential,
					 double chemical_potential,
					 double barionic_density)
{
    return regularized_thermodynamic_potential
  		   + NUM_COLORS * chemical_potential * barionic_density;
}

double Pressure(double regularized_thermodynamic_potential)
{
    return -regularized_thermodynamic_potential;
}
