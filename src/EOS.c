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

double UnidimensionalRootFinder(gsl_function * F, double lower_bound, double upper_bound, double abs_error, double rel_error, int max_iter);
double F0(double mass, double momentum);
double F2(double mass, double momentum);
double F_E(double mass, double momentum);
double ZeroedGapEquation(double mass, void * input);
double VacuumMassEquation(double mass, void * input);


double UnidimensionalRootFinder(gsl_function * F, double lower_bound, double upper_bound, double abs_error, double rel_error, int max_iter)
{
    // Setup root finding solver
    const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection; // Maybe this would be better: gsl_root_fsolver_brent
    
    gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, F, lower_bound, upper_bound);
	
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
	} while(GSL_CONTINUE == gsl_root_test_interval(x_lower, x_upper, abs_error, rel_error) && i <=  max_iter);

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
    
    //    printf("fermi momentum: %f\t scalar density: %f\n", param->fermi_momentum, scalar_density);
	
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
    double term = 2.0 * NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -2.0) * parameters.G_S * mass * F_diff / pow(M_PI, 2.0);
    
    return mass - parameters.bare_mass - term;
}

int WriteVacuumMassEquation(char * filename, double minimum_mass, double maximum_mass, int points_number){
    
    double m = 0;
    
    double step = (maximum_mass - minimum_mass) / (points_number - 1);
    
    FILE * f = fopen(filename, "w");
    
    if (NULL == f) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    while (m < points_number) {
        fprintf(f, "%20.15E\t%20.15E\n", m, VacuumMassEquation(m, NULL));
        m += step;
    }
    
    fclose(f);
    return 0;
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

int WriteZeroedGapEquation(char * filename, double minimum_mass, double maximum_mass, int points_number, double fermi_momentum){
    
    double m = 0;
    
    double step = (maximum_mass - minimum_mass) / (points_number - 1);
    
    FILE * f = fopen(filename, "w");
    
    if (NULL == f) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    gap_equation_input input;
    input.fermi_momentum = fermi_momentum;
    
    while (m < maximum_mass) {
        fprintf(f, "%20.15E\t%20.15E\n", m, ZeroedGapEquation(m, &input));
        m += step;
    }
    
    fclose(f);
    return 0;
}

double ThermodynamicPotential2(double mass,
                              double barionic_density,
                              double fermi_momentum,
                              double scalar_density,
                              double chemical_potential,
                              double vacuum_thermodynamic_potential)
{
    double kinectic = - NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (F2(mass, parameters.cutoff) - F2(mass, fermi_momentum)) / pow(M_PI, 2.0);
    double first_term = parameters.bare_mass * scalar_density;
    double second_term = - CONST_HBAR_C * parameters.G_S * pow(scalar_density, 2.0);
    double third_term = - NUM_COLORS * chemical_potential * barionic_density;
    
    return 2.0 * kinectic + first_term + second_term + third_term - vacuum_thermodynamic_potential;
}

double VacuumThermodynamicPotential2(double vacuum_mass, double barionic_density, double chemical_potential)
{
    double vacuum_scalar_density = NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (vacuum_mass / pow(M_PI, 2.0))
                                   * (F0(vacuum_mass, 0.0) - F0(vacuum_mass, parameters.cutoff));
    
    double kinectic = - NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (F2(vacuum_mass, parameters.cutoff) - F2(vacuum_mass, 0.0)) / pow(M_PI, 2.0);
    double first_term = parameters.bare_mass * vacuum_scalar_density;
    double second_term = - CONST_HBAR_C * parameters.G_S * pow(vacuum_scalar_density, 2.0);
    
    return 2.0 * kinectic + first_term + second_term;
}

double F2(double mass, double momentum)
{
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
    
    return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum + 2.0 * pow(momentum, 3.0)) * E
    + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}


double VacuumThermodynamicPotential(double vacuum_mass, double fermi_momentum)
{
    double F_E_diff = F_E(vacuum_mass, parameters.cutoff) - F_E(vacuum_mass, fermi_momentum);
    
    double first_term = - NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * F_E_diff / pow(M_PI, 2.0);
    double second_term = pow(vacuum_mass - parameters.bare_mass, 2.0) / (4.0 * parameters.G_S * CONST_HBAR_C);
    
    return first_term + second_term;
}

double ThermodynamicPotential(double mass,
                              double fermi_momentum,
                              double barionic_density,
                              double chemical_potential,
                              double vacuum_thermodynamic_potential)

{
    double F_diff = F_E(mass, parameters.cutoff) - F_E(mass, fermi_momentum); // From Buballa's Phys. Rep. 407 (2005), p. 233, 234, I guess that this should be from zero to cutoff, not from fermi momentum to cutoff (it's an integral, so the upper limit goes before the sign)
    
    double first_term = - NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (F_diff + chemical_potential * pow(fermi_momentum, 3.0)) / pow(M_PI, 2.0);
    double second_term = pow(mass - parameters.bare_mass, 2.0) / (4.0 * parameters.G_S * CONST_HBAR_C);
    
    return first_term + second_term - vacuum_thermodynamic_potential;
}

double F_E(double mass, double momentum)
{
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
    
    return (momentum * pow(E, 3.0) - 0.5 * pow(mass, 2.0) * momentum * E - 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass)) / 4.0;
}

double EnergyDensity(double thermodynamic_potential, double chemical_potential, double barionic_density)
{
    return thermodynamic_potential + NUM_COLORS * chemical_potential * barionic_density;
}

double Pressure(double thermodynamic_potential)
{
    return -thermodynamic_potential;
}