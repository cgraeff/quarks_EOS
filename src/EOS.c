//
//  EOS.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"

double F0(double mass, double momentum);
double F2(double mass, double momentum);
double F_E(double mass, double momentum);

typedef struct _multi_dim_gap_eq_param {
    double temperature;
    double barionic_density;
} multi_dim_gap_eq_param;

typedef struct _fermi_dirac_distrib_integrand{
    double mass;
    double chemical_potential;
    double temperature;
} fermi_dirac_distrib_integrand;

#define ROOT_FINDER_DIMENSION_2 2

int ZeroedGapAndRenormalizedChemicalPotentialEquations(const gsl_vector * x,
                                                       void * p,
                                                       gsl_vector * f);
void TwodimensionalRootFinder(gsl_multiroot_function * f,
                              double x_initial_guess,
                              double y_initial_guess,
                              double * return_x,
                              double * return_y);

double FermiDiracDistributionIntegral(double mass,
                                      double renormalized_chemical_potential);

double SomethingFermiDiracDistributionIntegral(double mass,
                                               double renormalized_chemical_potential);

double FermiDiracDistributionIntegralIntegrand(double momentum, void * parameters);
double SomethingFermiDiracDistributionIntegralIntegrand(double momentum,
                                                        void * params);

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);
double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);
double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit);



void CalculateMassAndRenormalizedChemicalPotentialSimultaneously(double temperature,
                                                                 double barionic_density,
                                                                 double * return_mass,
                                                                 double * return_renormalized_chemical_potential)
{
    multi_dim_gap_eq_param p;
    p.temperature = temperature;
    p.barionic_density = barionic_density;
    
    gsl_multiroot_function f;
    f.f = &ZeroedGapAndRenormalizedChemicalPotentialEquations;
    f.n = ROOT_FINDER_DIMENSION_2;
    f.params = (void *)&p;
    
    double x;
    double y;
    
    TwodimensionalRootFinder(&f,
                             parameters.mass_and_renor_chem_pot_solution_mass_guess,
                             parameters.mass_and_renor_chem_pot_solution_renor_chem_pot_guess,
                             &x,
                             &y);
    
    *return_mass = x;
    *return_renormalized_chemical_potential = y;
}

void TwodimensionalRootFinder(gsl_multiroot_function * f,
                              double x_initial_guess,
                              double y_initial_guess,
                              double * return_x,
                              double * return_y)
{
 	int status;
  	size_t iter = 0;

  	gsl_vector * initial_guess = gsl_vector_alloc(ROOT_FINDER_DIMENSION_2);

  	gsl_vector_set(initial_guess, 0, x_initial_guess);
  	gsl_vector_set(initial_guess, 1, y_initial_guess);

    const gsl_multiroot_fsolver_type * solver_type = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(solver_type,
                                                                 ROOT_FINDER_DIMENSION_2);
    
  	gsl_multiroot_fsolver_set(solver, f, initial_guess);

  	do {
		iter++;
      	status = gsl_multiroot_fsolver_iterate(solver);

      	//print_state(iter, s);

        if (status == GSL_EBADFUNC){
            printf("TwodimensionalRootFinder: Error: Infinity or division by zero.\n");
            exit(EXIT_FAILURE);
        }
        else if (status == GSL_ENOPROG){
            printf("TwodimensionalRootFinder: Error: Solver is stuck. Try a different initial guess.\n");
            exit(EXIT_FAILURE);
        }

        // Check if the root is good enough:
        // tests for the convergence of the sequence by comparing the last step dx with the
        // absolute error epsabs and relative error epsrel to the current position x. The test
        // returns GSL_SUCCESS if the following condition is achieved,
        //
        // |dx_i| < epsabs + epsrel |x_i|

        gsl_vector * x = gsl_multiroot_fsolver_root(solver); // current root
        gsl_vector * dx = gsl_multiroot_fsolver_dx(solver); // last step
        
      	status = gsl_multiroot_test_delta(dx,
                                          x,
                                          parameters.mass_and_renorm_chem_pot_solution_abs_error,
                                          parameters.mass_and_renorm_chem_pot_solution_rel_error);
        
    } while (status == GSL_CONTINUE
             && iter < parameters.mass_and_renor_chem_pot_solution_max_iter);

    gsl_vector * solution = gsl_multiroot_fsolver_root(solver);
    
  	*return_x = gsl_vector_get(solution, 0);
  	*return_y = gsl_vector_get(solution, 1);

  	gsl_multiroot_fsolver_free(solver);
  	gsl_vector_free (initial_guess);
}

int ZeroedGapAndRenormalizedChemicalPotentialEquations(const gsl_vector * x,
													   void * p,
													   gsl_vector * f)
{
	multi_dim_gap_eq_param * params = (multi_dim_gap_eq_param *)p;

   	const double mass = gsl_vector_get(x,0);
   	const double renormalized_chemical_potential = gsl_vector_get(x,1);

  	double integral = FermiDiracDistributionIntegral(mass,
													 renormalized_chemical_potential);

  	double integral_2 = SomethingFermiDiracDistributionIntegral(mass,
																renormalized_chemical_potential);

  	double zeroed_gap_eq = mass
  						   - parameters.bare_mass
  						   - 4.0 * NUM_COLORS * NUM_FLAVORS * mass * integral_2;
  	double zeroed_chemical_pot_eq = params->barionic_density
  									- 2.0 * NUM_COLORS * NUM_FLAVORS * integral;

   	gsl_vector_set (f, 0, zeroed_gap_eq);
   	gsl_vector_set (f, 1, zeroed_chemical_pot_eq);

    return GSL_SUCCESS;
}

double FermiDiracDistributionIntegral(double mass,
                                      double renormalized_chemical_potential)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &FermiDiracDistributionIntegralIntegrand;
    F.params = &p;
    
    return OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);
}

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit)
{
	gsl_integration_workspace * workspace =
        gsl_integration_workspace_alloc(parameters.mass_and_renor_chem_pot_solution_fermi_dirac_integ_max_interval_num);

    double integral;
    double abserr;
  	gsl_integration_qag(F,
                        lower_limit,
                        upper_limit,
                        parameters.mass_and_renor_chem_pot_solution_fermi_dirac_integ_abs_error,
                        parameters.mass_and_renor_chem_pot_solution_fermi_dirac_integ_rel_error,
                        parameters.mass_and_renor_chem_pot_solution_fermi_dirac_integ_max_sub_interval,
                        parameters.mass_and_renor_chem_pot_solution_fermi_dirac_integ_key,
                        workspace,
                        &integral,
                        &abserr);

	gsl_integration_workspace_free(workspace);

	return integral;
}

double FermiDiracDistributionIntegralIntegrand(double momentum, void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;
    
 	double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

  	double particle_term = FermiDiracDistributionForParticles(E,
															  p->chemical_potential,
															  p->temperature);
  	double antiparticle_term = FermiDiracDistributionForAntiparticles(E,
																	  p->chemical_potential,
																	  p->temperature);
  	return (particle_term - antiparticle_term) * pow(momentum, 2.0);
}

double SomethingFermiDiracDistributionIntegral(double mass,
											   double renormalized_chemical_potential)
{
	// Maybe it's possible to write the first term as the scalar density
	// at zero temperature

    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &SomethingFermiDiracDistributionIntegralIntegrand;
    F.params = &p;
    
    return OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);
}

double SomethingFermiDiracDistributionIntegralIntegrand(double momentum,
														void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;
    
    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double particle_term = FermiDiracDistributionForParticles(E,
                                                              p->chemical_potential,
                                                              p->temperature);
    double antiparticle_term = FermiDiracDistributionForAntiparticles(E,
                                                                      p->chemical_potential,
                                                                      p->temperature);
    return (1.0 + particle_term - antiparticle_term) * pow(momentum, 2.0) / E;
}

double FermiDiracDistributionForParticles(double energy,
										  double chemical_potential,
										  double temperature)
{
	return 1.0 / (1.0 + exp((energy - chemical_potential)/temperature));
}

double FermiDiracDistributionForAntiparticles(double energy,
											  double chemical_potential,
											  double temperature)
{
  	return 1.0 / (1.0 + exp((energy - chemical_potential)/temperature));
}

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
    
    // The 'if' statement does the work of a step function
    double arg = 0;
    if (pow(renor_chem_pot, 2.0) - pow(param->mass, 2.0) >= 0)
        arg = pow(renor_chem_pot, 2.0) - pow(param->mass, 2.0);

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
    
    return first_term + second_term + third_term;
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
