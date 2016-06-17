//
//  FiniteTemperatureEOS.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multiroots.h>

#include "Constants.h"
#include "Parameters.h"

#include "FiniteTemperatureEOS.h"

#define ROOT_FINDER_DIMENSION_2 2

int ZeroedGapAndBarionicDensityEquations(const gsl_vector * x,
                                                       void * p,
                                                       gsl_vector * f);
void TwodimensionalRootFinder(gsl_multiroot_function * f,
                              double x_initial_guess,
                              double y_initial_guess,
                              double * return_x,
                              double * return_y);

double FermiDiracDistributionFromDensityIntegralIntegrand(double momentum, void * parameters);
double FermiDiracDistributionIntegralFromGapEquationIntegrand(double momentum,
                                                              void * params);

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);
double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);
double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit);



void CalculateMassAndRenormalizedChemicalPotentialSimultaneously(double barionic_density,
                                                                 double * return_mass,
                                                                 double * return_renormalized_chemical_potential)
{
    multi_dim_gap_eq_param p;
    p.barionic_density = barionic_density;
    
    gsl_multiroot_function f;
    f.f = &ZeroedGapAndBarionicDensityEquations;
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

int ZeroedGapAndBarionicDensityEquations(const gsl_vector * x,
                                         void * p,
                                         gsl_vector * f)
{
    multi_dim_gap_eq_param * params = (multi_dim_gap_eq_param *)p;
    
   	const double mass = gsl_vector_get(x,0);
   	const double renormalized_chemical_potential = gsl_vector_get(x,1);
    
    double integral = FermiDiracDistributionFromDensityIntegral(mass,
                                                                renormalized_chemical_potential);
    
    double integral_2 = FermiDiracDistributionIntegralFromGapEquation(mass,
                                                                      renormalized_chemical_potential);
    
    double zeroed_gap_eq = mass
                           - parameters.bare_mass
                           - 2.0 * NUM_COLORS * NUM_FLAVORS * mass
                             * parameters.G_S * integral_2 / pow(M_PI * CONST_HBAR_C, 2.0);
    double zeroed_bar_dens_eq = params->barionic_density
                                - NUM_COLORS * NUM_FLAVORS * integral
                                  / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
    
   	gsl_vector_set (f, 0, zeroed_gap_eq);
   	gsl_vector_set (f, 1, zeroed_bar_dens_eq);
    
    return GSL_SUCCESS;
}

double FermiDiracDistributionFromDensityIntegral(double mass,
                                                 double renormalized_chemical_potential)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &FermiDiracDistributionFromDensityIntegralIntegrand;
    F.params = &p;
    
    return OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);
}

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit)
{
    gsl_integration_workspace * workspace =
    gsl_integration_workspace_alloc(parameters.fermi_dirac_integrals_max_interval_num);
    
    double integral;
    double abserr;
    gsl_integration_qag(F,
                        lower_limit,
                        upper_limit,
                        parameters.fermi_dirac_integrals_abs_error,
                        parameters.fermi_dirac_integrals_rel_error,
                        parameters.fermi_dirac_integrals_max_sub_interval,
                        parameters.fermi_dirac_integrals_integration_key,
                        workspace,
                        &integral,
                        &abserr);
    
    gsl_integration_workspace_free(workspace);
    
    return integral;
}

double FermiDiracDistributionFromDensityIntegralIntegrand(double momentum, void * params)
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

double FermiDiracDistributionIntegralFromGapEquation(double mass,
                                                     double renormalized_chemical_potential)
{
    // Maybe it's possible to write the first term as the scalar density
    // at zero temperature
    
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromGapEquationIntegrand;
    F.params = &p;
    
    return OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);
}

double FermiDiracDistributionIntegralFromGapEquationIntegrand(double momentum,
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
    return (1.0 - particle_term - antiparticle_term) * pow(momentum, 2.0) / E;
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
    return 1.0 / (1.0 + exp((energy + chemical_potential)/temperature));
}
