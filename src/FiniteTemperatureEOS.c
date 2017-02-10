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
#include "CommandlineOptions.h"
#include "FiniteTemperatureEOS.h"
#include "RootFinding.h"
#include "FermiDiracDistributions.h"

typedef struct _multi_dim_root_params{
    double barionic_density;
} multi_dim_root_params;

int MultiDimensionalRootFinderHelperFunction(const gsl_vector *x,
                                             void             *p,
                                             gsl_vector       *return_values);

double ZeroMassSpecialCaseHelperFunction(double x,
                                         void  *par);

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit);

double ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(double mass,
                                                                     double renormalized_chemical_potential);
double ThermodynamicPotentialForFiniteTemperatureFreeGasContributionIntegrand(double momentum, void * parameters);


double g(double t, double e, double c);

/*  Solution for gap equation, beta equilibrium and charge neutrality:
 *  Prototype:
 *          void SolveMultiRoots(double  barionic_density,
 *                              double *return_mass,
 *                              double *return_proton_fraction);
 *
 *  Purpose:
 *      This function will take parameters like initial guess, erros, tolerance
 *      and call the rootfinding functions in a proper way. It exists to handle
 *      mappings (see below) and to handle the special case of mass = 0 (also
 *      below). It will return the mass which solves the gap function and the
 *      proton fraction for which the beta equilibrium is respected.
 *
 *  Mappings:
 *      The variables for the root finding are assumed to cover the range
 *      (-\infty, +\infty), but that is not the case for the variables
 *      that we are trying to solve. Here both the mass $m$ and the renormalized
 *      chemical potential are such that
 *          $m \in [0, +\infty)$,
 *          $\mu_R \in [0, +\infty)$.
 *      To solve that, we use the mappings:
 *          $m = x^2$,
 *          $\mu_R = y^2$.
 *      The initial guesses must be transformed by inverting the relations
 *      above.
 *
 *  Zero mass special case:
 *      The zero mass case is important as most calculations will be
 *      performed at this particular case, which due to characteristics
 *      of the multidimensional root-finding algorithm, may be problematic to
 *      solve (it works most of the time, but sometimes calculations result in NaNs).
 *      This is due to problems in the calculation of derivatives of
 *      the function with respect to mass which arise from low variability
 *      of the function near zero mass.
 *
 *      This case, however is not the one reached at the start of calculations.
 *      In addition to that, once it is reached, all subsequent calculations are
 *      performed at approximatelly zero mass.
 *
 *      We take these characteristics into account and do the zero mass case
 *      with a special path, where we just assume mass = 0 (this effectivelly
 *      reduces the dimension of the system). This will avoid
 *      any calculation of potentially problematic derivatives. The special
 *      path is triggered by the condition
 *          mass < mass_tolerance
 *      which must be true. The tolerance should be adjusted in the setup of
 *      parameters.
 */

void SolveMultiRoots(double  barionic_density,
                     double *return_mass,
                     double *return_renorm_chem_pot)
{
    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;
    
    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.barionic_density = barionic_density;
    
    // Check for zero mass special case. As mass != 0 is the
    // case that appears first, it is implemented first.
    if (parameters.finite_temperature.mass_guess > ZERO_MASS_TOL){
        
        gsl_multiroot_function f;
        f.f = &MultiDimensionalRootFinderHelperFunction;
        f.n = dimension;
        f.params = (void *)&p;
        
        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);
        
        gsl_vector_set(initial_guess, 0, sqrt(parameters.finite_temperature.mass_guess));
        gsl_vector_set(initial_guess, 1, sqrt(parameters.finite_temperature.renor_chem_pot_guess));
        
        int status = MultidimensionalRootFinder(dimension,
                                                &f,
                                                initial_guess,
                                                parameters.finite_temperature.abs_error,
                                                parameters.finite_temperature.rel_error,
                                                parameters.finite_temperature.max_iter,
                                                return_results);

        if (status != 0){
            printf("Something is wrong with the rootfinding.\n");
            abort();
        }
        
        // Save results in return variables,
        // taking care of the mappinps
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_renorm_chem_pot = pow(gsl_vector_get(return_results, 1), 2.0);
        
        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);
        
        // Save solution as guess for next iteration
        parameters.finite_temperature.mass_guess = *return_mass;
        parameters.finite_temperature.renor_chem_pot_guess = *return_renorm_chem_pot;
        
        return;

    }
    else{ // Handle special case: Zero mass case

        gsl_function F;
        F.function = &ZeroMassSpecialCaseHelperFunction;
        F.params = &p;
        
        // Set root bounds observing the mappings
        double lower_bound = sqrt(parameters.finite_temperature.renor_chem_pot_lower_bound);
        double upper_bound = sqrt(parameters.finite_temperature.renor_chem_pot_upper_bound);
        
        // As we are left with just one variable and one equation to solve,
        // now an one-dimensional algorithm may be employed. Otherwise,
        // the dimension ought to be decreased by one an the multidimensional
        // employed again.
        double return_result;
        
        int status = UnidimensionalRootFinder(&F,
                                              lower_bound,
                                              upper_bound,
                                              parameters.finite_temperature.abs_error,
                                              parameters.finite_temperature.rel_error,
                                              parameters.finite_temperature.max_iter,
                                              &return_result);
        if (status != 0){
            printf("\nBounds do not straddle root.\n");
            abort();
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = 0.0;
        *return_renorm_chem_pot = pow(return_result, 2.0);
        return;
    }
}

double ZeroMassSpecialCaseHelperFunction(double  x,
                                         void   *par)
{
    const int dimension = 2;
    
    gsl_vector * input_values = gsl_vector_alloc(dimension);
    gsl_vector * return_values = gsl_vector_alloc(dimension);
    
    // Set mass = 0, which is our special case
    gsl_vector_set(input_values, 0, 0);
    
    // Pass value selected by the root finding routine
    gsl_vector_set(input_values, 1, x);
    
    MultiDimensionalRootFinderHelperFunction(input_values, par, return_values);
    
    double return_value = gsl_vector_get(return_values, 1);
    
    gsl_vector_free(input_values);
    gsl_vector_free(return_values);
    
    return return_value;
}

int MultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                             void               *p,
                                             gsl_vector         *return_values)
{
    multi_dim_root_params * params = (multi_dim_root_params *)p;
    double barionic_density = params->barionic_density;
    
    // Mappings:
    //      The variables for the root finding are assumed to cover the range
    //      (-\infty, +\infty), but that is not the case for the variables
    //      that we are trying to solve.
    //      To solve that, we use the mappings:
    //          $m = x^2$
    //          $\mu_R = y^2$
    //      The initial guesses must be transformed by inverting the relations
    //      above
    
   	const double mass = pow(gsl_vector_get(x, 0), 2.0);
   	const double renormalized_chemical_potential = pow(gsl_vector_get(x,1), 2.0);

    double zeroed_gap_eq = ZeroedGapEquationForFiniteTemperature(mass, renormalized_chemical_potential, NULL);
    double zeroed_bar_dens_eq = ZeroedBarionicDensityEquationForFiniteTemperature(mass,
                                                                                  renormalized_chemical_potential,
                                                                                  (void *)&barionic_density);
    
   	gsl_vector_set(return_values, 0, zeroed_gap_eq);
   	gsl_vector_set(return_values, 1, zeroed_bar_dens_eq);
    
    return GSL_SUCCESS;
}

double ZeroedGapEquationForFiniteTemperature(double mass, double renormalized_chemical_potential, void * params)
{
    double integral = FermiDiracDistributionIntegralFromGapEquation(mass,
                                                                    renormalized_chemical_potential);
    
    double constant = 2.0 * NUM_COLORS * NUM_FLAVORS * parameters.model.G_S / pow(M_PI * CONST_HBAR_C, 2.0);
    
    return mass - parameters.model.bare_mass - constant * mass * integral;
}

double ZeroedBarionicDensityEquationForFiniteTemperature(double mass,
                                                         double renormalized_chemical_potential,
                                                         void * barionic_density)
{
    double density = *((double *)barionic_density);
    double integral = FermiDiracDistributionFromDensityIntegral(mass,
                                                                renormalized_chemical_potential);
    
    double quarks_dens = NUM_COLORS * NUM_FLAVORS * integral / (pow(M_PI, 2.0));

    return 3.0 * density - quarks_dens /  pow(CONST_HBAR_C, 3.0);
}

double ThermodynamicPotentialForFiniteTemperature(double mass,
                                                  double chemical_potential,
                                                  double renormalized_chemical_potential)
{
    double first_term = ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(mass,
                                                                                      renormalized_chemical_potential);

    double second_term = pow(mass - parameters.model.bare_mass, 2.0) / (4.0 * parameters.model.G_S * CONST_HBAR_C);

    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.model.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0) //TODO: Check sign
        / (4.0 * parameters.model.G_V * CONST_HBAR_C);

    return first_term + second_term - third_term;
}

double ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(double mass,
                                                                     double renormalized_chemical_potential)
{
    therm_pot_free_gas_contrib_params p;
    p.mass = mass;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.finite_temperature.temperature;

    gsl_function F;
    F.function = &ThermodynamicPotentialForFiniteTemperatureFreeGasContributionIntegrand;
    F.params = &p;

    double integral = OnedimensionalIntegrator(&F, 0.0, parameters.model.cutoff);

    return - NUM_COLORS * NUM_FLAVORS * integral / pow(M_PI, 2.0);
}

double ThermodynamicPotentialForFiniteTemperatureFreeGasContributionIntegrand(double momentum, void * parameters)
{
    therm_pot_free_gas_contrib_params * p = (therm_pot_free_gas_contrib_params *)parameters;

    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    /* From docs:
     If x is nearly zero, then the common expression log(1 + x) will not be able to produce accurate
     results, as most (or all) of the information in x will be lost by addition.  Instead, use log1p(x) to
     perform the same computation without undue loss of accuracy. */
    double first_term = p->temperature * log1p(exp(-(energy - p->renormalized_chemical_potential)/p->temperature));
    double second_term = p->temperature * log1p(exp(-(energy + p->renormalized_chemical_potential)/p->temperature));

    return pow(momentum, 2.0) * (energy + first_term + second_term) / pow(CONST_HBAR_C, 3.0);
}

double Entropy(double mass, double temperature, double renormalized_chemical_potential)
{
	EntropyParameters params = parameters.finite_temperature.entropy;

    entropy_integrand_parameters p;
    p.mass = mass;
    p.temperature = temperature;
    p.renormalized_chemical_potential = renormalized_chemical_potential;

    gsl_function F;
    F.function = &EntropyIntegrand;
    F.params = &p;

    gsl_integration_workspace * workspace =
        gsl_integration_workspace_alloc(params.max_sub_interval);

	double integral;
	double abserr;
    gsl_integration_qag(&F,
                        params.lower_limit,
                        params.upper_limit,
                        params.abs_error,
                        params.rel_error,
                        params.max_sub_interval,
                        params.integration_key,
                        workspace,
                        &integral,
                        &abserr);

    gsl_integration_workspace_free(workspace);

    return NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * integral / pow(M_PI, 2.0);
}

// FIXME: Remove once I'm sure the other version is working ok
double EntropyIntegrandFromDerivative(double momentum, void * parameters)
{
    // This is the expression I obtained from $- \partial \omega / \partial T$

    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;

    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double fd_dist_part = FermiDiracDistributionForParticles(energy,
                                                             p->renormalized_chemical_potential,
                                                             p->temperature);

    double fd_dist_antipart = FermiDiracDistributionForAntiparticles(energy,
                                                                     p->renormalized_chemical_potential,
                                                                     p->temperature);

    double particles_term = - log1p(-fd_dist_part)
                            + (energy - p->renormalized_chemical_potential) * fd_dist_part / p->temperature;

    double antiparticles_term = - log1p(-fd_dist_antipart)
                                + (energy + p->renormalized_chemical_potential) * fd_dist_antipart / p->temperature;

    return pow(momentum, 2.0) * (particles_term + antiparticles_term);
}

double EntropyIntegrandArt(double momentum, void * parameters)
{
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;
    
    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double fd_dist_part = FermiDiracDistributionForParticles(energy,
                                                             p->renormalized_chemical_potential,
                                                             p->temperature);

    double fd_dist_antipart = FermiDiracDistributionForAntiparticles(energy,
                                                                     p->renormalized_chemical_potential,
                                                                     p->temperature);

    double particles_term = fd_dist_part * log(fd_dist_part)
                            + (1.0 - fd_dist_part) * log1p(-fd_dist_part);

    double antiparticles_term = fd_dist_antipart * log(fd_dist_antipart)
                                + (1.0 - fd_dist_antipart) * log1p(-fd_dist_antipart);

    return -pow(momentum, 2.0) * (particles_term + antiparticles_term);
}


double g(double t, double e, double c)
{
    double a = - (e - c)/t;
    return log1p(exp(a)) - a * exp(a) / (1.0 + exp(a));
}

double EntropyIntegrand(double momentum, void * parameters)
{
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;
    
    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double first_term = g(p->temperature, energy, p->renormalized_chemical_potential);
    double second_term = g(p->temperature, energy, -p->renormalized_chemical_potential);

    return pow(momentum, 2.0) * (first_term + second_term);
}

double PressureForFiniteTemperature(double regularized_thermodynamic_potential)
{
    return - regularized_thermodynamic_potential;
}

double EnergyForFiniteTemperature(double regularized_thermodynamic_potential,
                                  double chemical_potential,
                                  double barionic_density,
                                  double temperature,
                                  double entropy)
{
    return regularized_thermodynamic_potential
            + temperature * entropy
            + chemical_potential * NUM_COLORS * barionic_density;
}


