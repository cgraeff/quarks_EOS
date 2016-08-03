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

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit);

double ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(double mass,
                                                                     double renormalized_chemical_potential);
double ThermodynamicPotentialForFiniteTemperatureFreeGasContributionIntegrand(double momentum, void * parameters);


double g(double t, double e, double c);


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
        
        if (status == GSL_EBADFUNC){
            printf("TwodimensionalRootFinder: Error: Infinity or division by zero.\n");
            abort();
        }
        else if (status == GSL_ENOPROG){
            printf("TwodimensionalRootFinder: Error: Solver is stuck. Try a different initial guess.\n");
            abort();
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
    
    double zeroed_gap_eq = ZeroedGapEquationForFiniteTemperature(mass, renormalized_chemical_potential);
    double zeroed_bar_dens_eq = ZeroedBarionicDensityEquationForFiniteDensity(mass,
                                                                              renormalized_chemical_potential,
                                                                              params->barionic_density);

   	gsl_vector_set(f, 0, zeroed_gap_eq);
   	gsl_vector_set(f, 1, zeroed_bar_dens_eq);
    
    return GSL_SUCCESS;
}

double ZeroedGapEquationForFiniteTemperature(double mass, double renormalized_chemical_potential)
{
    double integral = FermiDiracDistributionIntegralFromGapEquation(mass,
                                                                    renormalized_chemical_potential);
    
    double constant = 2.0 * NUM_COLORS * NUM_FLAVORS * parameters.G_S / pow(M_PI * CONST_HBAR_C, 2.0);
    
    return mass - parameters.bare_mass - constant * mass * integral;
}

double ZeroedBarionicDensityEquationForFiniteDensity(double mass,
                                                     double renormalized_chemical_potential,
                                                     double barionic_density)
{
    double integral = FermiDiracDistributionFromDensityIntegral(mass,
                                                                renormalized_chemical_potential);
    
    double quarks_dens = NUM_COLORS * NUM_FLAVORS * integral / (pow(M_PI, 2.0));
    
    return 3.0 * barionic_density - quarks_dens /  pow(CONST_HBAR_C, 3.0);
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

  	double integral = OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);

  	return integral;
}

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit)
{
    gsl_integration_workspace * workspace =
    gsl_integration_workspace_alloc(parameters.fermi_dirac_integrals_max_interval_num);
    
    double integral = 0;
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
    
	double integral = OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);

    return integral;
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

double ThermodynamicPotentialForFiniteTemperature(double mass,
                                                  double chemical_potential,
                                                  double renormalized_chemical_potential)
{
    double first_term = ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(mass,
                                                                                      renormalized_chemical_potential);
    
    double second_term = pow(mass - parameters.bare_mass, 2.0)
    / (4.0 * parameters.G_S * CONST_HBAR_C);
    
    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0)
        / (4.0 * parameters.G_V * CONST_HBAR_C);
    
    // Trying + sign on third term like in T=0 case
    return first_term + second_term + third_term;
}

double ThermodynamicPotentialForFiniteTemperatureFreeGasContribution(double mass,
                                                                     double renormalized_chemical_potential)
{
    therm_pot_free_gas_contrib_params p;
    p.mass = mass;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &ThermodynamicPotentialForFiniteTemperatureFreeGasContributionIntegrand;
    F.params = &p;
    
    double integral = OnedimensionalIntegrator(&F, 0.0, parameters.cutoff);
    
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
    int interval_num = 1000;
    
    entropy_integrand_parameters p;
    p.mass = mass;
    p.temperature = temperature;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    
    gsl_function F;
    F.function = &EntropyIntegrand;
    F.params = &p;
    
    gsl_integration_workspace * workspace =
        gsl_integration_workspace_alloc(interval_num);
    
    double integral = 0;
    double abserr;
    double lower_limit = 0.0;
    double upper_limit = parameters.cutoff;
    double abs_error = 1.0E-3;
    double rel_error = 1.0E-3;
    int max_sub_interval = interval_num;
    int integration_key = GSL_INTEG_GAUSS61;
    
    gsl_integration_qag(&F,
                        lower_limit,
                        upper_limit,
                        abs_error,
                        rel_error,
                        max_sub_interval,
                        integration_key,
                        workspace,
                        &integral,
                        &abserr);
    
    gsl_integration_workspace_free(workspace);
    
    return NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * integral / pow(M_PI, 2.0);
}

// FIXME: Remove once I'm sure the other version is working ok
/*
double EntropyIntegrandFromDerivative(double momentum, void * parameters)
{
    // This is the expression I obtained from - \partial \omega / \partial T
    
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
*/

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

void MapFunction(double (*function)(double, double),
                 double min_x,
                 double max_x,
                 int x_num_pts,
                 double min_y,
                 double max_y,
                 int y_num_pts,
                 double tolerance,
                 gsl_vector * output_x,
                 gsl_vector * output_y,
                 int *num_points,
                 bool show_progress)
{
    double x_step = (max_x - min_x) / (double)(x_num_pts - 1);
    double y_step = (max_y - min_y) / (double)(y_num_pts - 2);
    
    int progress_counter = 0;
    int pts_counter = 0;
    double x = min_x;
    double y = min_y;
    
    for (int i = 0; i < x_num_pts; i++) {
        for (int j = 0; j < y_num_pts; j++){
            
            double f = function(x, y);
            
            if (fabs(f) <= tolerance){
                gsl_vector_set(output_x, pts_counter, x);
                gsl_vector_set(output_y, pts_counter, y);
                pts_counter++;
            }
            
            if (show_progress)
                if ((100 * progress_counter) % (x_num_pts * y_num_pts) == 0){
                    printf("\rMap %d %%", (100 * progress_counter) / (x_num_pts * y_num_pts));
                    fflush(stdout);
                }
            
            progress_counter++;
            y += y_step;
        }

        x += x_step;
    }
    
    *num_points = pts_counter;
    return;
}

// We assume that both maps are made using the same grid
// No steps are taken to verify the possibility of two intersections.
// In this routine, we will sweep the common region of the x axis,
// determining the point where the distance in the y axis between
// points in both maps is minimal.
// We take some steps to reduce multiple points with the same value of x to one point
// with y given by the mean.
// Once the point of minimal y distance is found, we use the previous one and the next one to
// calculate two linear fits, one for each maps, and the determine the intersection of the fits
//
// Separar essas funções em um arquivo à parte
void IntersectionPointOfTwoMaps(gsl_vector * map1_x,
                                gsl_vector * map1_y,
                                int map1_num_points,
                                gsl_vector * map2_x,
                                gsl_vector * map2_y,
                                int map2_num_points)
{
    //
    // Determine common range in x axis and rebuild
	// the index
	//
    
    double common_interval_start = fmax(gsl_vector_get(map1_x, 1),
                                        gsl_vector_get(map2_x, 1));
    double common_interval_ends = fmin(gsl_vector_get(map1_x, map1_num_points),
                                       gsl_vector_get(map2_x, map2_num_points));
    
    gsl_vector * common_map1_x = gsl_vector_alloc(map1_num_points);
    gsl_vector * common_map1_y = gsl_vector_alloc(map1_num_points);
    
    int common_pts_map1 = 0;
    for (int i = 0; i < map1_num_points; i++){
        double x = gsl_vector_get(map1_x, i);
        if (x >= common_interval_start){
            gsl_vector_set(common_map1_x, common_pts_map1, x);
            gsl_vector_set(common_map1_y, common_pts_map1, gsl_vector_get(map1_y, i));
            common_pts_map1++;
        }
        
        if (x > common_interval_ends)
            break;
    }

  // separar em função para não ficar repetido aí em baixo, é praticamente a mesma coisa
    gsl_vector * common_map2_x = gsl_vector_alloc(map2_num_points);
    gsl_vector * common_map2_y = gsl_vector_alloc(map2_num_points);
    
    int common_pts_map2 = 0;
    for (int i = 0; i < map2_num_points; i++){
        double x = gsl_vector_get(map2_x, i);
        if (x >= common_interval_start){
            gsl_vector_set(common_map2_x, common_pts_map2, x);
            gsl_vector_set(common_map2_y, common_pts_map2, gsl_vector_get(map2_y, i));
            common_pts_map2++;
        }
        
        if (x > common_interval_ends)
            break;
    }
    
    //
    // Reduce multiple points
    //

    int red_map1_count = 0;
    gsl_vector * red_map1_x = gsl_vector_alloc(common_pts_map1);
    gsl_vector * red_map1_y = gsl_vector_alloc(common_pts_map1);
    
    // Loop looking for repeated values of x. The multiple points
	// will be reduced to one by taking the mean of the various values for y.
    // Use first point as a reference.

    double x = gsl_vector_get(common_map1_x, 0);
    double y_sum = gsl_vector_get(common_map1_y, 0);
    int multiplicity = 1;

    for (int i = 1; i < common_pts_map1; i++){

        double x_i = gsl_vector_get(common_map1_x, i);

	  	// if the value of x is repeated, add the
		// value of y and increase the multiplicity
        if (x_i == x){
            y_sum += gsl_vector_get(common_map1_y, i);
            multiplicity++;
        }
        else{
            // x has changed, save it and corresponding mean y,
			// increase the number of reduced points
            gsl_vector_set(red_map1_x, red_map1_count, x);
            gsl_vector_set(red_map1_y, red_map1_count, y_sum / (double)multiplicity);
            red_map1_count++;
            
            // prepare next iteration
            x = x_i;
            y_sum = gsl_vector_get(common_map1_y, i);
            multiplicity = 1;
        }
    }
    // save last point
    gsl_vector_set(red_map1_x, red_map1_count, x);
    gsl_vector_set(red_map1_y, red_map1_count, y_sum / (double)multiplicity);
    red_map1_count++;
    
  	// separa o acima em uma função e repetir para o segundo mapa
    int red_map2_count = 0;
    gsl_vector * red_map2_x = gsl_vector_alloc(common_pts_map2);
    gsl_vector * red_map2_y = gsl_vector_alloc(common_pts_map2);

  	if (red_map1_count != red_map2_count){
		printf("I expect the number of points to be equal on both maps"
			   "as it should depend only in choice for the grid");
	    abort();
	}

  	//
	// Search for the value of x in which we have the least diference in
	// the values for y
	//

    double y_diff = fabs(gsl_vector_get(red_map1_y, 0) - gsl_vector_get(red_map2_y, 0));
    int least_y_diff_index = 0;

	for (int i = 1; i < red_map1_count; i++){

        double d = fabs(gsl_vector_get(red_map1_y, i) - gsl_vector_get(red_map2_y, i));

        if (d < y_diff){

          y_diff = d;
          least_y_diff_index = i;
        }
	}
    
    // usar os valores dos pontos associados aos índices anterior e posterior para
    // traçar duas retas e determinar a interseção entre elas.
    
    // descartar os vetores criados aqui

    // Retornar
    // os valores de x e y para os quais ocorre a interseção.
}
