//
//  ZeroTemperatureEOS.c
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
#include "ZeroTemperatureEOS.h"
#include "RootFinding.h"

double F0(double mass, double momentum);
double F_E(double mass, double momentum);


double GapEquationSolver(double fermi_momentum)
{
    // Prepare input for ZeroedGapEquation
    gap_equation_input input;
    input.fermi_momentum = fermi_momentum;
    
    // Prepare function to be passed to the root finding algorithm
    gsl_function F;
    F.function = &ZeroedGapEquation;
    F.params = &input;
    
    double root;
    int status = UnidimensionalRootFinder(&F,
                                          parameters.zero_temperature.lower_bound,
                                          parameters.zero_temperature.upper_bound,
                                          parameters.zero_temperature.abs_error,
                                          parameters.zero_temperature.rel_error,
                                          parameters.zero_temperature.max_iterations,
                                          &root);
    if (status == -1){
        return 0;
    }
    
    return root;
}

double ZeroedGapEquation(double mass, void * input)
{
	gap_equation_input * param = (gap_equation_input *)input;
	
	double scalar_density = ScalarDensity(mass, param->fermi_momentum);

	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.model.G_S * scalar_density;

	return mass + gap_1st_term - parameters.model.bare_mass;
}

double VacuumMassDetermination()
{
    // Prepare function to be passed to the root finding algorithm.
    // No parameters are needed.
    gsl_function F;
    F.function = &VacuumMassEquation;
    
    double root;
    int status = UnidimensionalRootFinder(&F,
                                          parameters.vacuum_mass_determination.lower_bound,
                                          parameters.vacuum_mass_determination.upper_bound,
                                          parameters.vacuum_mass_determination.abs_error,
                                          parameters.vacuum_mass_determination.rel_error,
                                          parameters.vacuum_mass_determination.max_iterations,
                                          &root);
    
    if (status == -1){
        return 0;
    }
    
    return root;
}

double VacuumMassEquation(double mass, void * input)
{
    double F_diff = F0(mass, parameters.model.cutoff) - F0(mass, 0.0);
    double term = 2.0 * NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -2.0)
  				  * parameters.model.G_S * mass * F_diff
  				  / pow(M_PI, 2.0);
    
    return mass - parameters.model.bare_mass - term;
}

double ZeroedRenormalizedChemicalPotentialEquation(double renor_chem_pot,
												   void * input)
{
  	renorm_chem_pot_equation_input * param = (renorm_chem_pot_equation_input *) input;

  	double c = 2.0 * parameters.model.G_V * NUM_COLORS * NUM_FLAVORS
  				   / (3.0 * pow(M_PI * CONST_HBAR_C, 2.0));
    
    // The 'if' statement does the work of a step function
    double arg = 0;
    if (pow(renor_chem_pot, 2.0) - pow(param->mass, 2.0) >= 0){
        arg = pow(renor_chem_pot, 2.0) - pow(param->mass, 2.0);
	}

	return renor_chem_pot - param->chemical_potential + c * pow(arg, 3.0 / 2.0);
}

double ScalarDensity(double mass, double fermi_momentum)
{
    if (mass == 0){
        return 0;
    }

	return NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0) * (mass / pow(M_PI, 2.0))
           * (F0(mass, fermi_momentum) - F0(mass, parameters.model.cutoff));
}

double F0(double mass, double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
	
	return (1.0 / 2.0) * (momentum * E - pow(mass, 2.0) * log((momentum + E) / mass));
}

double ThermodynamicPotential(double mass,
                              double fermi_momentum,
                              double chemical_potential,
                              double renormalized_chemical_potential)

{
    double F_diff = F_E(mass, parameters.model.cutoff) - F_E(mass, fermi_momentum);
    
    double first_term = - NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0)
                        * (F_diff + renormalized_chemical_potential * pow(fermi_momentum, 3.0) / 3.0)
                        / pow(M_PI, 2.0);
    double second_term = pow(mass - parameters.model.bare_mass, 2.0)
  						 / (4.0 * parameters.model.G_S * CONST_HBAR_C);
    
    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.model.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0)
                     / (4.0 * parameters.model.G_V * CONST_HBAR_C);
    
    return first_term + second_term - third_term;
}

double F_E(double mass, double momentum)
{
    if (mass == 0)
        return pow(momentum, 4.0) / 4.0;

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
