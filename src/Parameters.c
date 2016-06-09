//
//  Parameters.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

// Chosen parameter set (globally accessible)
Parameters parameters;
static Parameters parameters_sets_list[256] = {0};

void NewCopyOfDefaultParametersSet();

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.
    int index = 0;

  	// Start list with default for easier listing of available parameterizations
  	parameters_sets_list[0] = NewCopyOfDefaultParametersSet();
  	index++;

  	// START DECLARATION OF PARAMETERS SETS:
	//
  	// 	For declaration of parameters sets:
	// 		- copy from the default one,
	// 		- overwrite values as desired,
	// 		- append to the list of parameterizations

	Parameters p;

  	p = NewCopyOfDefaultParametersSet();
    p.parameters_set_identifier = "Buballa_1";
  	p.parameters_set_origin = "Set 1 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 650.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.14 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = 0.0;

  	// Append to list of parameterizations
  	parameters_sets_list[index] = p;
  	index ++;

  	////

  	p = NewCopyOfDefaultParametersSet();
    p.parameters_set_identifier = "Buballa_2";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 600.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.45 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = 0.0;

  	// Append to list of parameterizations
  	parameters_sets_list[index] = p;
  	index ++;

  	////

  	p = NewCopyOfDefaultParametersSet();
    p.parameters_set_identifier = "Buballa_3";
  	p.parameters_set_origin = "Set 3 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 570.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.84 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = 0.0;

  	// Append to list of parameterizations
  	parameters_sets_list[index] = p;
  	index ++;

	////

  	p = NewCopyOfDefaultParametersSet();
    p.parameters_set_identifier = "BuballaR_2";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports 407 (2005) 205-376 with G_V = 0";
    p.bare_mass = 5.6;                                             // MeV
    p.cutoff = 587.9;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = 0.0;

  	// Append to list of parameterizations
  	parameters_sets_list[index] = p;
  	index ++;

  	////

  	p = NewCopyOfDefaultParametersSet();
    p.parameters_set_identifier = "BuballaR_2";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports 407 (2005) 205-376 with G_V = G_S";
    p.bare_mass = 5.6;                                             // MeV
    p.cutoff = 587.9;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = p.G_S;

  	// Append to list of parameterizations
  	parameters_sets_list[index] = p;
  	index ++;

  	// END DECLARATION OF PARAMETERS SETS

  	// If asked to, print parameters sets and exit
  	if (options.list_available_parameterizations){
		for (int i = 0; i < index; i++){
		  	printf("Parameters set %s\n"
				   "\tOrigin: %s\n",
				   parameters_sets_list[i].parameters_set_identifier,
				   parameters_sets_list[i].parameters_set_origin);
		}
		exit(EXIT_SUCCESS);
	}
}

void SetParametersSet(char * parameters_set_identifier)
{
	  for (int i = 0; i < index; i++)
		if (!strcasecmp(parameters_sets_list[i].parameters_set_identifier,
				   		parameters_set_identifier)){
			parameters = parameters_sets_list[i];
		}
  		else{
		  	printf("Parameters set %s unrecognized.\n"
				   "Use -l to list available parameterizations.\n",
				   options.parameterization);
		  	exit(EXIT_FAILURE);
		}
}

void NewCopyOfDefaultParametersSet()
{
  	Parameters p;

  	p.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports 407 (2005) 205-376";
    p.parameters_set_identifier = "BuballaR_2";
    p.bare_mass = 5.6;
    p.cutoff = 587.9;

  	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / parameters.cutoff, 2.0);
    p.G_V = 0.0;

  	p.points_number = 1000;
    p.minimum_density = 0.01;
    p.maximum_density = 0.2625;

    p.gap_minimum_mass = 0.0;
    p.gap_maximum_mass = 1000.0;
    p.gap_points_number = 2000;

    // Lower bound: Low, but not zero. In zero f(M) = 0 if bare_mass is zero,
    // and that may be a problem for root finding;
    p.gap_eq_solver_lower_bound = 1.0E-3;
    p.gap_eq_solver_upper_bound = 1000.0;  // MeV (about three times the nucleon mass)
    p.gap_eq_solver_abs_error = 0.005;
    p.gap_eq_solver_rel_error = 5.0E-5;
    p.gap_eq_solver_max_iterations = 2000;

    p.vac_mass_det_lower_bound = 1.0E-3;   // Low, but not zero. Same reason as above.
    p.vac_mass_det_upper_bound = 1000.0;   // Near the value of the nucleon mass.
    p.vac_mass_det_abs_error = 0.005;
    p.vac_mass_det_rel_error = 5.0E-5;
    p.vac_mass_det_max_iterations = 2000;

  	p.renormalized_chemical_potential_lower_bound = 1.0E-3;
    p.renormalized_chemical_potential_upper_bound = 1000.0;
    p.renormalized_chemical_potential_abs_error = 0.005;
    p.renormalized_chemical_potential_rel_error = 5.0E-5;
    p.renormalized_chemical_potential_max_iter = 2000;

  	return p;
}

