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

#include "Parameters.h"
#include "CommandlineOptions.h"

Parameters parameters;

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.
    
    parameters.points_number = 1000;
    parameters.minimum_density = 0.01;
    parameters.maximum_density = 0.35;
    parameters.proton_fraction = 0.5;
    
    parameters.gap_minimum_mass = 0.0;
    parameters.gap_maximum_mass = 1000.0;
    parameters.gap_points_number = 2000;
    
    parameters.gap_eq_solver_lower_bound = 1.0E-3;   // Low, but not zero. In zero f(M) = 0 if bare_mass is zero, and that may be a problem for root finding;
    parameters.gap_eq_solver_upper_bound = 3000.0;	// MeV (about three times the nucleon mass)
    parameters.gap_eq_solver_abs_error = 0.05;
    parameters.gap_eq_solver_rel_error = 5.0E-4;
    parameters.gap_eq_solver_max_iterations = 1000;
    
    parameters.vac_mass_det_lower_bound = 1.0E-3;
    parameters.vac_mass_det_upper_bound = 3000.0;
    parameters.vac_mass_det_abs_error = 0.05;
    parameters.vac_mass_det_rel_error = 5.0E-4;
    parameters.vac_mass_det_max_iterations = 1000;
	

    // To make the use of different parameterizations easier, define a name in the enum
    // and assign values for the variables in a case of the switch below.
    
    typedef enum{
        Name_1,
        Name_2
    }Parameters_set;
    
    Parameters_set par = Name_1;
    
    if (NULL != options.parameterization){
    	// strcasecmp returns 0 if the strings are equal (in a case insensitive manner),
    	// but 0 is false, so it's necessary to add the ! (not) operator.
    	if (!strcasecmp("name_1", options.parameterization)){
    		par = Name_1;
    	}
    	else if (!strcasecmp("name_2", options.parameterization)){
    		par = Name_2;
    	}
    	else{
    		printf("You must choose from a defined parameterization.\n");
    		exit(EXIT_FAILURE);
    	}
    }
    
    switch (par) {
        case Name_1:
            parameters.G_S = 4.855;
            parameters.cutoff = 388.189;
            parameters.bare_mass = 0.0;
            break;
        case Name_2:
            parameters.G_S = 4.855;
            parameters.cutoff = 388.189;
            parameters.bare_mass = 0.0;
            break;
        default:
            printf("You must define the parameters of the parameterization.\n");
            exit(EXIT_FAILURE);
            break;
    }
}

