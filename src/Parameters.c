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

Parameters parameters;

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.
    
    parameters.points_number = 1000;
    parameters.minimum_density = 0.01;
    parameters.maximum_density = 0.2625;
    
    parameters.gap_minimum_mass = 0.0;
    parameters.gap_maximum_mass = 1000.0;
    parameters.gap_points_number = 2000;
    
    parameters.gap_eq_solver_lower_bound = 1.0E-3;  // Low, but not zero. In zero f(M) = 0 if bare_mass is zero,
                                                    // and that may be a problem for root finding;
    
    parameters.gap_eq_solver_upper_bound = 1000.0;  // MeV (about three times the nucleon mass)
    parameters.gap_eq_solver_abs_error = 0.005;
    parameters.gap_eq_solver_rel_error = 5.0E-5;
    parameters.gap_eq_solver_max_iterations = 2000;
    
    parameters.vac_mass_det_lower_bound = 1.0E-3;   // Low, but not zero. Same reason as above.
    parameters.vac_mass_det_upper_bound = 1000.0;   // Near the value of the nucleon mass.
    parameters.vac_mass_det_abs_error = 0.05;
    parameters.vac_mass_det_rel_error = 5.0E-4;
    parameters.vac_mass_det_max_iterations = 1000;
	

    // To make the use of different parameterizations easier, define a name in the enum
    // and assign values for the variables in a case of the switch below.
    
    typedef enum{
        Buballa_1,  // Set 1 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408
        Buballa_2,  // Set 2 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408
        Buballa_3,  // Set 3 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408
        D_1         // Set found in Debora's code
    }Parameters_set;
    
    Parameters_set par = Buballa_1; // Default value
    
    if (NULL != options.parameterization){
    	// strcasecmp returns 0 if the strings are equal (in a case insensitive manner),
    	// but 0 is false, so it's necessary to add the ! (not) operator.
    	if (!strcasecmp("Buballa_1", options.parameterization)){
    		par = Buballa_1;
    	}
    	else if (!strcasecmp("Buballa_2", options.parameterization)){
    		par = Buballa_2;
    	}
        else if (!strcasecmp("Buballa_3", options.parameterization)){
            par = Buballa_3;
        }
        else if (!strcasecmp("D_1", options.parameterization)){
            par = D_1;
        }
    	else{
    		printf("You must choose from a defined parameterization.\n");
    		exit(EXIT_FAILURE);
    	}
    }
    
    switch (par) {
        case Buballa_1:
            parameters.parameters_set_identifier = "Buballa_1";
            parameters.G_S = pow(CONST_HBAR_C, 2.0) * 2.14 / pow(650, 2.0); // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
            parameters.cutoff = 650;  // MeV
            parameters.bare_mass = 0.0; // MeV
            break;
        case Buballa_2:
            parameters.parameters_set_identifier = "Buballa_2";
            parameters.G_S = pow(CONST_HBAR_C, 2.0) * 2.45 / pow(600, 2.0); // Same as Buballa_1
            parameters.cutoff = 600;
            parameters.bare_mass = 0.0;
            break;
        case Buballa_3:
            parameters.parameters_set_identifier = "Buballa_3";
            parameters.G_S = pow(CONST_HBAR_C, 2.0) * 2.84 / pow(570, 2.0); // Same as Buballa_1
            parameters.cutoff = 570;
            parameters.bare_mass = 0.0;
            break;
        case D_1:
            parameters.parameters_set_identifier = "D_1";
            parameters.G_S = pow(CONST_HBAR_C, 2.0) * 2.44 / pow(570, 2.0); // Same as Buballa_1
            parameters.cutoff = 587.9;
            parameters.bare_mass = 5.6;
            break;
        default:
            printf("You must define the parameters of the parameterization.\n");
            exit(EXIT_FAILURE);
            break;
    }
}

