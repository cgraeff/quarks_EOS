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
#include <gsl/gsl_integration.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

// Chosen parameter set (globally accessible)
Parameters parameters;
static Parameters parameters_sets_list[256] = {0};
static int parameters_sets_list_count = 0;

Parameters NewCopyOfParametersSetFromTemplate();
void AppendParametersSetToList(Parameters a_set);

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.

  	// START DECLARATION OF PARAMETERS SETS:
	//
  	// 	For declaration of parameters sets:
	// 		- copy from the default one,
	// 		- overwrite values as desired,
	// 		- append to the list of parameterizations
    //
    // The first declared parameters set will be treated as the standard
    // case. To run others, use the "-p IDENTIFIER" commandline option, without
    // quotation marks) where IDENTIFIER is the content of the variable
    // parameters_set_identifier in the definitions bellow.

	Parameters p;

  	p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_1";
  	p.parameters_set_origin = "Set 1 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 650.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.14 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = 0.0;

  	AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_2";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 600.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.45 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = 0.0;

  	AppendParametersSetToList(p);

  	////

  	p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_3";
  	p.parameters_set_origin = "Set 3 from M. Buballa, Nucl. Phys. A 611 (1996) 393-408";
    p.bare_mass = 0.0;                                             // MeV
    p.cutoff = 570.0;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.84 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = 0.0;

  	AppendParametersSetToList(p);

	////

  	p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports 407 (2005) 205-376 with G_V = 0";
    p.bare_mass = 5.6;                                             // MeV
    p.cutoff = 587.9;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = 0.0;

  	AppendParametersSetToList(p);

  	////

  	p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV";
  	p.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports 407 (2005) 205-376 with G_V = G_S";
    p.bare_mass = 5.6;                                             // MeV
    p.cutoff = 587.9;                                              // MeV

	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = p.G_S;

  	AppendParametersSetToList(p);

  	// END DECLARATION OF PARAMETERS SETS
    
    // Verify that the set identifiers are unique
    for (int i = 0; i < parameters_sets_list_count; i++){
        for (int j = 0; j < parameters_sets_list_count; j++){
            
            if (i == j)
                continue;
            
            if(!strcasecmp(parameters_sets_list[i].parameters_set_identifier,
                           parameters_sets_list[j].parameters_set_identifier)){
                printf("Two parameters sets share the \"%s\" identifier, but it should be unique.\n",
                       parameters_sets_list[i].parameters_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

  	// If asked to, print parameters sets and exit
  	if (options.list_available_parameterizations){
		for (int i = 0; i < parameters_sets_list_count; i++){
		  	printf("Parameters set %s\n"
				   "\tOrigin: %s\n",
				   parameters_sets_list[i].parameters_set_identifier,
				   parameters_sets_list[i].parameters_set_origin);
		}
		exit(EXIT_SUCCESS);
	}
    
    // If there isn't at least one set, exit
    if (parameters_sets_list_count == 0){
        printf("There are not parameters set declared. Declare at least one set.\n");
        exit(EXIT_FAILURE);
    }
}

void SetParametersSet(char * parameters_set_identifier)
{
    // If the identifier is null, return default case
    if (parameters_set_identifier == NULL){
        parameters = parameters_sets_list[0];
        return;
    }
        
    for (int i = 0; i < parameters_sets_list_count; i++){
          
        Parameters p = parameters_sets_list[i];
                          
        if (!strcasecmp(p.parameters_set_identifier, parameters_set_identifier)){
			parameters = p;
            return;
		}
    }
    
    printf("Parameters set %s unrecognized.\n"
           "Use -l to list available parameterizations.\n",
           parameters_set_identifier);
    exit(EXIT_FAILURE);
}

Parameters NewCopyOfParametersSetFromTemplate()
{
  	Parameters p;

  	p.parameters_set_origin = "A standard parameters set. With theory parameters from Set 2 of M. Buballa, Physics Reports 407 (2005) 205-376";
    p.parameters_set_identifier = "Template";
    p.bare_mass = 5.6;
    p.cutoff = 587.9;

  	// The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.G_S = 2.44 * pow(CONST_HBAR_C / p.cutoff, 2.0);
    p.G_V = 0.0;

  	p.points_number = 1000;
    p.minimum_density = 0.01;
    p.maximum_density = 1.78;

    p.gap_minimum_mass = 0.0;
    p.gap_maximum_mass = 1000.0;
    p.gap_points_number = 2000;

    // Lower bound: Low, but not zero. In zero f(M) = 0 if bare_mass is zero,
    // and that may be a problem for root finding;
    p.gap_eq_solver_lower_bound = 1.0E-3;
    p.gap_eq_solver_upper_bound = 1000.0;  // MeV (about three times the nucleon mass)
    p.gap_eq_solver_abs_error = 1.0E-5;
    p.gap_eq_solver_rel_error = 5.0E-5;
    p.gap_eq_solver_max_iterations = 2000;

    p.vac_mass_det_lower_bound = 1.0E-3;   // Low, but not zero. Same reason as above.
    p.vac_mass_det_upper_bound = 1000.0;   // Near the value of the nucleon mass.
    p.vac_mass_det_abs_error = 1.0E-5;
    p.vac_mass_det_rel_error = 1.0E-5;
    p.vac_mass_det_max_iterations = 2000;

  	p.renormalized_chemical_potential_lower_bound = 1.0E-3;
    p.renormalized_chemical_potential_upper_bound = 1000.0;
    p.renormalized_chemical_potential_abs_error = 1.0E-5;
    p.renormalized_chemical_potential_rel_error = 1.0E-5;
    p.renormalized_chemical_potential_max_iter = 2000;
    
    p.temperature = 0.0; // (MeV)
    
    // In the following, the guesses are used for the first iteration.
    // Any subsequent iteration uses as guesses the values of mass
    // and renormalized chemical potential from the previous iteration.
    // This shall work well for stepping of a control variable.
    p.mass_and_renor_chem_pot_solution_max_iter = 4000;
    p.mass_and_renor_chem_pot_solution_mass_guess = 300.0; // (MeV)
    p.mass_and_renor_chem_pot_solution_renor_chem_pot_guess = 400.0; //(MeV)
    p.mass_and_renorm_chem_pot_solution_abs_error = 1.0E-5;
    p.mass_and_renorm_chem_pot_solution_rel_error = 1.0E-5;
    
    p.fermi_dirac_integrals_max_interval_num = 8000;
    p.fermi_dirac_integrals_integration_key = GSL_INTEG_GAUSS61;
    p.fermi_dirac_integrals_max_sub_interval = 8000;
    p.fermi_dirac_integrals_abs_error = 1.0E-10;
    p.fermi_dirac_integrals_rel_error = 1.0E-10;

  	return p;
}

void AppendParametersSetToList(Parameters a_set)
{
    // Append to list of parameterizations
    parameters_sets_list[parameters_sets_list_count] = a_set;
    parameters_sets_list_count++;
}

void PrintParametersToFile(FILE * file)
{
    fprintf(file, "\n--- PARAMETERS ---\n");
    fprintf(file, "parameters_set_origin = %s\n", parameters.parameters_set_origin);
    fprintf(file, "parameters_set_identifier = %s\n\n", parameters.parameters_set_identifier);
    
    fprintf(file, "G_S = %f\n", parameters.G_S);
    fprintf(file, "G_V = %f\n", parameters.G_V);
    fprintf(file, "cutoff = %f\n", parameters.cutoff);
    fprintf(file, "bare_mass = %f\n", parameters.bare_mass);
    fprintf(file, "temperature = %f\n\n", parameters.temperature);
    
    fprintf(file, "points_number = %d\n", parameters.points_number);
    fprintf(file, "minimum_density = %f\n", parameters.minimum_density);
    fprintf(file, "maximum_density = %f\n\n", parameters.maximum_density);
    
    fprintf(file, "gap_minimum_mass = %f\n", parameters.gap_minimum_mass);
    fprintf(file, "gap_maximum_mass = %f\n", parameters.gap_maximum_mass);
    fprintf(file, "gap_points_number = %f\n\n", parameters.gap_points_number);
    
    fprintf(file, "gap_eq_solver_max_iterations = %d\n", parameters.gap_eq_solver_max_iterations);
    fprintf(file, "gap_eq_solver_lower_bound = %f\n", parameters.gap_eq_solver_lower_bound);
    fprintf(file, "gap_eq_solver_upper_bound = %f\n", parameters.gap_eq_solver_upper_bound);
    fprintf(file, "gap_eq_solver_abs_error = %f\n", parameters.gap_eq_solver_abs_error);
    fprintf(file, "gap_eq_solver_rel_error = %f\n\n", parameters.gap_eq_solver_rel_error);
    
    fprintf(file, "vac_mass_det_max_iterations = %d\n", parameters.vac_mass_det_max_iterations);
    fprintf(file, "vac_mass_det_lower_bound = %f\n", parameters.vac_mass_det_lower_bound);
    fprintf(file, "vac_mass_det_upper_bound = %f\n", parameters.vac_mass_det_upper_bound);
    fprintf(file, "vac_mass_det_abs_error = %f\n", parameters.vac_mass_det_abs_error);
    fprintf(file, "vac_mass_det_rel_error = %f\n\n", parameters.vac_mass_det_rel_error);
    
    fprintf(file, "renormalized_chemical_potential_max_iter = %d\n", parameters.renormalized_chemical_potential_max_iter);
    fprintf(file, "renormalized_chemical_potential_lower_bound = %f\n", parameters.renormalized_chemical_potential_lower_bound);
    fprintf(file, "renormalized_chemical_potential_upper_bound = %f\n", parameters.renormalized_chemical_potential_upper_bound);
    fprintf(file, "renormalized_chemical_potential_abs_error = %f\n", parameters.renormalized_chemical_potential_abs_error);
    fprintf(file, "renormalized_chemical_potential_rel_error = %f\n\n", parameters.renormalized_chemical_potential_rel_error);
    
    fprintf(file, "mass_and_renor_chem_pot_solution_max_iter = %d\n", parameters.mass_and_renor_chem_pot_solution_max_iter);
    fprintf(file, "mass_and_renor_chem_pot_solution_mass_guess = %f\n", parameters.mass_and_renor_chem_pot_solution_mass_guess);
    fprintf(file, "mass_and_renor_chem_pot_solution_renor_chem_pot_guess = %f\n", parameters.mass_and_renor_chem_pot_solution_renor_chem_pot_guess);
    fprintf(file, "mass_and_renorm_chem_pot_solution_abs_error = %f\n", parameters.mass_and_renorm_chem_pot_solution_abs_error);
    fprintf(file, "mass_and_renorm_chem_pot_solution_rel_error = %f\n\n", parameters.mass_and_renorm_chem_pot_solution_rel_error);
    
    fprintf(file, "fermi_dirac_integrals_max_interval_num = %d\n", parameters.fermi_dirac_integrals_max_interval_num);
    fprintf(file, "fermi_dirac_integrals_integration_key = %d\n", parameters.fermi_dirac_integrals_integration_key);
    fprintf(file, "fermi_dirac_integrals_max_sub_interval = %d\n", parameters.fermi_dirac_integrals_max_sub_interval);
    fprintf(file, "fermi_dirac_integrals_abs_error = %f\n", parameters.fermi_dirac_integrals_abs_error);
    fprintf(file, "fermi_dirac_integrals_rel_error = %f\n\n", parameters.fermi_dirac_integrals_rel_error);
    

}

