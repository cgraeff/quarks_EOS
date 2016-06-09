//
//  Tests.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include "Tests.h"

int RunTests()
{
	WriteVacuumMassEquation("data/vacuum_mass_equation.dat", 0.0, 1000.0, 1000);

    WriteZeroedRenormalizedChemicalPotentialEquation("data/zeroed_renorm_chemical_pot_equation.dat",
                                                     0.0,
                                                     1000.0,
                                                     1000,
                                                     0.0,
                                                     0.0);

    //    double chemical_potential[4] = {0.0, 300.0, 368.6, 400.0};
    double chemical_potential[4] = {0.0, 430, 440, 444.3};
    for (int i = 0; i < 4; i++){

        double minimum_mass = 0.0;
        double maximum_mass = 1000.0;
        int points_number = 1000;

        double m = 0;

        double step = (maximum_mass - minimum_mass) / (points_number - 1);

        char filename[256];
        sprintf(filename, "data/therm_%d.dat", i);
        FILE * f = fopen(filename, "w");

        if (NULL == f) {
            printf("Could not open for writting.\n");
            perror("Reason");
            exit(EXIT_FAILURE);
        }

	  	// Prepare input for ZeroedGapEquation
    	renorm_chem_pot_equation_input input;
    	input.chemical_potential = chemical_potential[i];
        input.mass = m;

    	// Prepare function to be passed to the root finding algorithm
    	gsl_function F;
    	F.function = &ZeroedRenormalizedChemicalPotentialEquation;
    	F.params = &input;

        while (m <= maximum_mass) {

			double renormalized_chemical_potential = chemical_potential[i];

		    if (parameters.G_V != 0.0){

		    	renormalized_chemical_potential =
			  			UnidimensionalRootFinder(&F,
												 parameters.renormalized_chemical_potential_lower_bound,
												 parameters.renormalized_chemical_potential_upper_bound,
												 parameters.renormalized_chemical_potential_abs_error,
												 parameters.renormalized_chemical_potential_rel_error,
												 parameters.renormalized_chemical_potential_max_iter);
			}

            double fermi_momentum = 0;
            if (pow(renormalized_chemical_potential, 2.0) > pow(m, 2.0)){
                fermi_momentum = sqrt(pow(renormalized_chemical_potential, 2.0)
			  						  - pow(m, 2.0));
            }

            fprintf(f,
                    "%20.15E\t%20.15E\n",
                    m,
                    -vacuum_thermodynamic_potential
                    + ThermodynamicPotential(m,
					 					     fermi_momentum,
										     chemical_potential[i],
										     renormalized_chemical_potential));
            m += step;
        }

        fclose(f);
        printf("...\n");
    }

    // Write zeroed gap equation for every value of density (for checking purpouses)
    char filename[256];
    sprintf(filename, "data/gap/gap_dens_%d.dat", i);
    WriteZeroedGapEquation(filename,
						   parameters.gap_minimum_mass,
						   parameters.gap_maximum_mass,
						   parameters.gap_points_number,
						   fermi_momentum);
}
