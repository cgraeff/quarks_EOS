//
//  main.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>

#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"
#include "AuxiliaryFunctions.h"

int PerformCalculation();

int main(int argc, char * argv[])
{
    // ProcessCommandlineArguments();
    ParametersSetup();
    
    PerformCalculation();
    
    return 0;
}

/*  The EOS is solved in the function bellow. In this particular implementation, 
    the barionic density was chosen as the running parameter.
 */
int PerformCalculation(){
    
    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    
	// Define the density step. We subtract 1 from the number of points to
	// make sure that the last point corresponds to parameters.maximum_density
    double barionic_density = parameters.minimum_density;

    double density_step = (parameters.maximum_density - parameters.minimum_density) / (parameters.points_number - 1);

    for (int i = 0; i < parameters.points_number; i++){
        
        barionic_density += density_step;
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        
        double fermi_momentum = pow(3.0 * pow(M_PI, 2.0) * barionic_density / NUM_FLAVORS, 1.0 / 3.0);
        
        char filename[256];
        sprintf(filename, "data/gap/gap_dens_%d.dat", i);
        WriteGapEquation(filename, parameters.gap_minimum_mass, parameters.gap_maximum_mass, parameters.gap_points_number, fermi_momentum);
		
		double mass = GapEquationSolver(fermi_momentum);
        gsl_vector_set(mass_vector, i, mass);
        
        double scalar_density = ScalarDensity(mass, fermi_momentum);
        gsl_vector_set(scalar_density_vector, i, scalar_density);
    }
    
    // Write results
    WriteVectorsToFile("data/mass.dat",
                       "# barionic density, mass\n",
                       2,
                       barionic_density_vector,
                       mass_vector);
    WriteVectorsToFile("data/scalar_density.dat",
                       "# barionic density, mass \n",
                       2,
                       barionic_density_vector,
                       scalar_density_vector);
    
    // Free memory associated with vectors
    gsl_vector_free(barionic_density_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(scalar_density_vector);

    return 0;
}
