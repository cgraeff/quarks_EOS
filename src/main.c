//
//  main.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "Tests.h"
#include "ZeroTemperatureEOS.h"
#include "FiniteTemperatureEOS.h"
#include "AuxiliaryFunctions.h"

int SolveZeroTemperatureEOS();
int SolveFiniteTemperatureEOS();

int main(int argc, char * argv[])
{
    CommandlineOptionsParse(argc, argv);
    ParametersSetup();

    if (options.tests){
  		RunTests();
        exit(EXIT_SUCCESS);
    }

  	// If option -p is used, set parameters set accordingly,
	// otherwise, use default set
  	SetParametersSet(options.parameterization);
    
    // If the temperature was chosen using
    // commandline options, use it
    // (-1.0 is a place holder value)
    if (options.temp != -1.0)
        parameters.temperature = options.temp;
    
    if (parameters.temperature == 0){
        SolveZeroTemperatureEOS();
    }
    else if (parameters.temperature > 0){
        SolveFiniteTemperatureEOS();
    }
    else{
        printf("Values of temperature must be non-negative.\n");
        printf("(%f was provided).\n", parameters.temperature);
        exit(EXIT_FAILURE);
    }

    return 0;
}

/*  The EOS is solved in the function bellow. In this particular implementation, 
    the barionic density was chosen as the running parameter.
 */
int SolveZeroTemperatureEOS(){
    
    // Print name of parametrization
    if (options.verbose)
        printf("Calculation performed with %s parameters set.\n", parameters.parameters_set_identifier);
    
    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
    
    /*
     * Vacuum mass determination
     */
    
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");
    
    double vacuum_mass = VacuumMassDetermination();

    double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
    double bag_constant = ThermodynamicPotential(parameters.bare_mass, 0.0, 0.0, 0.0)
                          - vacuum_thermodynamic_potential;

    if (options.verbose){
        printf("\tVacuum mass: %f\n", vacuum_mass);
        printf("\tBag constant: %f\n", bag_constant);
    }

    /*
     * Main loop (on barionic density)
     */
    
	// Define the density step. We subtract 1 from the number of points to
	// make sure that the last point corresponds to parameters.maximum_density
    double barionic_density = parameters.minimum_density;

    double density_step = (parameters.maximum_density - parameters.minimum_density)
  						  / (parameters.points_number - 1);

    if (options.verbose)
		printf("Solving gap equation and equations of state ...\n");

    for (int i = 0; i < parameters.points_number; i++){
        
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        if (options.verbose){
            printf("\r\tBarionic density: %f", barionic_density);
            fflush(stdout);
        }

        // Determination of Fermi momentum
        double fermi_momentum = CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * barionic_density / NUM_FLAVORS, 1.0 / 3.0);

        // Solution of Gap Equation, determination of scalar density
		double mass = GapEquationSolver(fermi_momentum);
        gsl_vector_set(mass_vector, i, mass);

        double scalar_density = ScalarDensity(mass, fermi_momentum);
        gsl_vector_set(scalar_density_vector, i, scalar_density);

        // Determination of chemical potential
        double renormalized_chemical_potential = sqrt(pow(fermi_momentum, 2.0) + pow(mass, 2.0));
        gsl_vector_set(renormalized_chemical_potential_vector, i, renormalized_chemical_potential);

        double chemical_potential =	renormalized_chemical_potential
                                    + 2.0 * parameters.G_V * NUM_FLAVORS * NUM_COLORS * pow(fermi_momentum, 3.0)
                                      / (3.0 * pow(M_PI * CONST_HBAR_C, 2.0));
        gsl_vector_set(chemical_potential_vector, i, chemical_potential);

        double regularized_thermodynamic_potential = - vacuum_thermodynamic_potential
                                                     + ThermodynamicPotential(mass,
                                                                              fermi_momentum,
                                                                              chemical_potential,
                                                                              renormalized_chemical_potential);

        gsl_vector_set(thermodynamic_potential_vector,
					   i,
					   regularized_thermodynamic_potential);

        // Determination of pressure
        double pressure = Pressure(regularized_thermodynamic_potential);
        gsl_vector_set(pressure_vector, i, pressure);

        // Determination of energy density
        double energy_density = EnergyDensity(regularized_thermodynamic_potential,
											  chemical_potential,
											  barionic_density);
        gsl_vector_set(energy_density_vector, i, energy_density);
        
        barionic_density += density_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop does't use new line, we need one now
    
    gsl_vector * energy_density_per_particle_vector = VectorNewVectorFromDivisionElementByElement(energy_density_vector,
                                                                                                  barionic_density_vector);
    /*
     * Save results
     */
    
    if (options.verbose)
        printf("Saving results ...\n");

    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile("mass.dat",
                       "# barionic density, mass\n",
                       2,
                       barionic_density_vector,
                       mass_vector);

    WriteVectorsToFile("scalar_density.dat",
                       "# barionic density, scalar density \n",
                       2,
                       barionic_density_vector,
                       scalar_density_vector);
    
    WriteVectorsToFile("renormalized_chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       renormalized_chemical_potential_vector);

    WriteVectorsToFile("chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       chemical_potential_vector);

    WriteVectorsToFile("thermodynamic_potential.dat",
                       "# barionic density, thermodynamic_potential \n",
                       2,
                       barionic_density_vector,
                       thermodynamic_potential_vector);
    
    WriteVectorsToFile("thermodynamic_potential_vs_mass.dat",
                       "# mass, thermodynamic potential\n",
                       2,
                       mass_vector,
                       thermodynamic_potential_vector);

    if (options.dirs)
        SetFilePath("output/Rep_Tsue/data");

    WriteVectorsToFile("mass_renorm_chem_pot.dat",
                       "# mass, renormalized chemical potential \n",
                       2,
                       mass_vector,
                       renormalized_chemical_potential_vector);

    WriteVectorsToFile("mass_chem_pot.dat",
                       "# mass, chemical potential \n",
                       2,
                       mass_vector,
                       chemical_potential_vector);

    WriteVectorsToFile("pressure_renorm_chem_pot.dat",
                       "# renormalized chemical potential, pressure\n",
                       2,
                       renormalized_chemical_potential_vector,
                       pressure_vector);

    WriteVectorsToFile("pressure_chem_pot.dat",
                       "# chemical potential, pressure\n",
                       2,
                       chemical_potential_vector,
                       pressure_vector);

    WriteVectorsToFile("barionic_dens_renorm_chemical_potential.dat",
                       "# barionic density, renormalized chemical potential \n",
                       2,
                       barionic_density_vector,
                       renormalized_chemical_potential_vector);

    WriteVectorsToFile("barionic_dens_chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       chemical_potential_vector);

    if (options.dirs)
        SetFilePath("output/EOS/data/");
    
    WriteVectorsToFile("pressure.dat",
                       "# barionic density, pressure \n",
                       2,
                       barionic_density_vector,
                       pressure_vector);

    WriteVectorsToFile("energy_density.dat",
                       "# barionic density, energy density \n",
                       2,
                       barionic_density_vector,
                       energy_density_vector);

    WriteVectorsToFile("energy_density_per_particle.dat",
                       "# barionic density, energy density per particle \n",
                       2,
                       barionic_density_vector,
                       energy_density_per_particle_vector);
    
    SetFilePath(NULL);

    /*
     * Clean up
     */

    // Free memory associated with vectors
    gsl_vector_free(barionic_density_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(scalar_density_vector);
    gsl_vector_free(renormalized_chemical_potential_vector);
    gsl_vector_free(chemical_potential_vector);
    gsl_vector_free(thermodynamic_potential_vector);
    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);


    if (options.verbose)
    	printf("Done!\n");

    return 0;
}

int SolveFiniteTemperatureEOS(){

    // Print name of parametrization
    if (options.verbose)
        printf("Calculation performed with %s parameters set.\n", parameters.parameters_set_identifier);
    
    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
    
    /*
     * Vacuum mass determination
     */
    
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");
    
    double vacuum_mass = VacuumMassDetermination();
    
    double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
    double bag_constant = ThermodynamicPotentialForFiniteTemperature(parameters.bare_mass, 0.0, 0.0)
                          - vacuum_thermodynamic_potential;
    
    if (options.verbose){
        printf("\tVacuum mass: %f\n", vacuum_mass);
        printf("\tBag constant: %f\n", bag_constant);
        printf("\tTemperature: %f\n", parameters.temperature);
    }
    
    /*
     * Main loop (on barionic density)
     */
    
    // Define the density step. We subtract 1 from the number of points to
    // make sure that the last point corresponds to parameters.maximum_density
    double barionic_density = parameters.minimum_density;
    
    double density_step = (parameters.maximum_density - parameters.minimum_density)
    / (parameters.points_number - 1);
    
    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");
    
    for (int i = 0; i < parameters.points_number; i++){
        
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        
        if (options.verbose){
            printf("\r\tBarionic density: %f", barionic_density);
            fflush(stdout);
        }
        
        // Solution of Gap Equation, determination of scalar density
        double mass;
        double renormalized_chemical_potential;
        SolveMultiRoots(barionic_density,
                        &mass,
                        &renormalized_chemical_potential);
        gsl_vector_set(mass_vector, i, mass);
        gsl_vector_set(renormalized_chemical_potential_vector, i, renormalized_chemical_potential);
        
        double chemical_potential =	renormalized_chemical_potential
                                    + 2.0 * NUM_COLORS * parameters.G_V * CONST_HBAR_C * barionic_density;
        
        gsl_vector_set(chemical_potential_vector, i, chemical_potential);
        
        double thermodynamic_potential = ThermodynamicPotentialForFiniteTemperature(mass,
                                                                                    chemical_potential,
                                                                                    renormalized_chemical_potential);
        double regularized_thermodynamic_potential = thermodynamic_potential - vacuum_thermodynamic_potential;
        
        // Just the regularized thermodynamic potential will be saved
        // as it's the only one that is used in other calculations
        gsl_vector_set(thermodynamic_potential_vector,
                       i,
                       regularized_thermodynamic_potential);
        
        // Determination of pressure
        double pressure = PressureForFiniteTemperature(regularized_thermodynamic_potential);
        gsl_vector_set(pressure_vector, i, pressure);
        
        // Determination of energy density
        double entropy = Entropy(mass, parameters.temperature, renormalized_chemical_potential);
        double energy_density = EnergyForFiniteTemperature(regularized_thermodynamic_potential,
                                                           chemical_potential,
                                                           barionic_density,
                                                           parameters.temperature,
                                                           entropy);
        
        gsl_vector_set(energy_density_vector, i, energy_density);
        
        barionic_density += density_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop doesn't use new line, we need one now

    // Calculate energy per particle
    gsl_vector * energy_density_per_particle_vector = VectorNewVectorFromDivisionElementByElement(energy_density_vector,
                                                                                                  barionic_density_vector);
    /*
     * Save results
     */

    if (options.verbose)
        printf("Saving results ...\n");
    
    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile("mass.dat",
                       "# barionic density, mass\n",
                       2,
                       barionic_density_vector,
                       mass_vector);

    WriteVectorsToFile("renormalized_chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       renormalized_chemical_potential_vector);

    WriteVectorsToFile("chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       chemical_potential_vector);
    
    WriteVectorsToFile("thermodynamic_potential.dat",
                       "# barionic density, thermodynamic_potential \n",
                       2,
                       barionic_density_vector,
                       thermodynamic_potential_vector);

    WriteVectorsToFile("thermodynamic_potential_vs_mass.dat",
                       "# mass, thermodynamic potential\n",
                       2,
                       mass_vector,
                       thermodynamic_potential_vector);
    
    if (options.dirs)
        SetFilePath("output/EOS/data/");

    WriteVectorsToFile("pressure.dat",
                       "# barionic density, pressure \n",
                       2,
                       barionic_density_vector,
                       pressure_vector);

    WriteVectorsToFile("energy_density.dat",
                       "# barionic density, energy density \n",
                       2,
                       barionic_density_vector,
                       energy_density_vector);

    WriteVectorsToFile("energy_density_per_particle.dat",
                       "# barionic density, energy density per particle \n",
                       2,
                       barionic_density_vector,
                       energy_density_per_particle_vector);

    SetFilePath(NULL);

    /*
     * Clean up
     */

    // Free memory associated with vectors
    gsl_vector_free(barionic_density_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(thermodynamic_potential_vector);
    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    gsl_vector_free(energy_density_per_particle_vector);


    if (options.verbose)
        printf("Done!\n");
    
    return 0;
}
