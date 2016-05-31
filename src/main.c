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

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"
#include "AuxiliaryFunctions.h"

int PerformCalculation();

int main(int argc, char * argv[])
{
    CommandlineOptionsParse(argc, argv);
    ParametersSetup();
    
    PerformCalculation();
    
    return 0;
}

/*  The EOS is solved in the function bellow. In this particular implementation, 
    the barionic density was chosen as the running parameter.
 */
int PerformCalculation(){
    
    // Print name of parametrization
    if (options.verbose)
        printf("Calculation performed with %s parameters set.\n", parameters.parameters_set_identifier);
    
    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * fermi_momentum_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * vacuum_thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
    
    // Vacuum mass determination
    if (options.verbose)
        printf("Determining the vacuum mass ...\n");
    
    WriteVacuumMassEquation("data/vacuum_mass_equation.dat", 0.0, 1000.0, 1000);
    double vacuum_mass = VacuumMassDetermination();
    
    if (options.verbose)
        printf("\tVacuum mass: %f\n", vacuum_mass);
    
    // DEBUG
    {
        double minimum_mass = 0.0;
        double maximum_mass = 1000.0;
        int points_number = 1000;
        double chemical_potential = 410.0;
        
        double m = 0;
        
        double step = (maximum_mass - minimum_mass) / (points_number - 1);
        
        FILE * f = fopen("data/therm.dat", "w");
        
        if (NULL == f) {
            printf("Could not open for writting.\n");
            perror("Reason");
            exit(EXIT_FAILURE);
        }
        
        while (m <= maximum_mass) {

            double fermi_momentum = 0;
            if (pow(chemical_potential, 2.0) > pow(m, 2.0)){
                fermi_momentum = sqrt(pow(chemical_potential, 2.0) - pow(m, 2.0));
            }
            
            double barionic_density = pow(fermi_momentum, 3.0) * NUM_FLAVORS / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
            
            double vacuum_thermodynamic_potential = VacuumThermodynamicPotential(vacuum_mass, 0.0 * fermi_momentum);
            fprintf(f,
                    "%20.15E\t%20.15E\n",
                    m,
                    ThermodynamicPotential(m, fermi_momentum, barionic_density, chemical_potential, vacuum_thermodynamic_potential));
            m += step;
        }
        
        fclose(f);
        printf("...\n");
    }
    exit(0);

    
	// Define the density step. We subtract 1 from the number of points to
	// make sure that the last point corresponds to parameters.maximum_density
    double barionic_density = parameters.minimum_density;

    double density_step = (parameters.maximum_density - parameters.minimum_density) / (parameters.points_number - 1);

	if (options.verbose)
		printf("Solving gap equation and equations of state ...\n");

    for (int i = 0; i < parameters.points_number; i++){
        
        barionic_density += density_step;
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        if (options.verbose){
            printf("\r\tBarionic density: %f", barionic_density);
            fflush(stdout);
        }
        
        // Determination of Fermi momentum
        double fermi_momentum = CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * barionic_density / NUM_FLAVORS, 1.0 / 3.0);
        gsl_vector_set(fermi_momentum_vector, i, fermi_momentum);
        
        // Write zeroed gap equation for every value of density (for checking purpouses)
        char filename[256];
        sprintf(filename, "data/gap/gap_dens_%d.dat", i);
        WriteZeroedGapEquation(filename, parameters.gap_minimum_mass, parameters.gap_maximum_mass, parameters.gap_points_number, fermi_momentum);
		
        // Solution of Gap Equation, determination of scalar density
		double mass = GapEquationSolver(fermi_momentum);
        gsl_vector_set(mass_vector, i, mass);
        
        double scalar_density = ScalarDensity(mass, fermi_momentum);
        gsl_vector_set(scalar_density_vector, i, scalar_density);
        
        // Determination of chemical potential
        double chemical_potential = sqrt(pow(fermi_momentum, 2.0) + pow(mass, 2.0));
        gsl_vector_set(chemical_potential_vector, i, chemical_potential);
        
        // Determination of termodinamic potential
        double vacuum_thermodynamic_potential = VacuumThermodynamicPotential2(vacuum_mass, barionic_density, chemical_potential);
        double thermodynamic_potential = ThermodynamicPotential2(mass,
                                                                barionic_density,
                                                                fermi_momentum,
                                                                scalar_density,
                                                                chemical_potential,
                                                                vacuum_thermodynamic_potential);
 
 /*     double vacuum_thermodynamic_potential = VacuumThermodynamicPotential(vacuum_mass, fermi_momentum);
        double thermodynamic_potential = ThermodynamicPotential(mass,
                                                                fermi_momentum,
                                                                barionic_density,
                                                                chemical_potential,
                                                                vacuum_thermodynamic_potential);
*/
        gsl_vector_set(vacuum_thermodynamic_potential_vector, i, vacuum_thermodynamic_potential);
        gsl_vector_set(thermodynamic_potential_vector, i, thermodynamic_potential);
        
        // Determination of pressure
        double pressure = Pressure(thermodynamic_potential);
        gsl_vector_set(pressure_vector, i, pressure);
        
        // Determination of energy density
        double energy_density = EnergyDensity(thermodynamic_potential, chemical_potential, barionic_density);
        gsl_vector_set(energy_density_vector, i, energy_density);
        
    }
    
    // Write results
    if (options.verbose)
        printf("\n"                     // Printing in the loop does't use new line, so here we need one 
               "Saving results ...\n");

    WriteVectorsToFile("data/mass.dat",
                       "# barionic density, mass\n",
                       2,
                       barionic_density_vector,
                       mass_vector);
    
    WriteVectorsToFile("data/fermi_momentum.dat",
                       "# barionic density, fermi momentum\n",
                       2,
                       barionic_density_vector,
                       fermi_momentum_vector);
    
    WriteVectorsToFile("data/scalar_density.dat",
                       "# barionic density, scalar density \n",
                       2,
                       barionic_density_vector,
                       scalar_density_vector);
    
    WriteVectorsToFile("data/chemical_potential.dat",
                       "# barionic density, chemical potential \n",
                       2,
                       barionic_density_vector,
                       chemical_potential_vector);
    
    WriteVectorsToFile("data/vacuum_thermodynamic_potential.dat",
                       "# barionic density, vacuum thermodynamic potential \n",
                       2,
                       barionic_density_vector,
                       vacuum_thermodynamic_potential_vector);

    WriteVectorsToFile("data/thermodynamic_potential.dat",
                       "# barionic density, thermodynamic_potential \n",
                       2,
                       barionic_density_vector,
                       thermodynamic_potential_vector);

    WriteVectorsToFile("data/pressure.dat",
                       "# barionic density, pressure \n",
                       2,
                       barionic_density_vector,
                       pressure_vector);

    WriteVectorsToFile("data/energy_density.dat",
                       "# barionic density, scalar density \n",
                       2,
                       barionic_density_vector,
                       energy_density_vector);
    
    gsl_vector * energy_density_per_particle_vector = VectorNewVectorFromDivisionElementByElement(energy_density_vector, barionic_density_vector);
    WriteVectorsToFile("data/energy_density_per_particle.dat",
                       "# barionic density, energy density per particle \n",
                       2,
                       barionic_density_vector,
                       energy_density_per_particle_vector);
    
    WriteVectorsToFile("data/thermodynamic_potential_vs_mass.dat",
                       "# mass, thermodynamic potential\n",
                       2,
                       mass_vector,
                       thermodynamic_potential_vector);
    
    // Free memory associated with vectors
    gsl_vector_free(barionic_density_vector);
    gsl_vector_free(fermi_momentum_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(scalar_density_vector);
    gsl_vector_free(thermodynamic_potential_vector);
    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    
    if (options.verbose)
    	printf("Done!\n");

    return 0;
}
