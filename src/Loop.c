
#include <math.h>
#include <gsl/gsl_vector.h>

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "ZeroTemperatureEOS.h"
#include "FiniteTemperatureEOS.h"
#include "AuxiliaryFunctions.h"
#include "Derivatives.h"

/*  The EOS is solved in the function bellow. In this particular implementation,
    the barionic density was chosen as the running parameter.
 */
int SolveZeroTemperatureEOS(){

    // Print name of parametrization
    if (options.verbose){
        printf("Calculation performed with %s parameters set.\n",
			   parameters.parameters_set_identifier);
	}

	int points_number = parameters.loop_variable.points_number;

    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(points_number);

    /*
     * Vacuum mass determination
     */

    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double vacuum_mass = VacuumMassDetermination(); //TODO: move to a separate file?

    double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
    double bag_constant = ThermodynamicPotential(parameters.model.bare_mass, 0.0, 0.0, 0.0)
                          - vacuum_thermodynamic_potential;

    if (options.verbose){
        printf("\tVacuum mass: %f\n", vacuum_mass);
        printf("\tBag constant: %f\n", bag_constant);
    }

    /*
     * Main loop (on barionic density)
     */

	// Start with barionic density set at the minimum value
    double barionic_density = parameters.loop_variable.min_value;
    double density_step = Step(parameters.loop_variable.min_value,
							   parameters.loop_variable.max_value,
							   parameters.loop_variable.points_number);

    if (options.verbose)
		printf("Solving gap equation and equations of state ...\n");

    double barionic_density_at_zero_mass;
    double chemical_potential_at_zero_mass;
    bool chiral_restoration = false;

    for (int i = 0; i < parameters.loop_variable.points_number; i++){

        gsl_vector_set(barionic_density_vector, i, barionic_density);
        if (options.verbose){
            printf("\r\tBarionic density: %f", barionic_density);
            fflush(stdout);
        }

        // Determination of Fermi momentum
        double fermi_momentum = CONST_HBAR_C
								* pow(3.0 * pow(M_PI, 2.0) * barionic_density / NUM_FLAVORS,
									  1.0 / 3.0);

        // Solution of Gap Equation, determination of scalar density
        double mass = GapEquationSolver(fermi_momentum);
        gsl_vector_set(mass_vector, i, mass);

        double scalar_density = ScalarDensity(mass, fermi_momentum);
        gsl_vector_set(scalar_density_vector, i, scalar_density);

        // Determination of chemical potential
        double renormalized_chemical_potential = sqrt(pow(fermi_momentum, 2.0) + pow(mass, 2.0));
        gsl_vector_set(renormalized_chemical_potential_vector, i, renormalized_chemical_potential);

        double chemical_potential =	renormalized_chemical_potential
                                    + 2.0 * NUM_COLORS * parameters.model.G_V * CONST_HBAR_C * barionic_density;

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

        if (mass <= ZERO_MASS_TOL && chiral_restoration == false){
            barionic_density_at_zero_mass = barionic_density;
            chemical_potential_at_zero_mass = chemical_potential;
            chiral_restoration = true;
        }


        barionic_density += density_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop does't use new line, we need one now

    if (chiral_restoration)
        printf("Chiral restoration at density, chemical potential: %f, %f\n",
               barionic_density_at_zero_mass,
               chemical_potential_at_zero_mass);

    gsl_vector * energy_density_per_particle_vector =
		VectorNewVectorFromDivisionElementByElement(energy_density_vector,
                                            		barionic_density_vector);

	// Determination of chiral restoration for bare_mass != 0
	gsl_vector * mass_2nd_derivative_vector;
	gsl_vector * mass_derivative_vector;
	if (parameters.model.bare_mass != 0){

		mass_derivative_vector = gsl_vector_alloc(points_number);
		mass_2nd_derivative_vector = gsl_vector_alloc(points_number);

		for (int i = 0; i < points_number; i++){

			double derivative = DerivativeWithOrder(1,
													chemical_potential_vector,
													mass_vector,
													i);

			gsl_vector_set(mass_derivative_vector, i, derivative);
		}

		for (int i = 0; i < points_number; i++){

			// calculate $\partial^2 m / \partial \mu^2$
			double derivative = DerivativeWithOrder(2,
        	                   						chemical_potential_vector,
        	                   					    mass_vector,
        	                   					 	i);

			gsl_vector_set(mass_2nd_derivative_vector, i, derivative);
		}
	}

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

	if (parameters.model.bare_mass != 0){

		WriteVectorsToFile("der_mass_to_chem_pot.dat",
						   "# derivative of mass with respect"
						   "to the chemical potential\n"
						   "#\n# chemical_potential (MeV), derivative\n",
						   2,
						   chemical_potential_vector,
						   mass_derivative_vector);

		WriteVectorsToFile("sec_der_mass_to_chem_pot.dat",
						   "# second derivative of mass with respect"
						   "to the chemical potential\n"
						   "#\n# chemical_potential (MeV), derivative (MeV^{-1})\n",
						   2,
						   chemical_potential_vector,
						   mass_2nd_derivative_vector);
	}

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

    WriteVectorsToFile ("pressure_vs_energy_density.dat",
                        "# pressure (MeV/fm^3), energy density (MeV/fm^3)",
                        2,
                        pressure_vector,
                        energy_density_vector);

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
    if (options.verbose){
        printf("Calculation performed with %s parameters set.\n",
			   parameters.parameters_set_identifier);
	}

	int points_number = parameters.loop_variable.points_number;
    // Vectors to store results
    gsl_vector * barionic_density_vector = gsl_vector_alloc(points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(points_number);

    /*
     * Vacuum mass determination
     */

    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double vacuum_mass = VacuumMassDetermination(); // TODO: move to another file?

    double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
    double bag_constant = ThermodynamicPotentialForFiniteTemperature(parameters.model.bare_mass, 0.0, 0.0)
                          - vacuum_thermodynamic_potential;

    if (options.verbose){
        printf("\tVacuum mass: %f\n", vacuum_mass);
        printf("\tBag constant: %f\n", bag_constant);
        printf("\tTemperature: %f\n", parameters.finite_temperature.temperature);
    }

    /*
     * Main loop (on barionic density)
     */

    // Define the density step. We subtract 1 from the number of points to
    // make sure that the last point corresponds to parameters.maximum_density
    double barionic_density = parameters.loop_variable.min_value;
    double density_step = Step(parameters.loop_variable.min_value,
							   parameters.loop_variable.max_value,
							   parameters.loop_variable.points_number);

    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");

    for (int i = 0; i < points_number; i++){

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
                                    + 2.0 * NUM_COLORS * parameters.model.G_V * CONST_HBAR_C * barionic_density;

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
        double entropy = Entropy(mass,
								 parameters.finite_temperature.temperature,
								 renormalized_chemical_potential);
        double energy_density = EnergyForFiniteTemperature(regularized_thermodynamic_potential,
                                                           chemical_potential,
                                                           barionic_density,
                                                           parameters.finite_temperature.temperature,
                                                           entropy);

        gsl_vector_set(energy_density_vector, i, energy_density);

        barionic_density += density_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop doesn't use new line, we need one now

    // Calculate energy per particle
    gsl_vector * energy_density_per_particle_vector =
		VectorNewVectorFromDivisionElementByElement(energy_density_vector,
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
