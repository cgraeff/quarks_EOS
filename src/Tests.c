//
//  Tests.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include "Tests.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include "ZeroTemperatureEOS.h"
#include "FiniteTemperatureEOS.h"
#include "Parameters.h"
#include "CommandlineOptions.h"
#include "AuxiliaryFunctions.h"
#include "Constants.h"
#include "Maps.h"
#include "RootFinding.h"

double ZeroedGapEquationForFiniteTemperatureTest(double mass, void * p);

void RunTests()
{
#pragma mark Vacuum Mass Zeroed Equation

    if (true)
    {
        printf("\tVacuum Mass Zeroed Equation\n");
        SetFilePath("tests/vacuum-mass-equation/data/");

        SetParametersSet(options.parameterization);

		double min_mass = 0.0;
		double max_mass = 1000.0;
		int points_number = 1000;

    	double step = Step(min_mass, max_mass, points_number);

    	FILE * f = OpenFile("vacuum_mass_equation.dat");

		double m = 0;
    	while (m < points_number) {
        	fprintf(f, "%20.15E\t%20.15E\n", m, VacuumMassEquation(m, NULL));
        	m += step;
    	}

    	fclose(f);

        SetFilePath("tests/vacuum-mass-equation/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
				"The following tests were executed for %s parameterization:\n",
				parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tVacuum mass equation calculation.\n\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark ZeroedGapEquation
	if (true)
    {
        printf("\tVacuum Mass Zeroed Equation\n");
        SetFilePath("tests/vacuum-mass-equation/data/");

        SetParametersSet(options.parameterization);

		double min_mass = 0.0;
		double max_mass = 1000.0;
		int points_number = 1000;

    	double step = Step(min_mass, max_mass, points_number);
        
        double barionic_density[4] = {parameters.loop_variable.min_value,
                                      parameters.loop_variable.max_value * 1.0 / 3.0,
                                      parameters.loop_variable.max_value * 2.0 / 3.0,
                                      parameters.loop_variable.max_value};
        for (int i = 0; i < 4; i++){
            double fermi_momentum = CONST_HBAR_C
                                    * pow(3.0 * pow(M_PI, 2.0) * barionic_density[i] / NUM_FLAVORS,
                                          1.0 / 3.0);
            char filename[256];
            sprintf(filename, "zeroed_gap_equation_%d.dat", i);

    		FILE * f = OpenFile(filename);

    		gap_equation_input input;
    		input.fermi_momentum = fermi_momentum;

			double m = 0;
    		while (m < max_mass) {
        		fprintf(f, "%20.15E\t%20.15E\n", m, ZeroedGapEquation(m, &input));
        		m += step;
    		}
        }

        SetFilePath("tests/vacuum-mass-equation/");
        FILE * log_file = OpenFile("run.log");
        
        fprintf(log_file,
				"The following tests were executed for %s parameterization:\n",
				parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tZeroed gap equation calculation\n\n");
        
        PrintParametersToFile(log_file);
        
        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Zeroed Renormalized Chemical Potential
	if (true)
    {
        printf("\tVacuum Mass Zeroed Equation\n");
        SetFilePath("tests/vacuum-mass-equation/data/");

        SetParametersSet(options.parameterization);

		double mass = 0;
		double chemical_potential = 0;

		double min_renormalized_chemical_potential = 0.0;
		double max_renormalized_chemical_potential = 1000.0;
		int points_number = 1000;

    	double step = Step(min_renormalized_chemical_potential,
						   max_renormalized_chemical_potential,
						   points_number);

    	FILE * f = OpenFile("zeroed_renorm_chemical_pot_equation.dat");

    	renorm_chem_pot_equation_input input;
    	input.chemical_potential = chemical_potential;
    	input.mass = mass;

		double mu = 0;
    	while (mu < points_number) {
        	fprintf(f,
        	        "%20.15E\t%20.15E\n",
        	        mu,
        	        ZeroedRenormalizedChemicalPotentialEquation(mu, &input));
        	mu += step;
    	}

    	fclose(f);

        SetFilePath("tests/vacuum-mass-equation/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
				"The following tests were executed for %s parameterization:\n",
				parameters.parameters_set_identifier);
        fprintf(log_file,
				"\tZeroed renormalized chemical potential equation calculation\n\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611
    
    // Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611 (1996) 393-408
    // (the figure uses parameters of Set II from the article)
    if (true)
    {
        printf("\tReproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611\n");
        SetFilePath("tests/Buballa-Fig1-R/data/");
        
        SetParametersSet("Buballa_2");

        double chemical_potential[4] = {0.0, 350.0, 378.5, 410.0};

        double vacuum_mass = VacuumMassDetermination();

        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);

        for (int i = 0; i < 4; i++){

            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;

            double step = (maximum_mass - minimum_mass) / (points_number - 1);

            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);

            double m = 0;

            for (int j = 0; j < points_number; j++) {

                double renormalized_chemical_potential = chemical_potential[i];

                double fermi_momentum = 0;
                if (pow(renormalized_chemical_potential, 2.0) > pow(m, 2.0)){
                    fermi_momentum = sqrt(pow(renormalized_chemical_potential, 2.0) - pow(m, 2.0));
                }
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += step;
            }

            char filename[256];
            sprintf(filename, "Fig1_Buballa_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);
        }
        
        SetFilePath("tests/Buballa-Fig1-R/");
        FILE * log_file = OpenFile("run.log");
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file, "\tCalculation of vacuum mass, resultin in %f MeV.\n", vacuum_mass);
        fprintf(log_file,
                "\tCalculation of the thermodynamic potential as function of mass "
                "(Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611 (1996) 393-408).\n"
                "\tfiles: tests/data/Fig1_Buballa_*.dat\n");
        
        PrintParametersToFile(log_file);
        
        fclose(log_file);

        SetFilePath(NULL);
    }
    
#pragma mark Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports

    // Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = 0)
    if (true)
    {
        printf("\tReproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports\n");

        SetFilePath("tests/Buballa-Fig2.8-L/data/");

        SetParametersSet("BuballaR_2");
        
		int points_number = 1000;

        parameters.model.bare_mass = 0; // In the reference, the bare mass was
										// zero for this test

		double mass = 0;
        
        double chemical_potential[4] = {0.0, 300.0, 368.6, 400.0};

		double min_renormalized_chemical_potential = 0.0;
		double max_renormalized_chemical_potential = 1000;

		double min_mass = 0.0;
        double max_mass = 1000.0;

		double renorm_chem_pot_step = Step(min_renormalized_chemical_potential,
						   				   max_renormalized_chemical_potential,
						   				   points_number);

		double mass_step = Step(min_mass, max_mass, points_number);
        
        double vacuum_mass = VacuumMassDetermination();

        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);

        for (int i = 0; i < 4; i++){

            char filename_1[256];
            sprintf(filename_1, "ZeroedRenormalizedChemPotEquation_BR2L_%d.dat", i);

    		FILE * f = OpenFile(filename_1);

    		renorm_chem_pot_equation_input input;
    		input.chemical_potential = chemical_potential[i];
    		input.mass = mass;

			double mu = 0;
    		while (mu < points_number) {
        		fprintf(f,
           	     		"%20.15E\t%20.15E\n",
                		mu,
                		ZeroedRenormalizedChemicalPotentialEquation(mu, &input));
        		mu += renorm_chem_pot_step;
   	 		}

    		fclose(f);

            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);
            
            double m = 0;
            
            for (int j = 0; j < points_number; j++) {
                
                double renormalized_chemical_potential = chemical_potential[i];
                
                double fermi_momentum = 0;
                if (pow(renormalized_chemical_potential, 2.0) > pow(m, 2.0)){
                    fermi_momentum = sqrt(pow(renormalized_chemical_potential, 2.0)
                                          - pow(m, 2.0));
                }
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += mass_step;
            }
            
            char filename[256];
            sprintf(filename, "Fig2.8L_BuballaR_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);
  
        }
        
        SetFilePath("tests/Buballa-Fig2.8-L/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of the hermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:tests/data/Fig2.8L_BuballaR_*.dat\n");

        PrintParametersToFile(log_file);
        
        fclose(log_file);

        SetFilePath(NULL);
    }
    
#pragma mark Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports
    
    // Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = G_S)
    if (true)
    {
        printf("\tReproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports\n");
        SetFilePath("tests/Buballa-Fig2.8-R/data/");
        
        SetParametersSet("BuballaR_2_GV");
        
		int points_number = 1000;

        parameters.model.bare_mass = 0; // In the reference, the bare mass was
										// zero for this test
        
        double chemical_potential[4] = {0.0, 430.0, 440.0, 444.3};

        double mass[4] = {1.0, 100, 300, 700};

		double min_renorm_chemical_potential = 0.0;
		double max_renorm_chemical_potential = 1000.0;

		double step = Step(min_renorm_chemical_potential,
						   max_renorm_chemical_potential,
						   points_number);

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){
                char filename_1[256];
                sprintf(filename_1, "ZeroedRenormalizedChemPotEquation_BR2R_%d_%d.dat", i, j);

    			FILE * f = OpenFile(filename_1);

    			renorm_chem_pot_equation_input input;
    			input.chemical_potential = chemical_potential[i];
    			input.mass = mass[j];

				double mu = 0;
    			while (mu < points_number) {
        			fprintf(f,
                			"%20.15E\t%20.15E\n",
                			mu,
                			ZeroedRenormalizedChemicalPotentialEquation(mu, &input));
        			mu += step;
    			}

    			fclose(f);
            }
        }
        
        double vacuum_mass = VacuumMassDetermination();
        
        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
        
        for (int i = 0; i < 4; i++){
            
            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;
            
            double step = (maximum_mass - minimum_mass) / (points_number - 1);
            
            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);
            gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(points_number);
            
            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            renorm_chem_pot_equation_input input;
            F.function = &ZeroedRenormalizedChemicalPotentialEquation;
            F.params = &input;
     
            /*
             * Determine the renormalized chemical potential by root-finding
             * for each value of mass
             */
            
			RenormalizedChemicalPotentialIntegrationParameters params =
				parameters.renormalized_chemical_potential_integration;

            double m = 0;
            for (int j = 0; j < points_number; j++) {
                
                double renormalized_chemical_potential;

                // Prepare input for ZeroedRenormalizedChemicalPotentialEquation
                input.chemical_potential = chemical_potential[i];
                input.mass = m;
                
                // If the chemical potential is zero, the solution is zero.
                // Otherwise it needs to be calculated
                if (chemical_potential[i] == 0){
                    renormalized_chemical_potential = 0;
                }
                else{
                    int status = UnidimensionalRootFinder(&F,
                                                          params.lower_bound,
                                                          params.upper_bound,
                                                          params.abs_error,
                                                          params.rel_error,
                                                          params.max_iter,
                                                          &renormalized_chemical_potential);
                    if (status == -1)
                        renormalized_chemical_potential = 0;
                }
                
                gsl_vector_set(renormalized_chemical_potential_vector,
							   j,
							   renormalized_chemical_potential);
                
                double fermi_momentum = 0;
                if (pow(renormalized_chemical_potential, 2.0) > pow(m, 2.0)){
                    fermi_momentum = sqrt(pow(renormalized_chemical_potential, 2.0)
                                          - pow(m, 2.0));
                }
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += step;
            }
            
            char filename[256];
            sprintf(filename, "Fig2.8R_BuballaR_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);
            
            sprintf(filename, "renormalized_chemical_potential_%d.dat", i);
            WriteVectorsToFile(filename,
                               "#mass, renormalized_chemical_potential\n",
                               2,
                               mass_vector,
                               renormalized_chemical_potential_vector);
            
            gsl_vector_free(output);
            gsl_vector_free(mass_vector);
            gsl_vector_free(renormalized_chemical_potential_vector);
 
        }

        SetFilePath("tests/Buballa-Fig2.8-R/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of the thermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:tests/data/Fig2.8R_BuballaR_*.dat\n");
        fprintf(log_file,
                "\tZeroed renormalized chemical potential equation for selected parameters"
				" (as function of renormalized chemical potential)\n");
        fprintf(log_file,
                "\tRenormalized chemical potential (as function of mass)\n");
        fprintf(log_file,
                "\tVacuum mass: %f\n", vacuum_mass);
        fprintf(log_file,
                "\tVacuum thermodynamic potential: %f\n", vacuum_thermodynamic_potential);
        
        PrintParametersToFile(log_file);
        
        fclose(log_file);

        SetFilePath(NULL);
    }
    
#pragma mark Fermi-Dirac Distributions
    
    // Prints Fermi-Dirac distributions for selected values of temperature
    // and chemical potential as function of momentum
    if (true)
    {
        printf("\tFermi-Dirac Distributions\n");
        SetFilePath("tests/Fermi-Dirac-distributions/data/");
        
        int n_pts = 1000;
        
        double min_momentum = 0.0;
        double max_momentum = 400.0;
        double step = (max_momentum - min_momentum) / (n_pts - 1);
        
        double temperature[4] = {5, 10.0, 20.0, 30.0};
        double chemical_potential[4] = {50.0, 100.0, 200.0, 300.0};
        double mass[4] = {50.0, 100.0, 150.0, 200.0};
        
        gsl_vector * momentum_vector = gsl_vector_alloc(n_pts);
        gsl_vector * fp_vector = gsl_vector_alloc(n_pts);
        gsl_vector * fa_vector = gsl_vector_alloc(n_pts);
        
        for (int i = 0; i < 4; i++){ // temperature
            for (int j = 0; j < 4; j++){ // chemical potential
                for (int k = 0; k < 4; k++){ // mass
                    
                    double momentum = min_momentum;
                    
                    for (int l = 0; l < n_pts; l++){

                        double energy = sqrt(pow(momentum, 2.0) + pow(mass[k], 2.0));
                        
                        double fp = FermiDiracDistributionForParticles(energy,
																	   chemical_potential[j],
																	   temperature[i]);

                        double fa = FermiDiracDistributionForAntiparticles(energy,
																		   chemical_potential[j],
																		   temperature[i]);
                        
                        gsl_vector_set(momentum_vector, l, momentum);
                        gsl_vector_set(fp_vector, l, fp);
                        gsl_vector_set(fa_vector, l, fa);
                        
                        momentum += step;
                    }

                    char filename[256];
                    sprintf(filename, "FD_%d_%d_%d.dat", i, j, k);
                    
                    WriteVectorsToFile(filename,
                                       "# momentum, Fermi-Dirac distribution for particles, = for antiparticles\n",
                                       3,
                                       momentum_vector,
                                       fp_vector,
                                       fa_vector);
                }
            }
        }
        
        gsl_vector_free(momentum_vector);
        gsl_vector_free(fp_vector);
        gsl_vector_free(fa_vector);
        
        SetFilePath("tests/Fermi-Dirac-distributions/");
        FILE * log_file = OpenFile("run.log");
        
        fprintf(log_file,
                "Prints Fermi-Dirac distributions for selected values of temperature\n"
                "and chemical potential as function of momentum\n");
        
        PrintParametersToFile(log_file);

        fclose(log_file);
        
        SetFilePath(NULL);
    }
    
#pragma mark Fermi-Dirac Distribution Integrals
    
    // Tests integration of Fermi-Dirac distributions
    if (true)
    {
        printf("\tFermi-Dirac Distribution Integrals\n");
        SetFilePath("tests/Fermi-Dirac-distrib-integrals/data/");
        
        SetParametersSet("BuballaR_2");
        
        int num_points = 1000;
        
        double min_mass = 0.01;
        double max_mass = 1000;
        double step = (max_mass - min_mass) / ((double)(num_points - 1));
        
        double renormalized_chemical_potential[4] = {100.0, 200.0, 400.0, 600.0};
        double temperature[4] = {5.0, 10.0, 15.0, 20.0};

        gsl_vector * mass_vector = gsl_vector_alloc(num_points);
        gsl_vector * fermi_dirac_int_1 = gsl_vector_alloc(num_points);
        gsl_vector * fermi_dirac_int_2 = gsl_vector_alloc(num_points);
        
        for (int i = 0; i < 4; i++){
            
            parameters.finite_temperature.temperature = temperature[i];

            for (int j = 0; j < 4; j++){
                
                double m = min_mass;
                
                for (int k = 0; k < num_points; k++){
                    
                    double int_1 =
						FermiDiracDistributionFromDensityIntegral(m,
																  renormalized_chemical_potential[j]);

                    double int_2 =
						FermiDiracDistributionIntegralFromGapEquation(m,
																	  renormalized_chemical_potential[j]);
                
                    gsl_vector_set(mass_vector, k, m);
                    gsl_vector_set(fermi_dirac_int_1, k, int_1);
                    gsl_vector_set(fermi_dirac_int_2, k, int_2);
                    
                    m += step;
                    
                }
                
                char filename[256];
                sprintf(filename,
                        "fermi_dirac_distribution_from_density_integral_%d_%d.dat",
                        i,
                        j);
                
                WriteVectorsToFile(filename,
                                   "# mass, integral of fermi dirac dist from density\n",
                                   2,
                                   mass_vector,
                                   fermi_dirac_int_1);
                
                sprintf(filename,
                        "fermi_dirac_distribution_from_gap_eq_integral_%d_%d.dat",
                        i,
                        j);

                WriteVectorsToFile(filename,
                                   "# mass, integral of fermi dirac dist from gap eq\n",
                                   2,
                                   mass_vector,
                                   fermi_dirac_int_2);
            }
        }

        gsl_vector_free(mass_vector);
        gsl_vector_free(fermi_dirac_int_1);
        gsl_vector_free(fermi_dirac_int_2);
        
        SetFilePath("tests/Fermi-Dirac-distrib-integrals/");

        FILE * log_file = OpenFile("run.log");
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of integrals of fermi distributions as function of mass.\n");
        
        PrintParametersToFile(log_file);
        
        fclose(log_file);

        SetFilePath(NULL);
    }
    
#pragma mark Mass Gap Zeroed Equation for Selected Temperatures
    
    // writes gap equation as function of mass for selected temperatures
    //    What exactly is that? I don't remember ...
    //
    if (true)
	{
        printf("\tMass Gap Zeroed Equation for Selected Temperatures\n");
        SetFilePath("tests/zeroed-gap-equation/data/");
        
	  	SetParametersSet("BuballaR_2");

	  	int n_pts = 1000;

	  	double min_mass = 0.0;
	  	double max_mass = 1000.0;
	  	double step = (max_mass - min_mass) / (n_pts - 1);

	  	double temperature[4] = {5.0, 10.0, 15.0, 20.0};
	  	double renormalized_chemical_potential[4] = {100.0, 200.0, 400.0, 600.0};

	  	gsl_vector * m_vector = gsl_vector_alloc(n_pts);
	  	gsl_vector * zeroed_gap_vector = gsl_vector_alloc(n_pts);


	  	for (int i = 0; i < 4; i++){

		  	parameters.finite_temperature.temperature = temperature[i];

			for (int j = 0; j < 4; j++){

                double m = min_mass;

				for (int k = 0; k < n_pts; k++){

			  		double integ =
						FermiDiracDistributionIntegralFromGapEquation(m,
																	  renormalized_chemical_potential[j]);

			 		gsl_vector_set(zeroed_gap_vector,
								   k,
								   m - parameters.model.bare_mass - integ);

				  	gsl_vector_set(m_vector, k, m);

				  	m += step;
				}

				char filename[256];
			  	sprintf(filename, "zeroed_gap_eq_%d_%d.dat", i,	j);

				WriteVectorsToFile(filename,
								   "# mass, zeroed gap equation\n",
								   2,
								   m_vector,
								   zeroed_gap_vector);
			}
	  	}

	  	gsl_vector_free(m_vector);
	  	gsl_vector_free(zeroed_gap_vector);

        SetFilePath("tests/zeroed-gap-equation/");
        FILE * log_file = OpenFile("run.log");

		fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of zeroed gap eq for T != 0 as function of mass.\n");
        
        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
	}
    
#pragma mark Maps of Mass and Renormalized Chemical Potential Zeroed Equations
    
    // Calculates zeroed gap and barionic densities equations so we can see both
    // and have an insight of what's going on
    if (true)
    {
        printf("\tMaps of Mass and Renormalized Chemical Potential Zeroed Equations\n");
        SetFilePath("tests/maps/data/");
        
        SetParametersSet("BuballaR_2");

        const int num_densities = 10;
        const int num_temperatures = 10;

        const double barionic_density[10] = {0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 0.44};
        const double temperature[10] = {1.0, 3.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0};

        int mass_n_pts = 150;
        int renorm_chem_pot_n_pts = 150;
        double min_mass = 0.0;
        double max_mass = 600.0;
        double min_renormalized_chemical_potential = 0.0;
        double max_renormalized_chemical_potential = 600.0;

        double tolerance_dens = 0.05;
        double tolerance_gap = 0.5;

        gsl_vector * map_gap_x = gsl_vector_alloc(mass_n_pts * renorm_chem_pot_n_pts);
        gsl_vector * map_gap_y = gsl_vector_alloc(mass_n_pts * renorm_chem_pot_n_pts);
        int map_gap_num_points;

        gsl_vector * map_dens_x = gsl_vector_alloc(mass_n_pts * renorm_chem_pot_n_pts);
        gsl_vector * map_dens_y = gsl_vector_alloc(mass_n_pts * renorm_chem_pot_n_pts);
        int map_dens_num_points;

        for (int i = 0; i < num_temperatures; i++){ // Temperature

            parameters.finite_temperature.temperature = temperature[i];

            for (int j = 0; j < num_densities; j++){

                MapFunction(&ZeroedGapEquationForFiniteTemperature,
                            min_mass,
                            max_mass,
                            mass_n_pts,
                            min_renormalized_chemical_potential,
                            max_renormalized_chemical_potential,
                            renorm_chem_pot_n_pts,
                            tolerance_gap,
                            NULL,
                            false,
                            map_gap_x,
                            map_gap_y,
                            &map_gap_num_points);

                MapFunction(&ZeroedBarionicDensityEquationForFiniteTemperature,
                            min_mass,
                            max_mass,
                            mass_n_pts,
                            min_renormalized_chemical_potential,
                            max_renormalized_chemical_potential,
                            renorm_chem_pot_n_pts,
                            tolerance_dens,
                            (void *)&(barionic_density[j]),
                            false,
                            map_dens_x,
                            map_dens_y,
                            &map_dens_num_points);

                double x_intersection;
                double y_intersection;

                IntersectionPointOfTwoMaps(map_gap_x,
                                           map_gap_y,
                                           map_gap_num_points,
                                           map_dens_x,
                                           map_dens_y,
                                           map_dens_num_points,
                                           &x_intersection,
                                           &y_intersection);

                char filename[256];
                sprintf(filename, "map_gap_%d_%d.dat", i, j);
                WriteVectorsToFileUpToIndex(filename,
                                            "# Map of the region of zeroed gap equation that is near zero\n"
                                            "# mass, renormalized chemical potential\n",
                                            map_gap_num_points,
                                            2,
                                            map_gap_x,
                                            map_gap_y);

                sprintf(filename, "map_dens_%d_%d.dat", i, j);
                WriteVectorsToFileUpToIndex(filename,
                                            "# Map of the region of zeroed density gap equation that is near zero\n"
                                            "# mass, renormalized chemical potential\n",
                                            map_dens_num_points,
                                            2,
                                            map_dens_x,
                                            map_dens_y);

                sprintf(filename, "intersection_%d_%d.dat", i, j);
                FILE * file = OpenFile(filename);
                fprintf(file, "%20.15E\t%20.15E\n", x_intersection, y_intersection);
                fclose(file);

            }

        }

        gsl_vector_free(map_gap_x);
        gsl_vector_free(map_gap_y);
        gsl_vector_free(map_dens_x);
        gsl_vector_free(map_dens_y);

        SetFilePath("tests/maps/");
        FILE * log_file = OpenFile("run.log");
        
        fprintf(log_file,
                "Calculates zeroed gap and barionic densities equations so we can see both\n"
                "and have an insight of what's going on.\n");
        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Mass and Renormalized Chemical Potential for Finite Temperature

    // Prints mass and renormalized chemical potential calculation as function
    // of barionic density
    if (true)
    {
        printf("\tMass and Renormalized Chemical Potential for Finite Temperature\n");
        SetFilePath("tests/mass-renorm-chem-pot/data/");

        SetParametersSet("BuballaR_2");

        int n_pts = 100;

        double min_barionic_density = 0.1;
        double max_barionic_denstiy = 0.3;
        double step = (max_barionic_denstiy - min_barionic_density) / ((double)(n_pts - 1));

        double temperature[4] = {10.0, 15.0, 20.0, 25.0};

        gsl_vector * dens_vector = gsl_vector_alloc(n_pts);
        gsl_vector * mass_vector = gsl_vector_alloc(n_pts);
        gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(n_pts);

        for (int i = 0; i < 4; i++){

            double dens = min_barionic_density;
            parameters.finite_temperature.temperature = temperature[i];

            double mass;
            double renormalized_chemical_potential;

            for (int j = 0; j < n_pts; j++){

                SolveMultiRoots(dens,
                               &mass,
                               &renormalized_chemical_potential);

                gsl_vector_set(dens_vector, j, dens);
                gsl_vector_set(mass_vector, j, mass);
                gsl_vector_set(renormalized_chemical_potential_vector, j, renormalized_chemical_potential);

                dens += step;

            }

            char filename[256];
            sprintf(filename, "mass_and_renorm_chem_pot_%d.dat", i);

            WriteVectorsToFile(filename,
                               "# barionic density, mass, renormalized chemical potential\n",
                               3,
                               dens_vector,
                               mass_vector,
                               renormalized_chemical_potential_vector);
        }

        gsl_vector_free(dens_vector);
        gsl_vector_free(mass_vector);
        gsl_vector_free(renormalized_chemical_potential_vector);

        SetFilePath("tests/mass-renorm-chem-pot/");
        FILE * log_file = OpenFile("run.log");
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of mass and renormalized chemical potential as function of barionic density.\n");
        
        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }
    
#pragma mark Entropy calculation methods (transient)
    if (true)
    {
        printf("\tEntropy calculation methods\n");
        SetFilePath("tests/transient/data/");

        SetParametersSet("BuballaR_2");

        double mass[6] = {50.0, 100.0, 200.0, 300.0, 400.0, 500.0};
        double renormalized_chemical_potential[6] = {50.0, 100.0, 200.0, 300.0, 400.0, 500.0};
        double temperature[6] = {1.0, 5.0, 7.0, 10.0, 20.0, 35.0};

        int n_pts = 1000;
        double min_momentum = 0.0;
        double max_momentum = parameters.model.cutoff;
        double mom_step = (max_momentum - min_momentum) / (double)(n_pts - 1);

        // FIXME: Remove other methods of calculation once I'm sure the one adopted works fine
        gsl_vector * momentum_vector = gsl_vector_alloc(n_pts);
        gsl_vector * entropy_integrand_vector = gsl_vector_alloc(n_pts);
        gsl_vector * entropy_integrand_vector_deriv = gsl_vector_alloc(n_pts);
        gsl_vector * entropy_integrand_vector_art = gsl_vector_alloc(n_pts);

        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){
                for (int k = 0; k < 6; k++){
                    entropy_integrand_parameters par;
                    par.mass = mass[i];
                    par.renormalized_chemical_potential = renormalized_chemical_potential[j];
                    par.temperature = temperature[k];

                    double p = 0;
                    for (int l = 0; l < n_pts; l++){

                        double entropy_integrand = EntropyIntegrand(p, &par);
                        double entropy_integrand_deriv = EntropyIntegrandFromDerivative(p, &par);
                        double entropy_integrand_art = EntropyIntegrandArt(p, &par);

                        gsl_vector_set(momentum_vector, l, p);
                        gsl_vector_set(entropy_integrand_vector, l, entropy_integrand);
                        gsl_vector_set(entropy_integrand_vector_deriv, l, entropy_integrand_deriv);
                        gsl_vector_set(entropy_integrand_vector_art, l, entropy_integrand_art);

                        p += mom_step;
                    }

                    char filename[256];
                    sprintf(filename, "entropy_integrand_%d_%d_%d.dat", i, j, k);
                    
                    WriteVectorsToFile(filename,
                                       "# momentum, entropy integrand\n",
                                       2,
                                       momentum_vector,
                                       entropy_integrand_vector);

                    sprintf(filename, "entropy_integrand_deriv_%d_%d_%d.dat", i, j, k);

                    WriteVectorsToFile(filename,
                                       "# momentum, entropy integrand\n",
                                       2,
                                       momentum_vector,
                                       entropy_integrand_vector_deriv);

                    sprintf(filename, "entropy_integrand_art_%d_%d_%d.dat", i, j, k);

                    WriteVectorsToFile(filename,
                                       "# momentum, entropy integrand\n",
                                       2,
                                       momentum_vector,
                                       entropy_integrand_vector_art);

                }
            }
        }

        gsl_vector_free(momentum_vector);
        gsl_vector_free(entropy_integrand_vector);
        gsl_vector_free(entropy_integrand_vector_deriv);
        gsl_vector_free(entropy_integrand_vector_art);

        SetFilePath("tests/transient/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file, "Verifies the best routine to calculate entropy.\n");
        fprintf(log_file,
                "mass = {50.0, 100.0, 200.0, 300.0, 400.0, 500.0};\n"
                "renormalized chemical potential = {50.0, 100.0, 200.0, 300.0, 400.0, 500.0};\n"
                "temperature = {1.0, 5.0, 7.0, 10.0, 20.0, 35.0};\n");
        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Entropy
    // The entropy only depends on the kinectic term, so it does
    // have a 'closed' form. Here we use it to calculate the entropy
    // for a few parameters values just to see if we get a well
    // behaved solution.
    if (true)
    {
        printf("\tEntropy\n");
        SetFilePath("tests/entropy/data/");

        SetParametersSet("BuballaR_2");

        int n_pts = 1000;

        double temperature[4] = {1.0, 5.0, 10.0, 15.0};
        double renormalized_chemical_potential[4] = {50.0, 100.0, 200.0, 350.0};

        double mass_min = 0.0;
        double mass_max = 400.0;
        double mass_step = (mass_max - mass_min) / (double)(n_pts - 1);

        gsl_vector * mass_vector = gsl_vector_alloc(n_pts);
        gsl_vector * entropy_vector = gsl_vector_alloc(n_pts);

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){

                double mass = mass_min;
                for (int k = 0; k < n_pts; k++){
                    double entropy = Entropy(mass, temperature[i], renormalized_chemical_potential[j]);

                    gsl_vector_set(mass_vector, k, mass);
                    gsl_vector_set(entropy_vector, k, entropy);

                    mass += mass_step;
                }

                char filename[256];
                sprintf(filename, "entropy_%d_%d.dat", i, j);

                WriteVectorsToFile(filename,
                                   "# mass, entropy\n",
                                   2,
                                   mass_vector,
                                   entropy_vector);
            }
        }

        gsl_vector_free(mass_vector);
        gsl_vector_free(entropy_vector);

        SetFilePath("tests/entropy/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
                "The entropy only depends on the kinectic term, so it does\n"
                "have a 'closed' form. Here we use it to calculate the entropy\n"
                "for a few parameters values just to see if we get a well\n"
                "behaved solution.\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 2.7 from Buballa Physics Reports
    if (true)
    {
        printf("\tReproduce Fig. 2.7 from Buballa Physics Reports\n");
        SetFilePath("tests/Buballa-Fig2.7/data/");

        SetParametersSet("BuballaR_2");

        int n_pts = 1000;
        double temperature_min = 0.0;
        double temperature_max = 300.0;
        double temperature_step = (temperature_max - temperature_min) / (double)(n_pts - 1);

        gsl_vector * temperature_vector = gsl_vector_alloc(n_pts);
        gsl_vector * mass_vector = gsl_vector_alloc(n_pts);

        double temperature = temperature_min;
        for (int i = 0; i < n_pts; i++){

            parameters.finite_temperature.temperature = temperature;

            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            F.function = &ZeroedGapEquationForFiniteTemperatureTest;

            double mass;
            int status = UnidimensionalRootFinder(&F,
                                                  0.0,
                                                  500.0,
                                                  1.0E-7,
                                                  1.0E-7,
                                                  1000,
                                                  &mass);
            if (status == -1)
                mass = 0;

            gsl_vector_set(temperature_vector, i, temperature);
            gsl_vector_set(mass_vector, i, mass);

            temperature += temperature_step;
        }

        WriteVectorsToFile("mass_temperature.dat",
                           "# Reproduce Fig. 2.7 from Buballa, Physics Reports bare_mass not zero\n"
                           "# temperature, mass\n",
                           2,
                           temperature_vector,
                           mass_vector);

        // Repeat for bare_mass = 0
        parameters.model.bare_mass = 0;

        temperature = temperature_min;
        for (int i = 0; i < n_pts; i++){

            // FIXME: change the program to allow temperature as a variable
            // instead of as a parameter
            parameters.finite_temperature.temperature = temperature;

            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            F.function = &ZeroedGapEquationForFiniteTemperatureTest;

            // This case will not have solutions beyond T = 220 MeV
            double mass = 0;

            if (temperature < 220.0){
                int status = UnidimensionalRootFinder(&F,
                                                      0.1,
                                                      500.0,
                                                      1.0E-7,
                                                      1.0E-7,
                                                      6000,
                                                      &mass);
                if (status == -1)
                    mass = 0;
            }

            gsl_vector_set(temperature_vector, i, temperature);
            gsl_vector_set(mass_vector, i, mass);

            temperature += temperature_step;
        }

        WriteVectorsToFile("mass_temperature_bare_mass_zero.dat",
                           "# Reproduce Fig. 2.7 from Buballa, Physics Reports bare_mass = 0\n"
                           "# temperature, mass\n",
                           2,
                           temperature_vector,
                           mass_vector);

        gsl_vector_free(temperature_vector);
        gsl_vector_free(mass_vector);

        SetFilePath("tests/Buballa-Fig2.7/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file, "Reproduce Fig. 2.7 from Buballa Physics Reports\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Test root finding for renormalized chemical potential
    // After we reach the zero mass case, we use a onedimensional root finding
    // to determine the renormalized chemical potential. However, for some
    // reason the results seems to reach a maximum value. After this maximum,
    // the root finding will fail. One possibility is that we reach a maximum
    // value for the integral.
    if (true)
    {
        printf("\tTest root finding for renormalized chemical potential\n");
        SetFilePath("tests/renorm-chem-pot-transient/data/");

        SetParametersSet("Buballa_1");

        int num_pts = 1000;

        double mass = 0;
        double barionic_density[12] = {0.1, 0.5, 0.7, 1.0, 1.4, 2.0, 2.1, 2.2, 2.3, 2.35, 2.4, 2.5};

        // Range of chemical potential which will be tested
        double renorm_chem_pot_min = 0.0;
        double renorm_chem_pot_max = 3000;
        double step = (renorm_chem_pot_max - renorm_chem_pot_min) / (double)(num_pts - 1);

        gsl_vector * renorm_chem_pot_vector = gsl_vector_alloc(num_pts);
        gsl_vector * zeroed_dens_eq_vector = gsl_vector_alloc(num_pts);

        for (int i = 0; i < 12; i++){

            double mu_r = renorm_chem_pot_min;

            for (int j = 0; j < num_pts; j++){
                gsl_vector_set(renorm_chem_pot_vector, j, mu_r);

                double zeroed_dens_eq = ZeroedBarionicDensityEquationForFiniteTemperature(mass,
                                                                                          mu_r,
                                                                                          &(barionic_density[i]));
                gsl_vector_set(zeroed_dens_eq_vector, j, zeroed_dens_eq);

                mu_r += step;
            }

            char filename[256];
            sprintf(filename, "zeroed_dens_eq_bar_dens_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# renormalized chemical potential, zeroed density equation \n",
                               2,
                               renorm_chem_pot_vector,
                               zeroed_dens_eq_vector);
        }
        
        SetFilePath(NULL);
    }
}

double ZeroedGapEquationForFiniteTemperatureTest(double mass, void * p)
{
    // As the barionic_density is zero, so are the renormalized_chemical_potential
    // and the chemical_potential
    return ZeroedGapEquationForFiniteTemperature(mass, 0.0, NULL);
}
