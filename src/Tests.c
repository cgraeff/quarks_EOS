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
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include "EOS.h"
#include "Parameters.h"
#include "CommandlineOptions.h"
#include "AuxiliaryFunctions.h"
#include "Constants.h"

int WriteZeroedGapEquation(char * filename,
                           double minimum_mass,
                           double maximum_mass,
                           int points_number,
                           double fermi_momentum);

int WriteVacuumMassEquation(char * filename,
                            double minimum_mass,
                            double maximum_mass,
                            int points_number);

int WriteZeroedRenormalizedChemicalPotentialEquation(char * filename,
                                                     double minimum_renormalized_chemical_potential,
                                                     double maximum_renormalized_chemical_potential,
                                                     int points_number,
                                                     double chemical_potential,
                                                     double mass);
void RunTests()
{
    FILE * log_file = fopen("tests/tests.log", "w");
    if (log_file == NULL){
        printf("Could not open tests/tests.log file for writting!\n");
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    // Run tests for requested parameterization
    {
        SetParametersSet(options.parameterization);
    
        WriteVacuumMassEquation("tests/data/vacuum_mass_equation.dat", 0.0, 1000.0, 1000);

        WriteZeroedRenormalizedChemicalPotentialEquation("tests/data/zeroed_renorm_chemical_pot_equation.dat",
                                                         0.0,
                                                         1000.0,
                                                         1000,
                                                         0.0,
                                                         0.0);
        
        double barionic_density[4] = {parameters.minimum_density,
                                      parameters.maximum_density * 1.0 / 3.0,
                                      parameters.maximum_density * 2.0 / 3.0,
                                      parameters.maximum_density};
        for (int i = 0; i < 4; i++){
            double fermi_momentum = CONST_HBAR_C
                                    * pow(3.0 * pow(M_PI, 2.0) * barionic_density[i] / NUM_FLAVORS,
                                          1.0 / 3.0);
            char filename[256];
            sprintf(filename, "tests/data/zeroed_gap_equation_%d.dat", i);
            WriteZeroedGapEquation(filename,
                                   0.0,
                                   1000.0,
                                   1000,
                                   fermi_momentum);
        }
        
        
        fprintf(log_file, "The following tests were executed for %s parameterization:\n", parameters.parameters_set_identifier);
        fprintf(log_file, "\tVacuum mass equation was written in data/vacuum_mass_equation.dat\n");
        fprintf(log_file, "\tZeroed renormalized chemical potential equation was written in data/zeroed_renorm_chemical_pot_equation.dat\n");
        fprintf(log_file, "\tZeroed gap equation was written in tests/data/zeroed_gap_equation.dat\n");

    }
    
    fprintf(log_file, "\n");
    
    // Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611 (1996) 393-408
    // (the figure uses parameters of Set II from the article)
    {
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

            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            renorm_chem_pot_equation_input input;
            F.function = &ZeroedRenormalizedChemicalPotentialEquation;
            F.params = &input;

            double m = 0;

            for (int j; j < points_number; j++) {
                
                // Prepare input for ZeroedRenormalizedChemicalPotentialEquation
                input.chemical_potential = chemical_potential[i];
                input.mass = m;

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
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += step;
            }

            char filename[256];
            sprintf(filename, "tests/data/Fig1_Buballa_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential",
                               2,
                               mass_vector,
                               output);
        }
        
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tThermodynamic potential as function of mass "
                "(Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611 (1996) 393-408).\n"
                "\tfiles: data/Fig1_Buballa_*.dat\n");
    }
    
    fprintf(log_file, "\n");

    // Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = 0)
    {
        SetParametersSet("BuballaR_2");
        
        double chemical_potential[4] = {0.0, 300.0, 368.6, 400.0};
        
        double vacuum_mass = VacuumMassDetermination();
        
        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
        
        for (int i = 0; i < 4; i++){
            
            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;
            
            double step = (maximum_mass - minimum_mass) / (points_number - 1);
            
            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);
            
            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            renorm_chem_pot_equation_input input;
            F.function = &ZeroedRenormalizedChemicalPotentialEquation;
            F.params = &input;
            
            double m = 0;
            
            for (int j; j < points_number; j++) {
                
                // Prepare input for ZeroedRenormalizedChemicalPotentialEquation
                input.chemical_potential = chemical_potential[i];
                input.mass = m;
                
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
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += step;
            }
            
            char filename[256];
            sprintf(filename, "tests/data/Fig2.8L_BuballaR_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential",
                               2,
                               mass_vector,
                               output);
        }
        
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tThermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:data/Fig2.8L_BuballaR_*.dat\n");
    }
    
    fprintf(log_file, "\n");
    
    // Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = G_S)
    {
        SetParametersSet("BuballaR_2");
        
        double chemical_potential[4] = {0.0, 430.0, 440.0, 444.3};
        
        double vacuum_mass = VacuumMassDetermination();
        
        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass, 0.0, 0.0, 0.0);
        
        for (int i = 0; i < 4; i++){
            
            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;
            
            double step = (maximum_mass - minimum_mass) / (points_number - 1);
            
            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);
            
            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            renorm_chem_pot_equation_input input;
            F.function = &ZeroedRenormalizedChemicalPotentialEquation;
            F.params = &input;
            
            double m = 0;
            
            for (int j; j < points_number; j++) {
                
                // Prepare input for ZeroedRenormalizedChemicalPotentialEquation
                input.chemical_potential = chemical_potential[i];
                input.mass = m;
                
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
                
                double thermodynamic_potential = ThermodynamicPotential(m,
                                                                        fermi_momentum,
                                                                        chemical_potential[i],
                                                                        renormalized_chemical_potential);
                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output, j, thermodynamic_potential - vacuum_thermodynamic_potential);
                
                m += step;
            }
            
            char filename[256];
            sprintf(filename, "tests/data/Fig2.8R_BuballaR_%d.dat", i);
            
            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential",
                               2,
                               mass_vector,
                               output);
        }
        
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tThermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:data/Fig2.8R_BuballaR_*.dat\n");
    }

    fclose(log_file);
}

int WriteVacuumMassEquation(char * filename,
                            double minimum_mass,
                            double maximum_mass,
                            int points_number){
    
    double m = 0;
    
    double step = (maximum_mass - minimum_mass) / (points_number - 1);
    
    FILE * f = fopen(filename, "w");
    
    if (NULL == f) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    while (m < points_number) {
        fprintf(f, "%20.15E\t%20.15E\n", m, VacuumMassEquation(m, NULL));
        m += step;
    }
    
    fclose(f);
    return 0;
}

int WriteZeroedRenormalizedChemicalPotentialEquation(char * filename,
                                                     double minimum_renormalized_chemical_potential,
                                                     double maximum_renormalized_chemical_potential,
                                                     int points_number,
                                                     double chemical_potential,
                                                     double mass){
    
    double mu = 0;
    
    double step = (maximum_renormalized_chemical_potential - minimum_renormalized_chemical_potential)
    / (points_number - 1);
    
    FILE * f = fopen(filename, "w");
    
    if (NULL == f) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    renorm_chem_pot_equation_input input;
    input.chemical_potential = chemical_potential;
    input.mass = mass;
    
    while (mu < points_number) {
        fprintf(f,
                "%20.15E\t%20.15E\n",
                mu,
                ZeroedRenormalizedChemicalPotentialEquation(mu, &input));
        mu += step;
    }
    
    fclose(f);
    return 0;
}

int WriteZeroedGapEquation(char * filename,
                           double minimum_mass,
                           double maximum_mass,
                           int points_number,
                           double fermi_momentum){
    
    double m = 0;
    
    double step = (maximum_mass - minimum_mass) / (points_number - 1);
    
    FILE * f = fopen(filename, "w");
    
    if (NULL == f) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    gap_equation_input input;
    input.fermi_momentum = fermi_momentum;
    
    while (m < maximum_mass) {
        fprintf(f, "%20.15E\t%20.15E\n", m, ZeroedGapEquation(m, &input));
        m += step;
    }
    
    fclose(f);
    return 0;
}
