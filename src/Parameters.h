//
//  Parameters.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <stdbool.h>
#include <stdio.h>

#include "ZeroTemperatureEOS.h"
#include "FiniteTemperatureEOS.h"
#include "FermiDiracDistributions.h"
#include "Loop.h"

// Parameterization variables
typedef struct _ModelParameters{
	char * parameters_set_origin;	// Where the set was taken from
   	double G_S;             		// scalar-isoscalar coupling (fm^2)
   	double G_V;             		// vector-isoscalar coupling (fm^2)
   	double cutoff;          		// \Lambda (MeV)
   	double bare_mass;       		// (MeV)
} ModelParameters;

// Parameters for the variable for which we are
// performing the calculation of the EOS
typedef struct _LoopVariableParameters{
	int points_number; // Number of points in which the above range will be divide into
	double min_value;  // (fm^-3)
	double max_value;  // (fm^-3)
} LoopVariableParameters;

typedef struct _parameters
{
    char * parameters_set_identifier;

    // FIXME: this is for a test, take out
    // Parameters used when writing the zeroed gap equation
//    double gap_minimum_mass;
//    double gap_maximum_mass;
//    double gap_points_number;

	LoopVariableParameters loop_variable;
	ModelParameters model;
	ZeroTemperatureParameters zero_temperature;
	VacuumMassDeterminationParameters vacuum_mass_determination;
	RenormalizedChemicalPotentialIntegrationParameters renormalized_chemical_potential_integration;
	FiniteTemperatureParameters finite_temperature;
	FermiDiracIntegralsParameters fermi_dirac_integrals;

} Parameters;

extern Parameters parameters;

void ParametersSetup(void);
void SetParametersSet(char * parameters_set_identifier);

void PrintParametersToFile(FILE * file);

#endif /* Parameters_h */
