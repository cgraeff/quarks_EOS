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
#include "Loop.h"

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
        parameters.finite_temperature.temperature = options.temp;
    
    if (parameters.finite_temperature.temperature == 0){
        SolveZeroTemperatureEOS();
    }
    else if (parameters.finite_temperature.temperature > 0){
        SolveFiniteTemperatureEOS();
    }
    else{
        printf("Values of temperature must be non-negative.\n");
        printf("(%f was provided).\n",
			   parameters.finite_temperature.temperature);
        exit(EXIT_FAILURE);
    }

    return 0;
}


