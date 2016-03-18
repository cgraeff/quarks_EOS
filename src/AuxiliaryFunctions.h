//
//  AuxiliaryFunctions.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef AuxiliaryFunctions_h
#define AuxiliaryFunctions_h

#include <gsl/gsl_vector.h>

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...);
int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...);

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator);

int WriteZeroedGapEquation(char * filename, double minimum_mass, double maximum_mass, int points_number, double fermi_momentum);


#endif /* AuxiliaryFunctions_h */
