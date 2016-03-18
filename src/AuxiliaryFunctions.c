//
//  AuxiliaryFunctions.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "AuxiliaryFunctions.h"
#include "EOS.h"

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...)
{
	FILE * output = fopen(filename, "w");
    
    if (NULL == output) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
	
	fprintf(output, "%s", header);
	
	va_list arg_list;
	va_start(arg_list, vectors_count);
	
	gsl_vector * vectors[vectors_count];
	
	for (int i = 0; i < vectors_count; i++){
		gsl_vector * v = va_arg(arg_list, gsl_vector *);
		vectors[i] = v;
	}
	
	va_end(arg_list);
	
	if (vectors_count > 1)
		for (int i = 0; i < vectors_count - 1; i++)
			if ((vectors[i])->size != (vectors[i + 1])->size) {
				
				printf("ERROR: Vectors have different sizes.\n");
				exit(EXIT_FAILURE);
			}
	
	for (int i = 0; i < (vectors[0])->size; i++) {
		for (int j = 0; j < vectors_count; j++) {
			
			double x = gsl_vector_get(vectors[j], i);
			
			fprintf(output, "%20.15E", x);
			
			if (j != vectors_count - 1)
				fprintf(output, "\t");
		}
		
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	return 0;
}

int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...)
{
	FILE * output = fopen(filename, "w");
    
    if (NULL == output) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
	
	fprintf(output, "%s", header);
	
	va_list arg_list;
	va_start(arg_list, vectors_count);
	
	gsl_vector * vectors[vectors_count];
	
	for (int i = 0; i < vectors_count; i++) {
		vectors[i] = va_arg(arg_list, gsl_vector *);
	}
	
	va_end(arg_list);
	
	if (vectors_count > 1)
		for (int i = 0; i < vectors_count - 1; i++)
			if ((vectors[i])->size != (vectors[i + 1])->size) {
				
				printf("ERROR: Vectors have different sizes.\n");
				exit(EXIT_FAILURE);
			}
	
	for (int i = 0; i < (vectors[0])->size; i++) {
		
		fprintf(output, "%d\t", i);
		
		for (int j = 0; j < vectors_count; j++) {
			
			fprintf(output, "%20.15E", gsl_vector_get(vectors[j], i));
			
			if (j != vectors_count - 1)
				fprintf(output, "\t");
		}
		
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	return 0;
}

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator)
{
	if (numerator->size != denominator->size) {
		printf("The vectors must have the same size!\n");
		exit(EXIT_FAILURE);
	}
	
	gsl_vector * v = gsl_vector_alloc(numerator->size);
	
	for (int i = 0; i < numerator->size; i++){
		double value = gsl_vector_get(numerator, i) / gsl_vector_get(denominator, i);
		gsl_vector_set(v, i, value);
	}

	return v;
}

int WriteZeroedGapEquation(char * filename, double minimum_mass, double maximum_mass, int points_number, double fermi_momentum){
    
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
    
    while (m < points_number) {
        fprintf(f, "%20.15E\t%20.15E\n", m, ZeroedGapEquation(m, &input));
        m += step;
    }
    
    return 0;
}
