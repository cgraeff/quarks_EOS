//
//  CommandlineOptions.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-20.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef CommandlineOptions_h
#define CommandlineOptions_h

#include <stdbool.h>

typedef struct _options{
	// List options and flags that will be acessible
	// during the execution
	bool verbose;
  	bool skip_tests;
  	bool list_available_parameterizations;
	char * parameterization;
} Options;

extern Options options;

int CommandlineOptionsParse(int argc, char * argv[]);
void CommandlineOptionsPrintUsage(char * prog_name);

#endif /* CommandlineOptions.h */
