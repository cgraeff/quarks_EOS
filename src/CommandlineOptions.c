//
//	CommandlineOptions.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-20.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "CommandlineOptions.h"

// Default values for options and flags that will be acessible
// during the execution (specified in order of declaration).
Options options = {true, false, false, NULL};

int CommandlineOptionsParse(int argc, char * argv[])
{
	// Short options and long options must be declared in the next two variables.
	// In short_options, each short option is declared in a string by the character that
	// will invoke the option. If the option takes an argument, a colon (:) must be
	// follow the character. When a character contained in this string is found in the
	// options, the number that corresponds to this caracter is returned by getopt_long
	// (that is "a" -> 'a', ...).
	// EXAMPLE:
	// char * short_options = "a:bvuh";
	//
	// In long_options, an option struct is defined with: a long name, a declaration
	// of existence (or not) of argument for the option, an address of an int 
	// variable to be written with value 'val', and the value 'val'.
	// If instead of an address, NULL is given, then when getopt_long is invoked,
	// it will just return 'val' uppon finding the long option. If we just use for 'val'
	// the same character as in short_options, but numerically (with '' instead of ""),
	// the return value will be the same as in short_options, so both cases can be
	// covered by the same case in the switch.
	// EXAMPLE:
	// static struct option long_options[] = {{"a_option", required_argument, NULL, 'a'},
	//									      {"b_flag", no_argument, NULL, 'b'}};
	//
	// BEWARE: if an option takes many arguments, the spaces must be escaped, otherwise
	// the arguments after the first will be misinterpreted as unknown, or unclaimed.
	// This particular implementation will stop if there are any unprocessed arguments.
	
	char * short_options = "p:lquh";

	int opt;
	while ((opt = getopt(argc, argv, short_options)) != -1){

		// If an option have an argument, it is accessed through 'optarg'
		switch (opt){
			case 'p':
				options.parameterization = optarg;
				break;
			case 'l':
		  		options.list_available_parameterizations = true;
		  		break;
			case 'q':
				options.verbose = false;
				break;
			case 'u':
				CommandlineOptionsPrintUsage(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case 'h':
				CommandlineOptionsPrintUsage(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
		  		if (optopt == 'p'){
          			fprintf (stderr,
							 "Option -%c requires an argument. Use -h for help.\n",
							 optopt);
				}else if (isprint (optopt)){
          			fprintf (stderr,
							 "Unknown option `-%c'.  Use -h for help.\n",
							 optopt);
				}else{
          			fprintf (stderr,
							 "Unknown option character `\\x%x'.  Use -h for help.\n",
							 optopt);
				}
				exit(EXIT_FAILURE);
				break;
			default:
				printf("Use -h to print a list of accepted options.\n");
				exit(EXIT_FAILURE);
		}
	}


	// Print any remaining command line arguments (not options)
	// and let the user know that they are invalid
	if (optind < argc){
		printf ("%s: invalid arguments -- ", argv[0]);
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
		printf("Use -h or -u to print a list of accepted options.\n");
		exit(EXIT_FAILURE);
	}
	
	return 0;
}

void CommandlineOptionsPrintUsage(char * prog_name)
{
	printf("Usage: %s [options]\n", prog_name);
	printf("Options:\n"
		   "\t-p: Chooses a builtin parameterization.\n"
		   "\t-l: Lists availeable builtin parameterizations.\n"
		   "\t-f: Chooses a parameterization file.\n"
		   "\t-t: Saves a template parameterization file in current dir.\n"
		   "\t-q: Supress information (may be useful in scripts).\n"
		   "\t-u: Prints this message.\n"
		   "\t-h: Same as -u.\n");
	
	return;
}

