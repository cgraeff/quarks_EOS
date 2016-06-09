//
//	CommandlineOptions.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-20.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stddef.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>

#include "CommandlineOptions.h"

// Default values for options and flags that will be acessible
// during the execution (specified in order of declaration).
Options options = {true, false, false, NULL};

int CommandlineOptionsParse(int argc, char * argv[])
{
	// Short options must be declared in the next variable.
	// In short_options, each short option is declared in a string by the character that
	// will invoke the option. If the option takes an argument, a colon (:) must be
	// follow the character. When a character contained in this string is found in the
	// options, the number that corresponds to this caracter is returned by getopt()
	// (that is "a" -> 'a', ...).
	// EXAMPLE:
	// char * short_options = "a:bvuh";
	//
	// BEWARE: if an option takes many arguments, the spaces must be escaped, otherwise
	// the arguments after the first will be misinterpreted as unknown, or unclaimed.
	// This particular implementation will stop if there are any unprocessed arguments.
	
	char * short_options = "p:lqsuh";

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
            case 's':
                options.skip_tests = true;
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

