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
Options options = {false, NULL};

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
	
	char * short_options = "p:vuh";
	static struct option long_options[] = {{"parameterization", required_argument, NULL, 'p'},
										   {"verbose", no_argument, NULL, 'v'},
										   {"usage", no_argument, NULL, 'u'},
										   {"help", no_argument, NULL, 'v'},
										   {NULL, 0, NULL, 0}};
	int opt;
	while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1){
	
		// If an option have an argument, it is accessed through 'optarg'
		switch (opt){
			case 'p':
				options.parameterization = optarg;
				break;
			case 'v':
				options.verbose = true;
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
				// getopt_long should print an error
				printf("Use --usage or -u to print a list of accepted options.\n");
				exit(EXIT_FAILURE);
				break;
			default:
				printf("Use --usage or -u to print a list of accepted options.\n");
				exit(EXIT_FAILURE);
		}
	}


	// Print any remaining command line arguments (not options)
	if (optind < argc){
		printf ("%s: invalid arguments -- ", argv[0]);
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
		printf("Use --usage or -u to print a list of accepted options.\n");
		exit(EXIT_FAILURE);
	}
	
	return 0;
}

void CommandlineOptionsPrintUsage(char * prog_name)
{
	printf("Usage: %s [options]\n", prog_name);
	printf("Options:\n"
		   "\t--parameterization, -p: Chooses a parameterization.\n"
		   "\t--verbose, -v: Prints information to show progress.\n"
		   "\t--usage, -u: Prints this message.\n"
		   "\t--help, -h: Same as --usage.\n");
	
	return;
}

