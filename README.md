# Quarks EOS

This code calculates the Equations of State for quarks following M. Buballa, NJL-model analysis of dense quark matter, Physics Reports 407 (2005) 205-376.

## Requisites

To build and run this code, a C compiler and the GSL (GNU Scientific Library) must be installed on the system. The location of GSL headers and libraries is assumed to be `/usr/local/include` and `/usr/local/lib` on Mac OS. On Linux it should just work.

## Build and run instructions

Basic build (generate `eos`executable):
* Build with `make` in the top dir;

Running:
* When executed, `eos` will calculate the equations of state with the default parameterization;
* The following options are available:
 * `-p par`: uses `par` parameters set;
 * `-t val`: uses `val` for temperature value. Must be a `double` value; 
 * `-l`: list available parameters set;
 * `-q`: quiet (supress information written to standard out);
 * `-d`: write results using a dir structure;
 * `-a`: run tests. Automatically sets `-d`;
 * `-u`: prints usage;
 * `-h`: also prints usage;

For easier running and testing, the following commands are provided:
* Run with `make run` (implies `-d`);
* Remove product files with `make clean`;
* Arguments may be passed with `make ARGS="-p Set" run`, where `-p Set` stand for example arguments;
* Plot results with `make graph`
* Run tests with `make tests` (it is a shortcut to `make ARGS="-a"` with the default parameterization);
* Plot tests with `make tgraph`;

## Code structure

### Parameters

Parameters sets must be declared in `ParametersSetup()` in `Parameters.c`. Each parameters set is cloned from a reference set with `NewCopyOfParametersSetFromTemplate()` then the variables which differ from the template are set. Each parameters set must have a unique identifier. Finally, the set must be appendend to the list using `AppendParametersSetToList()`.

To choose a parameters set, the function `SetParametersSet()` must be used. This is used before the calculation of the main results and if no parameters set is explicitly requested in the command line, the first set declared is used. Each test in `Tests()` (in `Tests.c`) should declare a set.

Once set using `SetParametersSet()`, use `parameters.a_parameter` to access the parameters values, where `a_parameter` is one of the parameters declared in the struct in `Parameters.h`.

### Commandline options

The commandline options set in `CommandlineOptions.*` are globally available using `options.an_option`.

### Tests

Tests should be declared in the function `RunTests()` in `Tests.c`. This function is executed when the executable is called with `-a` option.