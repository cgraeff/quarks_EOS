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
* Run tests with `make tests` (it is a shortcut to `make ARGS="-t"` with the default parameterization);
* Plot tests with `make tgraph`;

## Code structure

TBD
