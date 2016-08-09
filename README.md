# Quarks EOS

This code calculates the Equations of State for quarks following M. Buballa, NJL-model analysis of dense quark matter, Physics Reports 407 (2005) 205-376.

## Requisites

To build and run this code, a C compiler and the GSL (GNU Scientific Library) must be installed on the system. The location of GSL headers and libraries is assumed to be `/usr/local/include` and `/usr/local/lib` on Mac OS. On Linux it should just work.

## Build and run instructions

* The directory structure is important, as they are not created automatically;
* Build with `make`;
* Run with `make run`;
* Arguments may be passed with `make ARGS="-p Set" run`, where `-p Set` stand for example arguments; Set stand for the name of a parameterization identifier (see `Parameters.c` for the available sets);
* Plot results with `make graph`
* Run tests with `make tests` (it is a shortcut to `make ARGS="-t"` with the default parameterization);
* Plot tests with `make graph`;
* Remove product files with `make clean`

## Code structure
