* `WriteIndexedVectorsToFile`:
   * Verify if the function works for just one vector;

* `WriteVectorsToFile`:
   * Verify if the function works for just one vector;

* `UnidimensionalRootFinder`:
    * Verify if the values of the function at `lower_bound` and `upper_bound` have different signs.
    This may be important if we want to run for densities past the chiral restoration;

* `ParametersSetup`:
    * Choose better values for `parameters.vac_mass_det_lower_bound` and `parameters.vac_mass_det_upper_bound`;

* `EnergyDensity`:
    * The chemical potential to be used in this equation is not the renormalized one.
    Since `G_V` is taken to be zero, this should not be a problem for now, but may be in the future,
    when `G_V` is not zero. (`mu_R = mu - 2 * G_V * barionic_density`

* Fix thermodynamic potential calculation

* Change the program to allow temperature as a variable instead of as a parameter

* The zero temperature calculation is completely separated from the finite temperature case.
The ideal way would be to each calculation have a special case for zero temperature.

* Calculate scalar density in finite temperature case. Make plot. Also do that for
Fermi momentum (the integral which gives origin to it; If it makes sense)
