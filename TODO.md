* `WriteIndexedVectorsToFile`:
   * Verify if the function works for just one vector;

* `WriteVectorsToFile`:
   * Verify if the function works for just one vector;

* `UnidimensionalRootFinder`:
    * Verify if the values of the function at `lower_bound` and `upper_bound` have different signs. This may be important if we want to run for densities past the chiral restoration;

* `ParametersSetup`:
    * Choose better values for `parameters.vac_mass_det_lower_bound` and `parameters.vac_mass_det_upper_bound`;

* `EnergyDensity`:
    * The chemical potential to be used in this equation is not the renormalized one. Since `G_V` is taken to be zero, this should not be a problem for now, but may be in the future, when `G_V` is not zero. (`mu_R = mu - 2 * G_V * barionic_density`
