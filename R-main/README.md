The files included in folder `/rubias/R-main/` allow the reproduction of
all data, statistics, and graphs presented in the publication
\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_. They also include numerous graphs for
further visualization of our model validation data, which were excluded
from the publication for brevity.

`cross-validation.R` includes the first two methods of validation listed
in the publication: reproduction of the Leave-One-Out Cross-Validation
in Hasselman *et al.* (2015) with and without parametric bootstrap
correction, and the application of Monte Carlo Cross-Validation to the
same dataset. `coalescent_sim.R` documents the third validation method,
the creation of simulated genotypes from the coalescent, followed by
Leave-One-Out Cross-Validation.

Sourcing these independent files will create (among other files from
intermediate steps) two data.frames for each validation method. One,
`( )_rho_data`, contains true and estimated reporting unit proportions
(rho) created during the corresponding cross-validation. The second,
`( )_rho_dev` contains summary statistics describing the deviation of
each estimation method from its corresponding true rho value. The
validation is identified by the prefix, with `coal`, `Hass`, and `mc`
representing coalescent, Leave-One-Out, and Monte Carlo
cross-validations, respectively. For quick access, these output data are
stored within the main `rubias` folder as `cjfas_data.RData`.

Sourcing will also create a multitude of rough graphs which visualize
this data. `.lg` graphs plot the true vs. estimated rho value *a la*
Figure 3, while `.mse`, `.m.bias`, and `.prop.bias` show mean squared
error, mean bias, and mean proportional bias for each reporting unit,
with and without parametric bootstrapping. Prefixes `c`, `h`, and `mc`
represent coalescent, Leave-One-Out, and Monte Carlo cross-validations.

The remaining two files of the `R-main` folder use the data created by
the previous two; however, they begin by loading `cjfas_data.RData`, and
so can be run independent of the original simulation scripts.
`proportion_sig_test` includes the statistics and rough figures
documenting the relationship between the residual for any given
simulation, $\\tilde{\\rho}\_r - \\rho^\\mathrm{sim}\_r$, and the
structural metric $\\frac{N\_C}{P} - \\rho^\\mathrm{sim}\_r$.
`cjfas_graphs.R` contains the exact scripts used to generate the
[ggplot2](http://ggplot2.tidyverse.org/) figures found in
\_\_\_\_\_\_\_\_\_\_.
