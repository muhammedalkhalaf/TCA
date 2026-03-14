TCA Package - R Documentation Files (*.Rd)
============================================

Generated: 2026-03-14

All required .Rd man page files for the TCA package have been created in this directory.

Files Generated:
================

1. TCA-package.Rd
   - Package-level documentation
   - Overview of main features and functions
   - Citation information for Wegner, Lieb, Smeekes (2025)

2. tca_systems_form.Rd
   - Function: tca_systems_form(Phi0, As, h, order = NULL, Psis = NULL)
   - Builds B and Omega matrices for systems form representation

3. tca_analyze.Rd
   - Function: tca_analyze(from, B, Omega, intermediates, K, h, order, mode = "overlapping", var_names = NULL)
   - Main TCA analysis with support for multiple decomposition modes

4. tca_decompose_binary.Rd
   - Function: tca_decompose_binary(from, B, Omega, var_idx, K, h, order)
   - Binary decomposition into through/not-through components

5. tca_validate_additivity.Rd
   - Function: tca_validate_additivity(from, B, Omega, K, h, order, var_names = NULL, verbose = TRUE)
   - Validates decomposition accuracy

6. tca_from_var.Rd
   - Function: tca_from_var(var_model, from, intermediates, h = 20, order = NULL, mode = "overlapping", identification = "cholesky", Phi0 = NULL)
   - Convenience wrapper for VAR model integration

7. plot_tca.Rd
   - Function: plot_tca(x, target = NULL, type = NULL, title = NULL, colors = NULL)
   - Visualization of channel contributions (ggplot2-based)

8. print.tca_result.Rd
   - S3 Method: print.tca_result(x, target = NULL, ...)
   - Pretty printing of tca_result objects

All .Rd files follow standard R documentation format with:
- Title
- Usage signatures
- Argument descriptions
- Return value documentation
- Detailed descriptions
- Examples (where applicable)
- References to Wegner, Lieb, Smeekes (2025)
