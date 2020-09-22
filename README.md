# trfp

Muon g-2 field "trfp" Python module for analyzing trolley/fixed probe correlations.



## Dependencies

This module requires:

1. numpy

2. scipy

3. pandas

4. Simon Corrodi's gm2 package (for interfacing with root), [found here](https://bitbucket.org/corrodis/gm2/src/dev/).

## Overview

Current as of 24 Sept 2019.

1. Data is read directly from the tier 1 ROOT files using S. Corrodi's gm2 package. This can either be done in real time or ahead of time, saving the data into pandas DataFrames. Import using the `root_to_pandas` function in `analysis_helper.py`. This creates an "Interpolation DataFrame."

2. The Interpolation DataFrame is turned into a "Moment DataFrame" using the trolley and fixed probe change of basis matrices. This is in the `calc_moment_df` function in `analysis_helper.py`. At this point, blinds can also be applied from the file `blinds.txt` (see `data_run_analysis.ipynb` for example).

3. The trolley moment dataframes are processed to remove the trolley image from the fixed probe stations. This is done using `trolley_footprint_replacement` in `helper_function_candidates.py`.

4. The baseline sync values for the trolley and fixed probes are calculated using `trolley_run_station_average` in `helper_function_candidates.py`.

5. The fixed probe moment dataframes are transformed into virtual trolley measurement dataframes using the sync values and the Jacobian matrices. This is in `vtm_calc` in`analysis_helper.py`.

6. The resulting vtm dataframe can be averaged in any way desired. Simple time-binning and naive azimuthal averaging are demonstrated in `9day_analysis.ipynb`.

## Sept 2020 Updates

Run 1 final cleanup occurred at the beginning of September 2020 (currently ongoing). The goal is to clean up old files, sort them, and collapse them into relatively few active analysis scripts. These can then be ported into `gm2fieldpy` and `trfp` can be committed and saved for posterity, but leave active development.