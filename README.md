# VitalRatesTest

This repository includes the code and data for carrying out the analysis in the manuscript titled: Inter-specific Variability in Demographic Processes Affects Abundance-Occupancy Relationships (Sen and Akcakaya, 2022). Here, we applied a recently introduced mark-recapture method, called RD-pop, to 42 bird species and used results from 17 species to investigate macroecological patterns related to abundace and intrinsic growth rate of populations. See https://github.com/bilgecansen/CJS-pop for details on RD-pop.

We fit CJS-pop to mark-recapture data obtained from Mapping Avian Productivity and Survivorship (MAPS) program. Unfortunately, we don't have yet permission to share MAPS data. We share the code that prepares the data and runs the analysis. The script **wrangle_chdata.R** prepares the MAPS data for model fitting. RD-pop model specifications can be found in **jags_script_cjs.R**. **models_cjspop_weather.R** runs CJS-pop using JAGS; while this code is written for the SeaWulf Cluster in Stony Brook University it can easily be adapted to run in other clusters. **summarize_cjspop.R** selects species suitable for macroecological analysis, and **summarize_dem.R** calculates relevant demographic parameters such as intrinsic growth rate and abundance. We share these results in **results_dem.R**. **macro_eco_plots.R** caries out the macroecological analysis (jags code for this segment is in **models_jags.R**) and it can reproduce all the graphs reported in Sen and Akcakaya (2022). **macro_eco_plots.R** requires only **results_dem.R** as data to run.

Şen B, Akçakaya HR (2022) Inter-specific variability in demographic processes affects abundance-occupancy relationships. Oecologia. https://doi.org/10.1007/s00442-021-05085-5
