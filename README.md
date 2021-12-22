# Future world population projections 
combines SSP population projections and RCP urban fractions to grid population for the future 

Author: Niklas Boke Olen 
niklasbokeolen@gmail.com                     



 R - Files need to be run in the following order


1. world_pop_prepare_data.R // some preparation and harmonization of data
1. CMPI6_worldpop_2010_unique.R   /// create the unique dataset for year 2010
2. CMIP6_population_gridding_world_aurora.R  /// main population gridding -- runs in parallell specified by cores parameter
