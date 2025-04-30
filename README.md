# G336.01 scripts

These scripts are for processing and plotting the results of the combined ALMA 12m C8 and C9 data.

## Dependencies

The scripts have been run with the following packages:

- CASA (v6.6.3) (with the exception of the calibrated data, which uses the version specified by the weblog)
- [tile_plotter](https://github.com/folguinch/tile_plotter) (v0.2.0)
- [line_little_helper](https://github.com/folguinch/line_little_helper) (v0.2.0)
- [astro_source](https://github.com/folguinch/astro_source) (v0.3.1)
- [GoContinuum](https://github.com/folguinch/GoContinuum) (v3.0.0)
- [velocity_tools](https://github.com/jpinedaf/velocity_tools/tree/v1.1) (v1.1)
- [yclean]() (v2.3.0)

Other dependences are usually installed by the codes above, and include `numpy` (v1.26.4), `astropy` (v6.0.1),
`spectral-cube` (v0.6.5) and `matplotlib` (v3.9.2)

## Scripts for published results

The following scripts were used to produce the results in Olguin et al. (2025):

- `G336_c8c9_workfkow.ipynb`: Steps for processing and extracting the final visibilities.
- `G336.018-0.827_selfcal_c8c9.py`: Steps followed for self-calibration. Results for each step are listed in `G336.018-0.827_selfcal_c8c9.md`.
- `G336.018-0.827_clean_automasking_briggs.py`: Generates the continuum images (including pbcor).
- `G336.018-0.827_continuum_results.py`: Fit Gaussian to the central source and measure streamer fluxes.
- `G336.01-0.82_streamline_model.py`: Compute streamer models.
- `G336.01-0.82_feria_model_update.py`: Compute FERIA model update.
- `G336.01-0.82_feria_model_fix_header.py`: Fix flux unit in model pv map header and get peak flux (to use in plots).
- `continuum_peak_distances.py`: Determine the distance between the central source and streamers continuum peaks.
- `source_pipeline.py`: Process the data cubes for CH3OH and SO to obtain moment and pv maps.
- `keplerian_pvmap_northern.py`: Generates and plots the pv map for the blue-shifted streamer.
- `plots.py`: Generate most of the paper plots.

To produce the contsub and line-free continuum vibilities we used `GoContinuum` with the configuration file `config/reduction/G336.018-0.827.cfg`.

