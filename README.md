# rads

The Radar Altimeter Database System (RADS) was developed by the Delft Institute for Earth-Oriented Space Research, the NOAA Laboratory for Satellite Altimetry, and EUMETSAT. Apart from actual altimeter data, RADS provides here a suite of applications and subroutines that simplify the reading, editing and handling of data from various radar altimeters. Although the actual content and layout of the underlying data products do not have to be identical for all altimeters, the user interface is. Also, the data base is easily expandable with additional data and/or additional corrections without impact to the user interface, or to the software in general. In fact, only in very few cases the software will need to be adjusted and recompiled, in even fewer cases adjustments to the actual tools will be required.

## Development
The RADS database and code has gone through various generations. NetCDF datasets were introduced in version 3, while the software more or less stayed the same as in version 2. This new version is a complete rewrite of the code, making it much easier to configure, handle, and expand, and for the first time taking all the advantages of the underlying NetCDF data and the linear algebra provided by Fortran 90.

## Documentation
There are two manuals that cover the use of the RADS software and the contents of the RADS data base:
* [RADS User Manual](https://github.com/remkos/rads/raw/master/doc/manuals/rads4_user_manual.pdf)
* [RADS Data Manual](https://github.com/remkos/rads/raw/master/doc/manuals/rads4_data_manual.pdf)

## Requirements
The only requirements to compile and run the code are:
* A unix type environment (Linux, Mac OS X, etc.).
* The make command.
* A fortran compiler (preferrably gfortran, but others like f90, f95, xlf90, xlf95 and ifort are known to work as well)
* The NetCDF library (version 4.3 or later) and together with the Fortran 90 interface and module file.
* Optionally, the git program.
* For downloading and synchronising the data base: the rsync program.

## Distribution
Tarballs are released on [GitHub 'releases' page](https://github.com/remkos/rads/releases).

## Bug reports and feature requests
Please submit your bug reports or feature requests, and track existing ones, on the [GitHub 'issues' page](https://github.com/remkos/rads/issues). To add content will need to [sign up to GitHub](https://github.com/join?source=header-home).

## Version History
Following is a history of releases on [GitHub](https://github.com/remkos/rads/releases), newest to oldest.
This does not include explanations of changes to the code that generates the data base.

### v4.5.2 (20 Oct 2023)
* Documentation
  * replaced name of rads2asc4 by rads2asc; same for other *4 executables.
* `rads.xml`:
` * updated `ref_frame_offset*` for JA3 (#180)
  * introduced `ref_frame_offset*` for S6A (#180)
  * prepared for SWOT Calval and Science Orbits with cycle numbers shifted -300 and +300 respectively to make them time ordered.
  * added all necessary information for SWOT nadir (sw)
  * adaptations for S6A baseline F08 and prepared for baseline F09
  * introduce `flag_manoeuvre` properfly (#177)
  * update SSB model info for S3A/
  * increase cycle limit for SRAL
  * fix wrong scale on `ssb_3d`, `ssb_3d_adaptive`, `ssb_adaptive`
  * added `*_nr` variables for 6a
  * updated
* `radsstat`:
  * allow to compute statistics for the difference between two satellites/missions
  * introduced `-c/N`, `--mean-only`, and `--groups` options
  * do not add attribute `coordinates` to variables.
  * do not write out phase in ASCII header (#188)
  * extend info in global attributes
* Increased number of characters allowed per option to 640.
* Updated manuals
* Updated `config.sub` and `config.guess`
* Corrected flag bit 0 for S6A (#192)

### v4.4.0 (27 Apr 2022)
* `rads.xml`:
  * prepared for Jason-3 Interleaved Orbit (Phase B)
  * align lower limit of `range_numval_ku_adaptive` with `range_numval_ku` (for j3)
  * added `mss_dtu21` (#174)
  * fix formatting of flag mask of 6a, j2, j3
  * add `*_adaptive` variables for j3 (#173)
  * update `ref_frame_offset` for j3 (#178)
  * correct `standard_name` for `mean_wave_period`
  * remove `ssb_tran2012` from S6A (it was never included in the data files)
  * change limits on sig0 for j3 as result of no longer biasing sig0
  * add `mss_comb15`
  * update limits for `range_numval_ku` (6a), `range_numval_ku_plrm` (3a 3b c2)
  * Update attributes of `ssb_cls`, `ssb_cls_c`, `ssb_mle3` for JA3 and S6A
  * Removed `ssb_hyb` from configuration of c2
  * including SARAL GDR-F wave model
* `radsstat`:
  * improve matching time stamps
  * fix problems (#171)
  * introduced `--no-stddev`
  * introduce check on collinearity as well as `--force` and `--dt` options
  * can now also do differences (along collinear tracks)
  * updated documentation
* `radscolin`:
  * added `--eqtime` option
  * option `--force`: also apply to rogue equator crossing times
  * use inclination to check if missions are collinear
* `radsxolist`: activated `--dual-asc` and `--dual-des` options
* Add MFWAM fields to S3A, S3B, S6A
* Updated manuals
* Updated `config.sub` and `config.guess`

### v4.3.7 (10 Mar 2021)
* Switch CryoSat-2 (c2) data over from NOAA-generated IGDR to NOP/IOP/GOP products
* Updated `config.sub` and `config.guess`
* `rads.xml` prepared for Sentinel-6A (6a)
* Switch SARAL (sa) data over to GDR-F; Change default ssb for SARAL to `ssb_tran2019` (which is what GDR-F does)
* Switch Jason data production over to GDR-F

### v4.3.6 (15 Aug 2019)
* `rads.xml`: Updated `ref_frame_offset` for JA3 (raise by 1 mm)
* `rads.xml`: Added `ref_frame_offset` for S3A/S3B PLRM
* `rads.xml`: Updated JA2 parameters since coming back on-line on 22-May-2019
* `rads.xml`: Provide and use GDR-F orbits for CryoSat-2 (alt_gdrf)
* `rads.xml`: Specified various subphases for sa/b in order to get correct prediction of equator times and longitudes. This solved issue #152.
* Added `conf/ntc_only.xml` and set latency to 2 (NTC) by default.
* Bug fix: Simple modification of syntax for option parsing check to avoid segfault when compiled with default intel-fc/12.1.9.293.
* New option introduced in `radsstat` that echoes to stdout fullpath to each pass file checked for data.
* Using `-L|--limits` on an alias now sets limits of all aliased variables, instead of only the first one.
* Documentation updates.
* Bug fix: Update equator prediction for longitude (not NaN) and using different phases.
* Added command `radspassesindex` as in RADS3.

### v4.3.5 (2 May 2019)
* `rads.xml`: Added `ref_frame_offset_plrm`.
* `rads.xml`: Properly use `wet_tropo_rad_plrm` instead of `wet_tropo_rad` for PLRM data.
* `rads.xml`: Use `alt_gdrf` for S3A and S3B.
* `rads.xml`: Added `latency` variable.
* Properly deal with longitude rollovers in computing means or differences in `radscolin`.
* Remove insistence that 0-dimension variable (constant) is a double.
* Bug fix: invalid values in RPN data notation may not have always worked correctly.
* Restrict the name of the time dimension to the first word in `S%time%info%dataname`.

### v4.3.4 (2 Apr 2019)
* Bug fix in `rads_def_var`.
* Small documentation update.
* Allow multiple mission phases with the same name (needed for Sentinel-3B and Jason-2 after geodetic phase rewind).
* Updated Sentinel-3B mission phases.
* Implemented internal tides.
* Removed support for FES2012 tide model.
* Bug fix in `rads_add_tide`.
* Removed DTU MSS13 from standard models provided.
* Added `topo_srtm15plus` to alias for `topo` ahread of removing `topo_srtm30plus`.

### v4.3.3 (19 Nov 2018)
* Added new optional argument `varid` to `rads_def_var` routine (fixed bug issue #139).

### v4.3.2 (7 Nov 2018)
* New mission definitions for Sentinel-3B.
* Do not store output varid in `rads_varinfo` struct, but determine as needed (fixed issue #126 and #138).
* Update data manual on Sentinel-3A and -3B (fixed issue #131).
* Specify separate `flag_meanings` for each mission (fixed issue #132).

### v4.3.1 (4 Sep 2018)
* New mission phase definition for Jason-2 Phase D.
* New mission definitions for Sentinel-3B.
* Added optional argument "deflate" to `nf90_def_axis`.
* Replaced iqsort by more stable version of quicksort (fixed issue #127).
* No more duplicate tracks in radsxogen when generating single- and dual-satellite xovers (fixed issue #129).
* Documented the --reject-on-nan=all option in radscolin4.

### v4.3.0 (30 May 2018)
* Added support for processing of Sentinel-3B data.
* Added routine `rads_set_phase`.
* Added \<end_time\> specifications on last mission phase of terminated missions.
* Added license file.
* Added `geoid_xgm2016` specification in support of issue #119.
* Added `topo_strm15plus` specification in support of issue #120.
* Fixed issue #121.

### v4.2.11 (16 Mar 2018)
* Updated information on Sentinel-3B orbit and mission phases.
* Provide better "histogram limits" from radsvar.
* Changed format for all time fields in 1985 seconds, so they do not overrun the maximum number of characters.
* Changed the reference frame offset for Sentinel-3A from 27 mm to 2 mm.
* Polishing of manuals

### v4.2.10 (1 Jan 2018)
* Copyright updated to 2018
* Increased number of cycles for CryoSat-2

### v4.2.9 (30 Nov 2017)
* Replaced freeunit() by getlun()
* Added standard_name and axis attributes when writing lat/lon grids
* Added valid_min, valid_max, grid_step attributes to output grids. Removed actual_range.
* Many changes to the developer code to produce the [Nov 2017 updates](https://github.com/remkos/rads/milestone/6?closed=1); see there for more information.
* Changes to configuration file (rads.xml) in accordance with the [Nov 2017 updates](https://github.com/remkos/rads/milestone/6?closed=1).

### v4.2.8 (12 Sep 2017)
* Fixed information on Jason-2 Phase C in `rads.xml`

### v4.2.7 (7 Sep 2017)
* Updates for Jason-2 Phase C ("tango")
* Added definition of qual_alt_rain_ice and qual_rad_rain_ice to Sentinel-3.

### v4.2.6 (27 Aug 2017)
* New features in radscolin: --diff-no-coord and --diff1.
* Numerous updates for the generation of RADS data for Sentinel-3, Jason-2, and Jason-3.
* New variables for modern mean sea surfaces and tide models.
* Updates to user and data manuals.

### v4.2.5 (6 Oct 2016)
* Updates to manuals (e.g., issue #99).
* Numerous updates for the generation of RADS data for Sentinel-3 and Jason-3 (developer tools only).
* Bug fix: nesting of \<if\> tags could go wrong when reading XML files.
* Installation: use nc-config if nf-config is not supplied.
* Adjustments for SARAL Drift Mission.
* Preparation for Jason-2 Interleaved Mission.

### v4.2.4 (6 June 2016)
* Bug fix (issue #97), which prevented radsxogen from running properly. This error was introduced in the master on 1 May 2016 (v4.2.2).
* Enhancements to the user manual.

### v4.2.3 (29 May 2016)
* Bug fix (issue #96), which prevented installation of symbolic links within shells other than bash.

### v4.2.2 (24 May 2016)
* Updated Jason-1 SLCCI orbits to Version 11 for ERS-1, ERS-2, Envisat, Jason-1, Jason-2, and TOPEX/Poseidon.
* Changed reference frame offset for CryoSat-2 Baseline C data from -26 to -34 mm.
* Major additions to the user manual.
* Changes and additions to configuration files for Jason-3 and Sentinel-3A.
* All help output to be sent to standard output (not standard error).
* rads2grd: --fmt replaced by --line-format, variable name in NetCDF now always ${var}_mean and ${var}_stddev, even when only one variable is used.
* Updated global attributes to be written in new NetCDF files.
* --quality_flag replaced by --quality-flag (with backward compatibility).
* Made adaptations for Jason-1 GDR-E (to be released later).
* Change format of time field to f17.6 in advance of 2016-09-09 when extra digit will be added.
* Bug fix (issue #89).
* Bug fix: updated missing building dependencies.

### v4.2.1 (22 Mar 2016)
* Moved directory "include" to "src/include" to avoid clash with default location for "include".
* Updated information on Sentinel-3A.
* Made clearer that -- can be used ahead of FILENAME argument to make sure that the option list is terminated; this addressed issue #85.
* When opening XML file, skip lines until the first line that contains '<'; this fixed issue #87.
* rads2asc: added -c and -p as alternatives to -sc and -sp.
* rads2asc, rads2nc, rads2adr: support --output OUTNAME/ for storing passfiles in directory OUTNAME.
* radscolin: when using -r# (not -rn or -r0) require at least # tracks, not the minimum of # and the available number of tracks.

### v4.2.0 (8 Mar 2016)
* Tuning for processing of Jason-3 data.
* Added \<plot_range\> XML tag, which can be used for generating the right plotting range.
* Added new functionalities to radsvar (i.e. create output that can be included in batch files).
* Replaced 'units = "count"' by 'units = "1"' in line with udunits standards.
* On -S option now allow satellite names with more than 2 characters (e.g. ja3, jason3).
* Introduced the possibility to read RADS data from multiple directories using a tag like \<data branch=".ext"\>.
* Allow "single-line" XML tags with "var" option like \<limits var="sla" sat="j3"\>-0.4 0.4\</limits\>
* Added s3combine to merge and split Sentinel-3 granules into pass files.
* Made some changes to rads_gen_s3, rads_gen_c2, ogdrsplit.

### v4.1.5 (23 Feb 2016)
* Added SARAL GDR-E orbits, which are now the default.
* Further updates for Jason-3.

### v4.1.4 (12 Feb 2016)
* Implement processing of Jason-3 and Sentinel-3 data.
* Implement satellite codes 'j3' (Jason-3), '3a' (Sentinel-3A), '3b' (Sentinel-3B).
* Minor bug fixes.

### v4.1.3 (1 Dec 2015)
* Added Jason-2 GDR-E orbits, which are now the default.

### v4.1.2 (13 Oct 2015)
* Bug fixes for 'radscolin --diff -r0' and for 'rads2grd -c'.
* Add 'configure' to distribution.

### v4.1.1 (27 Aug 2015)
* Bug fixes.
* Improved warnings for erroneous command line options.
* Subroutines documented with robodoc headers.
* Start of RADS User Manual.
* Installing documents with 'make install'.
* Removed 'verbose' argument from rads_init.

### v4.1.0 (19 Aug 2015)
* Numerous updates to the database: removed obsolete variables, updated several, introduced new variables.
* New handling of "long options". For example: can now use either --var=VAR or --var VAR (without = symbol).

### v4.0.3 (10 Aug 2015)
* Improvements to rads_add_surface; new landmask based on GSHHG 2.3.4 introduced.
* Fixed scaling of water_vapor_rad in rads_gen_j1; affected only cycles 1-283.
* Updated data manual with GDR-D/E references.

### v4.0.2 (25 Jun 2015)
* Prepared for GCC and GFortran version 5.
* rads_add_ww3_222 did not properly identify 6-hourly time step for recent grids.
* Improved module dependencies in Makefiles.
* Introduced GDR-E orbit for CryoSat-2.

### v4.0.1 (29 Apr 2015)
* Updated the way the version number is determined, using either the version tag in configure.ac or "git describe --tags" when using a git repository.
* Removed all the SVN $Id$ and $Revision$ tags.
* Updates for future Jason-3 processing.
* Updates for processing of CryoSat-2 Baseline C.

### v4.0.0 (19 Mar 2015)
* Import from GoogleCode into GitHub.

