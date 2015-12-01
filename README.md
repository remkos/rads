# rads

The Radar Altimeter Database System (RADS) was developed by the Delft Institute for Earth-Oriented Space Research and the NOAA Laboratory for Satellite Altimetry. Apart from actual altimeter data, RADS provides here a suite of applications and subroutines that simplify the reading, editing and handling of data from various radar altimeters. Although the actual content and layout of the underlying data products do not have to be identical for all altimeters, the user interface is. Also, the data base is easily expandable with additional data and/or additional corrections without impact to the user interface, or to the software in general. In fact, only in very few cases the software will need to be adjusted and recompiled, in even fewer cases adjustments to the actual tools will be required.

## Development
The RADS database and code has gone through various generations. NetCDF datasets were introduced in version 3, while the software more or less stayed the same as in version 2. This new version is a complete rewrite of the code, making it much easier to configure, handle, and expand, and for the first time taking all the advantages of the underlying netCDF data and the linear algebra provided by Fortran 90.

## Requirements
The only requirements to run the code are:
* A unix type environment (Linux, Mac OS X, etc.)
* A fortran compiler (preferrably gfortran, but others are known to work as well)
* The netCDF fortran library (version 4 or later preferred).
* Potentially in the future: CMake to provide the configuration before compilation

## Distribution
The code is currently still under major development. When released, it will be made available as a tarball, as well as through git.

## Version History
Following is a history of releases on GitHub, newest to oldest.

### v4.1.3
* Added Jason-2 GDR-E orbits, which are now the default.

### v4.1.2
* Bug fixes for 'radscolin --diff -r0' and for 'rads2grd -c'.
* Add 'configure' to distribution.

### v4.1.1
* Bug fixes.
* Improved warnings for erroneous command line options.
* Subroutines documented with robodoc headers.
* Start of RADS User Manual.
* Installing documents with 'make install'.
* Removed 'verbose' argument from rads_init.

### v4.1.0
* Numerous updates to the database: removed obsolete variables, updated several, introduced new variables.
* New handling off "long options". For example: can now use either --var=VAR or --var VAR (without = symbol).

### v4.0.3
* Improvements to rads_add_surface; new landmask based on GSHHG 2.3.4 introduced.
* Fixed scaling of water_vapor_rad in rads_gen_j1; affected only cycles 1-283.
* Updated data manual with GDR-D/E references.

### v4.0.2
* Prepared for GCC and GFortran version 5.
* rads_add_ww3_222 did not properly identify 6-hourly time step for recent grids.
* Improved module dependencies in Makefiles.
* Introduced GDR-E orbit for CryoSat-2.

### v4.0.1
* Updated the way the version number is determined, using either the version tag in configure.ac or "git describe --tags" when using a git repository.
* Removed all the SVN $Id$ and $Revision$ tags.
* Updates for future Jason-3 processing.
* Updates for processing of CryoSat-2 Baseline C.

### v4.0.0
* Import from GoogleCode into GitHub.

