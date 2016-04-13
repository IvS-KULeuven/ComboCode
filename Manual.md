# Welcome to the ComboCode User Manual
## Introduction
### What is this manual?
This manual is meant as a guide to running ComboCode and its two numerical codes. This is not an all-inclusive, comprehensive manual for the ComboCode capabilities. However, in-depth, up-to-date documentation for the Python package is available on <a href="https://IvS-KULeuven.github.io/ComboCode">on GitHub</a>, as a collection of doc-strings. The source code is also available there. Together with the in-depth documentation and the cookbooks provided in this manual for extracting information and using the additional modules, you should be able to use the package to its full extent. For additional questions or remarks, please contact <a href="https://github.com/robinlombaert">R. Lombaert</a>. 

For best results, it is recommended to use both radiative-transfer (RT) codes embedded in ComboCode: GASTRoNOoM (line RT) and MCMax (continuum RT), authored and maintained by L. Decin and M. Min, respectively. Both authors have to be contacted for use of the RT codes, also as part of ComboCode, and their contact details can be requested from <a href="https://github.com/robinlombaert">R. Lombaert</a>.

### Is this manual finished?
No! This is very much a work in progress. Any comments, questions or clarifications can be requested from <a href="https://github.com/robinlombaert">R. Lombaert</a>. We will try to keep the manual as up-to-date as possible, but for a while yet, it will not be complete. Bear with us!

## Goals of the ComboCode package
ComboCode is a Python based package designed to work with radiative-transfer codes and the data they are meant to model. 
The radiative transfer is usually calculated for cool stellar winds of evolved stars, such as AGB stars.
The functionality includes:
* <b>Modeling</b>:
    - Currently works with GASTRoNOoM for gas radiative transfer, and with MCMax for dust radiative transfer
    - Allows output of one code to be used as input for another
    - Automatic gas line selection based on available data and line listing
    - Databases for modeling output and easy parameter space searches
    - Interaction with a supercomputer cluster built into the databases
* <b>Data</b>: 
    - Management of data files associated with radio data, SEDs and spectroscopic data
    - Fitting routines for resolved emission lines
    - Statistical analysis for samples and individual sources, and both resolved and unresolved emission lines

## Running ComboCode
In what follows, you will set up your folder structure (much of which is done automatically) first. Then the ComboCode inputfile is described, and, finally, a few simple steps to run a model are given.

### Folder setup
The folder setup for running ComboCode can be split up into three parts: the folders specific to ComboCode, and the two folders for the GASTRoNOoM and MCMax RT codes that must be installed separately. 

Firstly, the folder setup that comes with ComboCode is set up when installing the git repository. Only the usr/ folder must be installed manually, as this folder contains the local user-specific settings. This is described in the <a href="https://github.com/IvS-KULeuven/ComboCode/blob/master/README.md">README</a> document. For completeness, these folders include:
<ol>
<li>cc -- Contains the Python modules.</li>
<li>aux -- Contains auxiliary files that ComboCode requires, and they come as part of the repository. These files are not to be changed by the user.</li>
<li>usr -- Contains files that ComboCode requires but are not part of the repository. The settings are user specific.</li>
<li>usr.dist -- The blueprint for usr/ that is part of the repository. This folder is to be copied to usr/ upon installation (see below).</li></ol>

In the usr/ folder, a file called Path.dat is located. This file manages all other folder locations not installed by git. Take your time choosing the folder locations in this file, and checking whether these folders exist. Typically, you will want a GASTRoNOoM, MCMax, and Data home folder. Most of the other folders can be chosen at will but usually belong in one of these three home folders. The default Path.dat in usr.dist makes suggestions for all of these. Finally, the instructions for the syntax of the folders is explained at the top of the Path.dat file.

Some of the folders contain data and files that must be provided for you by either the authors of MCMax or GASTRoNOoM, or by other users. These include (given by their path keys in Path.dat): 
<ol>
<li>ivsdata -- Data required by the cc.ivs module. They are available <a href="http://ster.kuleuven.be/~robinl/cc/ivsdata.tar.gz"> here</a> for download. These should not be changed by the user.</li>
<li>atm -- Medium-resolution spectra of model stellar atmospheres. These are included as a subfolder of the ivsdata and the path should be adjusted accordingly.</li>
<li>gdata -- Data required by GASTRoNOoM, including collision rates and spectroscopic descriptions of the molecules. These can be changed depending on user-specific needs but should only be done in dialogue with L. Decin. A set of molecules is available <a href="http://ster.kuleuven.be/~robinl/cc/GASTRoNOoM_data.zip">here</a>. Contact R. Lombaert or L. Decin for literature references or additional molecular species.</li>
<li>mobs -- Observation output files defining the settings for ray tracing of the spectra, visibilities, etc. These can be changed by the user at will. Example observation output files are available <a href="http://ster.kuleuven.be/~robinl/cc/Observation_Files.tar.gz">here</a>. Instructions on the use of these files is available at the <a href="https://sites.google.com/site/manualmcmax/project-definition">MCMax home page</a>.</li>
<li>mopac -- The dust opacities used by MCMax. These can be changed by the user at will, and each opacity file to be used in MCMax must be listed in usr/Dust.dat. A limited set of example opacity files used by MCMax are available <a href="http://ster.kuleuven.be/~robinl/cc/opacities.tar.gz">here</a>. The literature references for the dust properties are listed in the .info files. R. Lombaert can be contacted for additional dust species.</li>
<li>dphot -- The output folder for raw photometric data automatically downloaded by ComboCode when asking for data type 'Photometric_IvS' in Sed.dat, see below.</li>
</ol>

### User files in cc/usr
The contents of the cc/usr.dist folder must be copied to cc/usr/. The default settings will work fine, but can be changed depending on the needs of the user.  
<ol>
<li>Data.dat -- Contains wavelength regions for determining spectral RMS in separate bands of unresolved line observations such as PACS and SPIRE.</li>
<li>Dust.dat -- Describes dust species in terms of specific density, molar weight, destruction temperature and coefficients (Kama et al. 2009), and the subfolder/filename of the opacity file. It doesn't matter if .opac or .particle is listed. ComboCode will select the relevant file if available in the same location. List .topac in case temperature-dependent opacities are wanted. The dust name tag (SPECIES_SHORT) is the handle used in the ComboCode inputfile to refer to these dust species.</li>
<li>Indices.dat -- Lists the different sets of collision rates for several molecular species. The files must be correctly linked in the gdata folder discussed above when calculating new GASTRoNOoM models. When only reading existing models from the databases, Indices.dat allows comparing results for different sets of collision rates without having to link to the correct files in gdata. Only applicable for molecular species listed in Molecule.dat with USE_INDICES_DAT=1.</li>
<li>Molecule.dat -- Description of molecular species as used by ComboCode and GASTRoNOoM. Except the NAME_PLOT columns, these should not be changed unless new molecular species are added.</li>
<li>Path.dat -- Lists the folder locations used by ComboCode, see above.</li>
<li>Sed.dat -- The types of SED data to be included in plotting dust spectra. The DATA_TYPES column gives the name tag included in the dsed folder, see Data Management below. The ABS_ERR gives the absolute flux calibration uncertainty for spectroscopic data, and is ignored for photometric data. Photometric_IvS is a fixed data type used for automatically downloading photometry for a given source.</li>
<li>Star.dat -- Stars and their properties (long, lat and vlsr). The STAR_NAME is the name tag used in the data files and the ComboCode inputfile.</li>
<li>Telescope.dat -- Radio, submillimeter and far-infrared telescopes and/or instruments. Lists the telescope size and absolute flux calibration uncertainty for each. The name tag must be included in the filenames in the dradio folder in case of radio observations, see below under Data Management.</li>
</ol>

### The ComboCode inputfile

### How do I run ComboCode?





## Data management

### Resolved Molecular emission (Radio)

### Unresolved molecular emission (Infrared)

### Spectral energy distribution
Filename convention.



## Model management
### Combined dust and gas radiative transfer
ComboCode is an interface that provides access to two numerical RT codes for dust and gas respectively. The way these codes are linked through ComboCode is illustrated in the schematic below. 
![](aux/flow_chart_codes.pdf)

### Reading and using model output

### Plotting model output

### Database management




## Statistical methods

### Measuring goodness-of-fit

### Trend analysis




## Additional modules

### Line profile fitting

### Plotting line lists




