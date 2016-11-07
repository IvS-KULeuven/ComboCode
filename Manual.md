# Welcome to the ComboCode User Manual
## 1. Introduction
<a name="what">
### What is this manual?
</a>
This manual is meant as a guide to running ComboCode and its two numerical codes. This is not an all-inclusive, comprehensive manual for the ComboCode capabilities. However, in-depth, up-to-date documentation for the Python package is available on <a href="https://robinlombaert.github.io/ComboCode">on GitHub</a>, as a collection of doc-strings. The source code is also available there. Together with the in-depth documentation and the cookbooks provided in this manual for extracting information and using the additional modules, you should be able to use the package to its full extent. For additional questions or remarks, please <a href="https://github.com/robinlombaert">contact me</a>. 

For best results, it is recommended to use both radiative-transfer (RT) codes embedded in ComboCode: GASTRoNOoM (line RT) and MCMax (continuum RT), authored and maintained by L. Decin and M. Min, respectively. Both authors have to be contacted for use of the RT codes, also as part of ComboCode, and their contact details can be requested <a href="https://github.com/robinlombaert">from me</a>. A current work-in-progress is the implementation of the MCP & ALI codes used at Chalmers University of Technology, Sweden, currently maintained by H. Olofsson & M. Maercker.

This manual has been written by Robin Lombaert and revised by Marie Van de Sande.

<a name="finished">
### Is this manual finished ?
</a>
No! This is very much a work in progress. Any comments, questions or clarifications can be requested <a href="https://github.com/robinlombaert">by contacting me</a>. We will try to keep the manual as up-to-date as possible, but for a while yet, it will not be complete. Bear with us!

### Acknowledgments
ComboCode was originally developed at <b>KU Leuven, Belgium</b>. Currently, development and maintenance continue at <b>Chalmers University of Technology, Sweden</b>.

Developers: 
- <a href="https://github.com/robinlombaert">Robin Lombaert</a>, Chalmers University of Technology, Sweden (main author)
- <a href="https://github.com/MarieVdS">Marie van de Sande</a>, KU Leuven, Belgium (contributing author)

### The numerical codes currently included in ComboCode
* <b> GASTRoNOoM:</b> L. Decin (KU Leuven, Belgium) 
* <b> MCMax:</b> M. Min (UvA, the Netherlands)
* <b> MCP & ALI:</b> H. Olofsson & M. Maercker (Chalmers, Sweden) -- WiP 

### The numerical modules written in ComboCode
* <b> EnergyBalance: </b> R. Lombaert, H. Olofsson & M. Maercker (Chalmers, Sweden)

For the use of the EnergyBalance module, please contact one of the three authors listed here.

### Table of Contents
1. <a href="Manual.md#1-introduction">Introduction</a>
    - <a href="Manual.md#what">What is this manual?</a>
    - <a href="Manual.md#finished">Is this manual finished?</a>
    - <a href="Manual.md#acknowledgments">Acknowledgments</a>
    - <a href="Manual.md#the-numerical-codes-currently-included-in-combocode">The numerical codes currently included in ComboCode</a>
    - <a href="Manual.md#the-numerical-modules-written-in-combocode">The numerical modules written in ComboCode</a>
2. <a href="Manual.md#2-goals">Goals</a>
3. <a href="Manual.md#3-setting-up-combocode">Setting up ComboCode</a>
    - <a href="Manual.md#folder-setup">Folder setup</a>
    - <a href="Manual.md#userfiles">User files in cc/usr</a>
4. <a href="Manual.md#4-running-combocode">Running ComboCode</a>
    - <a href="Manual.md#inputfile">The ComboCode inputfile</a>    
    - <a href="Manual.md#how">How do I run ComboCode?</a>
    - <a href="Manual.md#object">The ComboCode object and its model contents</a>
5. <a href="Manual.md#5-data-management">Data Management</a>
    - <a href="Manual.md#resolved">Resolved molecular emission (radio)</a>    
    - <a href="Manual.md#unresolved">Unresolved molecular emission (infrared)</a>
    - <a href="Manual.md#spectral-energy-distribution">Spectral energy distribution</a>
6. <a href="Manual.md#6-model-management">Model Management</a>
    - <a href="Manual.md#combined-dust-and-gas-radiative-transfer">Combined dust and gas radiative transfer</a>
    - <a href="Manual.md#reading-and-using-model-output">Reading and using model output</a>
        1. <a href="Manual.md#spec">Molecular spectroscopy and radiative-transfer output</a>
        2. <a href="Manual.md#raytrace">Modeled line profiles and ray-tracing output for GASTRoNOoM</a>
        3. <a href="Manual.md#dust">Dust opacities for MCMax</a>
        4. <a href="Manual.md#general">General model input/output for GASTRoNOoM and MCMax</a>
    - <a href="Manual.md#plotting-model-output">Plotting model output</a>
7. <a href="Manual.md#7-database-management">Database Management</a>
    - <a href="Manual.md#types-of-databases">Types of databases</a>
    - <a href="Manual.md#database-locking">Database locking</a>
    - <a href="Manual.md#cleaning-your-databases">Cleaning your databases</a>
8. <a href="Manual.md#8-the-energybalance">The EnergyBalance</a>
    - <a href="Manual.md#setting-up-and-running-the-energybalance">Setting up and running the EnergyBalance</a>
        1. <a href="Manual.md#EBinput">Defining input for the EnergyBalance</a>
        2. <a href="Manual.md#EBiter">Iterating the EnergyBalance</a>
    - <a href="Manual.md#profilers">Profilers</a>
    - <a href="Manual.md#iteration-with-ali">Iteration with ALI</a>
9. <a href="Manual.md#9-statistical-methods">Statistical Methods</a>
    - <a href="Manual.md#measuring-goodness-of-fit">Measuring goodness of fit</a>
    - <a href="Manual.md#trend-analysis">Trend analysis</a>
10. <a href="Manual.md#10-additional-modules">Additional Modules</a>
    - <a href="Manual.md#line-profile-fitting">Line profile fitting</a>
    - <a href="Manual.md#plotting-line-lists">Plotting line lists</a>

## 2. Goals
ComboCode is a Python based package designed to work with radiative-transfer codes and the data they are meant to model. 
The radiative transfer is usually calculated for cool stellar winds of evolved stars, such as AGB stars.
The functionality includes:
* <b>Modeling</b>:
    - Currently works with GASTRoNOoM for gas radiative transfer, and with MCMax for dust radiative transfer
    - Allows output of one code to be used as input for another
    - Automatic gas line selection based on available data and line listing
    - Databases for modeling output and easy parameter space searches
    - Interaction with a supercomputer cluster built into the databases
    - Reading and using molecular spectroscopy, collision rates and level populations
    
* <b> Theory</b>:
    - Independent calculation of the energy balance of a stellar wind
    - Iteration with the ALI numerical radiative-transfer code for consistent molecular level populations
    
* <b>Data</b>: 
    - Management of data files associated with radio data, SEDs and spectroscopic data
    - Fitting routines for resolved emission lines
    - Statistical analysis for samples and individual sources, and both resolved and unresolved emission lines

## 3. Setting up ComboCode
In what follows, you will set up your folder structure first (much of which is done automatically). This concerns the entire ComboCode package, including the ComboCode and EnergyBalance classes.

### Folder setup
The folder setup for running ComboCode can be split up into three parts: the folders specific to ComboCode, and the two folders for the GASTRoNOoM and MCMax RT codes that must be installed separately. 

Firstly, the folder setup that comes with ComboCode is set up when installing the git repository. Only the usr/ folder must be installed manually, as this folder contains the local user-specific settings. This is described in the <a href="https://github.com/robinlombaert/ComboCode/blob/master/README.md">README</a> document. For completeness, these folders include:

<ol>
<li>cc -- Contains the Python modules.</li>
<li>aux -- Contains auxiliary files that ComboCode requires, and they come as part of the repository. These files are not to be changed by the user.</li>
<li>usr -- Contains files that ComboCode requires but are not part of the repository. The settings are user specific.</li>
<li>usr.dist -- The blueprint for usr/ that is part of the repository. This folder is to be copied to usr/ upon installation (see below).</li></ol>

In the usr/ folder, a file called Path.dat is located. This file manages all other folder locations not installed by git. Take your time choosing the folder locations in this file, and checking whether these folders exist. Typically, you will want a GASTRoNOoM, MCMax, and Data home folder (given by their small-case equivalents in Path.dat). Most of the other folders can be chosen at will but usually belong in one of these three home folders. The default Path.dat in usr.dist makes suggestions for all of these. Finally, the instructions for the syntax of the folders is explained at the top of the Path.dat file.

Some of the folders contain data and files that must be provided for you by either the authors of MCMax or GASTRoNOoM, or by other users. These include (given by their path keys in Path.dat):  

<ol>
<li>ivsdata -- Data required by the cc.ivs module. They are available <a href="http://ster.kuleuven.be/~robinl/cc/ivsdata.tar.gz"> here</a> for download (but you likely already downloaded these when installing the requirements of ComboCode). These should not be changed by the user.</li>
<li>atm -- Medium-resolution spectra of model stellar atmospheres. These are included as a subfolder of the ivsdata and the path should be adjusted accordingly.</li>
<li>gdata -- Data required by GASTRoNOoM, including collision rates and spectroscopic descriptions of the molecules, and the default dust opacity profile. These can be changed depending on user-specific needs but should only be done in dialogue with L. Decin. A set of molecules is available <a href="http://ster.kuleuven.be/~robinl/cc/GASTRoNOoM_data.zip">here</a>. Contact R. Lombaert or L. Decin for literature references or additional molecular species.</li>
<li>mobs -- Observation output files defining the settings for ray tracing of the spectra, visibilities, etc. These can be changed by the user at will. Example observation output files are available <a href="http://ster.kuleuven.be/~robinl/cc/Observation_Files.tar.gz">here</a>. Instructions on the use of these files is available at the <a href="https://sites.google.com/site/manualmcmax/project-definition">MCMax home page</a>.</li>
<li>mopac -- The dust opacities used by MCMax. These can be changed by the user at will, and each opacity file to be used in MCMax must be listed in usr/Dust.dat. A limited set of example opacity files used by MCMax are available <a href="http://ster.kuleuven.be/~robinl/cc/opacities.tar.gz">here</a>. The literature references for the dust properties are listed in the .info files. R. Lombaert can be contacted for additional dust species.</li>
<li>dphot -- The output folder for raw photometric data automatically downloaded by ComboCode when asking for data type 'Photometric_IvS' in Sed.dat, see below.</li>
</ol>

<a name="userfiles">
### User files in cc/usr
</a>
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

## 4. Running ComboCode
In this section, the ComboCode inputfile is described, and a few simple steps to run a model are given. Note that this concerns the ComboCode class, which runs the dust and gas radiative-transfer codes iteratively. If you are interested in the EnergyBalance class, which runs independently from the ComboCode class, you can skip this section.

### Combined dust and gas radiative transfer
The ComboCode class is an interface that provides access to two numerical RT codes for dust and gas respectively. The way these codes are linked through ComboCode is illustrated in the schematic below. Note that this schematic currently does not include iteration between energy balance and line RT in GASTRoNOoM, as this functionality is not yet implemented in ComboCode. Currently, the EnergyBalance class implemented in the ComboCode package is also not included as part of the ComboCode class.

![](https://github.com/robinlombaert/ComboCode/blob/master/aux/flow_chart_codes.png?raw=true)

<a name="inputfile">
### The ComboCode inputfile
</a>
WiP. See the inputfile for detailed description of the input parameters. A limited overview of the syntax is given in <a href="https://github.com/robinlombaert/ComboCode/blob/master/aux/info_inputComboCode.dat">the info_inputComboCode.dat</a> document.

<a name="how">
### How do I run ComboCode?
</a>
ComboCode is ran easily from a Python or iPython shell. Open one, and take the following steps, for an arbitrary input filename: 

    >>> #-- Import ComboCode into the Python shell:
    >>> from cc import ComboCode
    >>>
    >>> #-- Create a ComboCode object for the default input file: 
    >>> filename = '/Users/user_name/ComboCode/aux/inputComboCode.dat'
    >>> c1m = ComboCode.ComboCode(filename)
    >>>
    >>> #-- Start the ComboCode session:
    >>> c1m.startSession()

The c1m.startSession() command is essentially the body of the modeling calculation that is done. In what follows running this command is referred to as a "CC session". You cannot run this command twice. If you want to re-run a given filename, always create the ComboCode() object first (step 1) and then start the session (step 2).

It is possible to have multiple such CC sessions running concurrently in different shells. <a href="Manual.md#cleaning-your-databases">The databases make sure no conflicts can happen</a> between models requested in one session and another session, in case they are identical. 

Note that ComboCode, GASTRoNOoM and MCMax all print output to the shell. ComboCode-specific output typically concerns information about whether model calculation was successful, some information about the model output, etc. and is denoted with asterisks '*'. GASTRoNOoM and MCMax print output in their own formats and are preceded by an indication of which code is about to be called, and followed by "DONE!", so it is generally easy to tell which code is giving you information at any given time. 

Running the same ComboCode inputfile again after successful calculation will reload the models without calculating them, in which case no GASTRoNOoM or MCMax output will be printed to the shell. Whenever a model is calculated or is found in the database, the model ID is printed out. These IDs can differ between the MCMax and the three subcodes of GASTRoNOoM, so make sure you are working with the correct model output. Using these IDs you can retrieve your model output from the PATH\_GASTRONOOM and PATH\_MCMAX subfolders. Several methods and classes are available and described in Section 6 to read and use this information.

<a name="object">
### The ComboCode object and its model contents
</a>
The ComboCode object contains all of the information used when setting up and calculating models, as well loading data and doing the statistical analysis. All of these are accessible through properties of the ComboCode object. Below follow some examples, assuming that c1m is the ComboCode object, and the c1m.startSession() call has finished.

First, the model information is contained in a list of Star() objects. These objects essentially function as a dictionary, and contain keyword-based information about the models. Not only information given in the inputfile can be retrieved this way, but also information derived from it. So in case a keyword is not present in the dictionary, the Star() object will attempt to create it first, and return the value if successful. Otherwise, Star() objects will behave as a dictionary. 

    >>> #-- Retrieve the first model Star() object
    >>> s = c1m.star_grid[0]
    >>>
    >>> #-- This contains L_STAR and T_STAR from the CC inputfile. Ask for R_STAR
    >>> for k in ['T_STAR','L_STAR','R_STAR']:
    >>>     print s.has_key(k)
    >>>     print s[k]
    >>> 
    >>> #-- Retrieve the model IDs associated with this Star()
    >>> mid = s['LAST_MCMAX_MODEL']
    >>> gid = s['LAST_GASTRONOOM_MODEL']

Star() objects also contain several methods to read model output, such as temperature profiles for gas and dust, create the appropriate filename for a given dust species, or create the filename for the GASTRoNOoM thermodynamical output from the *cooling* subcode. Refer to the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.objects.Star.Star-class.html">online documentation</a> for more information on the Star class. Typically, you will be interested in any get\* or read\* methods. All of the calcKEY methods are for creating entries in the dictionary, so those values are available through **s['KEY']**. Any other methods are specifically called when needed by the package itself, so should not be called by the user, especially when the method name is preceded by \_\_.
    
    >>> #-- Retrieve the filename that contains the dust properties of amorphous carbon
    >>> fnd = s.getDustFn('AMCCDEPREI')
    >>>
    >>> #-- Retrieve a list of all dust species in the model
    >>> all_species = s.getDustList()
    >>>
    >>> #-- Read the radial density profile for carbon in AU, and average out the theta coordinate
    >>> rd = s.getDustRad('AMCCDEPREI',unit='AU')
    >>> rho = s.getDustDensity('AMCCDEPREI',avg_theta=1)
    >>>
    >>> #-- Get the CO and H2 number density. First set the filename arguments for the cool1 file for 12CO
    >>> kwargs = {'ftype':1,'mstr':'12C16O'}
    >>> rg = s.getGasRad(unit='cm',**kwargs)
    >>>
    >>> #-- The cool1 file contains the number density for its molecule, as well as H2
    >>> nmol = s.getNumberDensity(molecule=1,**kwargs)
    >>> nh2 = s.getNumberDensity(molecule=0,**kwargs)
    
The Star() objects also contain the information specific to molecules and transitions. Both are contained in their own type of object as well (Molecule() and Transition() respectively), and offer a slew of methods to get the molecule or transition properties, retrieve model information such as calculated line profiles, and compare with data. The objects have a string representation and can be used to check between objects for equality. For more information on the MlineReader and SphinxReader objects, see Section 6.

    >>> #-- Retrieve the list of Molecule() objects, or retrieve a specific molecule
    >>> molecs = s['GAS_LIST']
    >>> for m in molecs: print m
    >>> some_molec = s.getMolecule('12C16O') 
    >>>
    >>> #-- If used in context of GASTRoNOoM RT, the RadiatReader object is available (see Section 6)
    >>> rr = some_molec.radiat
    >>> print rr['levels']
    >>>
    >>> #-- And if the model was already calculated, the *mline* model ID is here:
    >>> mline_id = some_molec.getModelId()
    >>>
    >>> #-- Or the MlineReader object:
    >>> some_molec.readMline()
    >>> ml = some_molec.mline
    >>>
    >>>
    >>> #-- Retrieve the list of Transition() objects, or retrieve all transitions for a specific molecule
    >>> trans = s['GAS_LINES']
    >>> for t in trans: print t
    >>> some_trans = s.getTransitions('12C16O')
    >>> t1 = some_trans[0]
    >>>
    >>> #-- The Transition object is the workhorse of line modeling. 
    >>> trans_id = t1.getModelId()
    >>> frequency = t1.getFrequency()
    >>> Eup = t1.getEnergyUpper(unit='K')
    >>>
    >>> #-- Make a dictionary that contains all relevant info (is also used in database)
    >>> tdict = t1.makeDict()
    >>>
    >>> #-- Get the sphinx (final GASTRoNOoM subcode for ray tracing) filename that contains the line profile
    >>> fnsph = t1.makeSphinxFilename(number=2,include_path=1)
    >>>
    >>> #-- Get the SphinxReader object that read the sphinx output files, or the line strength
    >>> t1.readSphinx()
    >>> sph = t1.sphinx
    >>> ls = t1.getIntIntIntSphinx(units='si',cont_subtract=1)
    >>>
    >>> #-- If taken from ComboCode object, will contain data as well
    >>> dls_tmb = t1.getIntTmbData(units='tmb')
    >>> dls_si = t1.getIntTmbData(units='si')
    >>>
    >>> #-- Calculate the loglikelihood for the model compared to the data
    >>> lll = t1.getLoglikelihood(vmin=-10.0,vmax=10.0)

Finally, the ComboCode object also contains multiple dictionaries for each of the requested STAR\_NAME input key values. These contain either data or statistics objects for each of the given star names. The example here does not contain data, and has STAR\_NAME=model. But assuming the inputfile was set up for, e.g., 'rdor', the data and statistics can be retrieved as follows. 

    >>> #-- If PACS is on, retrieve the Pacs() object for R Dor and list the data files and the line fit results:
    >>> pacs = c1m.pacs['rdor']
    >>> print pacs.data_filenames 
    >>> print pacs.linefit
    >>> 
    >>> #-- Statistics can be calculated for the comparison between PACS measured line strengths and the model
    >>> stats = c1m.pacsstats['rdor']
    >>>
    >>> #-- Get the line strength model/data ratios for unblended lines for the first model
    >>> pacs_id = c1m.star_grid[0]['LAST_PACS_MODEL']
    >>> wav = stats.getRatios(sel_type='int_ratios',data_type='central_mwav',this_id=pacs_id)
    >>> ratio = stats.getRatios(sel_type='int_ratios',data_type='int_ratios',this_id=pacs_id)
    >>> error = stats.getRatios(sel_type='int_ratios',data_type='int_ratios_err',this_id=pacs_id)
    
Note that the data and statistics objects are usually prepared with some basic form of data and statistics settings as long as the respective keywords in the inputComboCode.dat file are turned on (such as PACS=1, STATISTICS=1). For more details on what is possible, check the <a href="http://robinlombaert.github.io/ComboCode/">online documentation</a> for ComboCode.py, and the respective data and statistics objects. 

## 5. Data management
Filename conventions!

<a name="resolved">
### Resolved molecular emission (radio)
</a>

<a name="unresolved">
### Unresolved molecular emission (infrared)
</a>

### Spectral energy distribution

## 6. Model management
Most of what is described here concerns the I/O output used in the entire package of ComboCode. As an example, the independent EnergyBalance module makes use of several of the spectroscopy readers described in this section.

### Reading and using model output
ComboCode includes several tools to read and use model input and output data, depending on the involved numerical codes (GASTRoNOoM mostly, but some functionality for MCP/ALI is given as well). All of the so-called **Reader** objects are found in the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.tools.readers-module.html">cc.tools.readers module</a>. This includes: 
- TxtReader, FitsReader: radio data (single resolved emission lines) in the form of fits and txt files (see above)
- SphinxReader: GASTRoNOoM-sphinx output with the ray-tracing results and the predicted line emission profiles
- MlineReader: GASTRoNOoM-mline output with the radiative-transfer model results, including line opacities, source function, scattering integral, level populations, and the molecular spectroscopy used for the model
- CollisReader, RadiatReader: the collisional rate data and molecular spectroscopy input, respectively, for GASTRoNOoM
- PopReader: level populations calculated with MCP/ALI
- LamdaReader: the collisional rate data and molecular spectroscopy from Lamda-format files as input for MCP/ALI
- KappaReader: reading and interpolating dust opacities used for MCMax

![](https://github.com/robinlombaert/ComboCode/blob/master/aux/flow_chart_readers.png?raw=true)

With the exception of KappaReader and LineList (entirely stand-alone), all classes inherit from the Reader base class, which itself functions as a dictionary. Instances of these classes thus can be treated as a dictionary that contains the information from the files. Instances of almost all classes can be created by passing them the filename of the output/input file. The only exceptions to this are KappaReader and RadiatReader (see below). Here follows an example that prints out the molecular excitation levels included in an mline radiative-transfer model: 
    
    >>> #-- Import and read
    >>> from cc.tools.readers import MlineReader
    >>> ml = MlineReader.MlineReader('ml1model_2016-07-26h15-22-11_12C16O.dat')
    >>> print ml['level']

Several helper functions are available for each class to retrieve specific information and create plots. Moreover, classes that read the same information for different codes (such as molecular spectroscopy for GASTRoNOoM and the Lamda format) employ the same syntax and methods to retrieve that information regardless of source, so they can be used interchangeably. 

<a name="spec">
#### 1. Molecular spectroscopy and radiative-transfer output
</a>
Molecular spectroscopy files are used in two formats: one for GASTRoNOoM (nonstandard format) and one for MCP/ALI (which is in the format of the Leiden Atomic and Molecular Database: <a href="http://home.strw.leidenuniv.nl/~moldata/">Lamda</a>). Most of the spectroscopy information is also contained in the mline output for GASTRoNOoM (this excludes collision rates).
    
RadiatReader, MlineReader and LamdaReader all inherit from SpectroscopyReader and MolReader, and thus share its methods. CollisReader also inherits from SpectropscopyReader, sharing the transition index methods. MlineReader inherits the capabilities of PopReader for level populations. LamdaReader inherits from CollisReader for collision rates. A few examples: 

    >>> #-- Import and read
    >>> from cc.tools.readers import CollisReader, LamdaReader, MlineReader
    >>> lr = LamdaReader.LamdaReader('co@rovib-yang.dat')
    >>> ml = MlineReader.MlineReader('ml1model_2016-07-26h15-22-11_12C16O.dat')
    >>> cr = CollisReader.CollisReader('12C16O_collis.dat')
    >>>
    >>> #-- Retrieve all collisional transition indices for a given lower level index
    >>> lr.getTI(llow=1,itype='coll_trans')
    >>> cr.getTI(llow=1)
    >>>
    >>> #-- Retrieve the collision rates for a given collisional transition (both ways 
    >>> #   apply to both objects), and the temperature grid
    >>> lr.getRates(index=23)
    >>> lr.getTemp()
    >>> cr.getRates(llow=1,lup=3)
    >>>
    >>> #-- Retrieve the spectroscopic transition index for a lower and upper index
    >>> lr.getTI(llow=1,lup=2,itype='trans')
    >>> ml.getTI(llow=1,lup=2)
    >>>
    >>> #-- Retrieve the energy of a given level
    >>> lr.getLEnergy(llow=12,unit='K')
    >>> ml.getLEnergy(llow=12,unit='erg')
    >>>
    >>> #-- Retrieve the level weights of a list of indices in array form
    >>> level_indices = [2,3,4,5,6,7,8]
    >>> weights = lr.getWeights(index=level_indices)
    >>>
    >>> #-- Retrieve the level populations. Then make a plot of all level populations.
    >>> ml.getPop(index=20)
    >>> ml.plotPop()
    >>>
    >>> #-- Create interpolators for level populations (similarly for collision rates in CollisReader)
    >>> ml.setInterp(itype='spline',k=1,ext=3)
    >>> interp = ml.getInterp(index=1)
    >>>
    >>> #-- Retrieve the modeled source function, line opacity at line center and scattering integral
    >>> #   as a function of transition index (dict key), and impact parameter in cm.
    >>> p = ml.getP()
    >>> for dtype in ['sf','lo','si']:
    >>>     print(ml[dtype])

       
There are many more methods available, described in the docstrings of each class and its methods, see <a href="http://robinlombaert.github.io/ComboCode/">the online documentation</a>. GASTRoNOoM-mline produces three output files: ml1, ml2 and ml3. Only ml1 (spectroscopy and circumstellar properties) and ml3 (radiative-transfer output for all transitions, and level populations) are read by MlineReader. It doesn't matter which of the two filenames are passed to MlineReader; both files will be read. Ml2 files give an overview of the iteration steps. MCP/ALI contrarily give less radiative-transfer output by default, but do give level populations (.pop files) to be read with PopReader.

Note that RadiatReader requires additional input, namely the number of transitions and levels included. Refer to the MOLECULE definition in the ComboCode input for GASTRoNOoM, where ny=ny_up+ny_low and nline are given. This is not needed for the MlineReader.

    >>> #-- Import and read
    >>> from cc.tools.readers import RadiatReader
    >>> rr = RadiatReader.RadiatReader('12C16O_radiat.dat',nline=240,ny=122)
        
Finally, a module is available to read other generic line lists from online databases such as
    
- Jet Propulsion Laboratory: <a href="http://spec.jpl.nasa.gov/">JPL</a>
- Cologne Database for Molecular Spectroscopy: <a href="https://www.astro.uni-koeln.de/cdms">CDMS</a>
    
This module does not inherit from the Reader base class and has its own functionality. 

    >>> #-- Import and read. Note the additional options.
    >>> from cc.tools.readers import LineList
    >>> ll = LineList.LineList('12C16O_CDMS.dat',x_min=100.,x_max=2000.,unit='GHz')
    >>> 
    >>> #-- Retrieve the linelist
    >>> lines = ll.getLineList()
    >>>
    >>> #-- Make Transition() objects for additional functionality
    >>> from cc.modeling.objects import Molecule
    >>> mm = Molecule.Molecule('12C16O',linelist=1)
    >>> transitions = ll.makeTransitions(molecule=mm,telescope='PACS')
    >>> for t in transitions: print(t)
    
Line lists from the Lamda online database can of course be read with the LamdaReader. 

<a name="raytrace">
#### 2. Modeled line profiles and ray-tracing output for GASTRoNOoM
</a>
The third subcode of GASTRoNOoM, named sphinx, provides the ray tracing and calculated the line emission profiles. The output consists of two files: sph1 and sph2, the former giving intensities as a function of impact parameter, the latter giving the intrinsic and beam-convolved line profiles in several units as a function of velocity. It doesn't matter which of the two filenames are passed to SphinxReader; both files will be read. Several methods are available to return the information. 

    >>> #-- Import and read
    >>> from cc.tools.readers import SphinxReader
    >>> fn = 'sph1model_2016-07-26h14-37-10_12C16O_vup0_jup1_vlow0_jlow0_JCMT_OFFSET0.00.dat'
    >>> sph = SphinxReader.SphinxReader(fn)
    >>>
    >>> #-- Retrieve the impact parameter grid and the normalized intensity
    >>> p = sph.getImpact()
    >>> intens = sph.getNormalizedIntensity()
    >>>
    >>> #-- Retrieve the intrinsic line profile
    >>> vel_intrinsic = sph.getVelocityIntrinsic()
    >>> lp_int_cgs = sph.getLPIntrinsic(cont_subtract=1)
    >>>
    >>> #-- Retrieve the beam-convolved profile in Kelvin
    >>> vel = sph.getVelocity()
    >>> lp_K = sph.getLPIntrinsic(cont_subtract=0)
        
Other methods are available for other types of information. The sph object also functions as a dictionary, so can be easily checked. 

<a name="dust">
#### 3. Dust opacities for MCMax
</a>
The MCMax dust opacities can be read with the KappaReader. This is the only Reader object that doesn't inherit from the Reader base class and has its own functionality. It depends heavily on the usr/Dust.dat file, that contains all the dust information associated with the opacity files and their name tags. 

    >>> #-- Create a KappaReader object.
    >>> from cc.tools.readers import KappaReader
    >>> kr = KappaReader.KappaReader()
    >>>
    >>> #-- Read the information of a dust species (e.g. amorphous carbon, included in usr.dist/Dust.dat)
    >>> species = 'AMCCDEPREI'
    >>> kr.readKappas(species)
    >>>
    >>> #-- Retrieve the information.
    >>> w = kr.getWavelength(species)
    >>> kappas = kr.getKappas(species,index=0)
    >>> ext_eff = kr.getExtEff(species,index=1)
    >>>
    >>> #-- Create an interpolator (spline), pass extra args k and ext to the spline.
    >>> interp = kr.interpolate(species=species,index=0,k=3,ext=0)
        
The opacity files contain three types of information as a function of wavelength, given by the index: 0: extinction, 1: absorption, 2: scattering.

<a name="general">
#### 4. General model input/output for GASTRoNOoM and MCMax
</a>
Some information is available through general methods that retrieve specifically requested information, rather than working through a Reader object. Several methods are available to read any type of multiple-column-based or 1-column-based model input/output. Examples are getGastronoomOutput, getInputData and getKeyData in the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.tools.io.DataIO-module.html">cc.tools.io.DataIO module</a>. See link for more information on the methods and how to use them. Some examples. 

    >>> #-- Import modules
    >>> from cc.tools.io import DataIO
    >>> import cc.path
    >>> import os
    >>> 
    >>> #-- Retrieve the radius and temperature for a GASTRoNOoM model (multi column output file)
    >>> fgr_file = 'coolfgr_allmodel_2016-07-26h14-37-10.dat'
    >>> rad = DataIO.getGastronoomOutput(filename=fgr_file,keyword='RADIUS',return_array=1)
    >>> temp = DataIO.getGastronoomOutput(filename=fgr_file,keyword='TEMP',return_array=1)
    >>>
    >>> #-- Retrieve the radius and temperature for a MCMax model (1 column output file)
    >>> dens_file = 'denstemp.dat'
    >>> incr = 150 # NRAd -- number of radial points
    >>> rad = DataIO.getKeyData(incr=incr,keyword='RADIUS',filename=dens_file)
    >>> incr = 450 # NRAD*NTHETA -- number of radial times theta points
    >>> temp = DataIO.getKeyData(incr=incr,keyword='TEMPERATURE',filename=dens_file)
    >>> 
    >>> #-- Read a plot cfg file as a dictionary
    >>> inputfile = os.path.join(cc.path.aux,'plot_config_example.cfg')
    >>> dd = DataIO.readDict(filename=inputfile,convert_lists=1,convert_floats=1,\
    >>>                      comment_chars=['#','!'])
    
Especially the MCMax and GASTRoNOoM model output can also be more easily accessed through the Star() objects that represent each model in the model grid. They are accessible from the c1m.star_grid list (see above for more details).

### Plotting model output
With the many reader classes and methods, the model output is always easily accessible for the user. However, ComboCode does provide some built-in plotting modules to plot modeled emission lines, SEDs, circumstellar profiles for temperature, density, etc. These modules can be used to make publication-quality plots through the use of configuration files (or cfg files). Most plotting methods have the option to pass such a cfg file. 

The most obvious way to plot model results (along with the data if available) is through the PLOT/CFG keywords listed at the bottom of the inputComboCode.dat file. Turning one of these on will produce a plot at the end of the CC session, and print out the filename of the plot. If you have loaded some models with the ComboCode object, and then wish to remake the plots after changing some parameters in the cfg files, then you can run c1m.runPlotManager() and the plots will be made anew. 

An example cfg file is available in cc/aux/plot\_config\_example.cfg. Any argument of the Plotting2.plotCols (or Plotting2.plotTiles) can be added in the cfg file. It converts lists and strings for you, through the use of the DataIO.readDict method. You can experiment with that method to see what is possible in terms of input. But most settings of the plots can be adapted this way. Some specific plotting methods also allow passing arguments through the cfg file and are used before calling the plotting method itself. 

An overview of all possible plot settings is available at <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.plotting.Plotting2-module.html#plotCols"> Plotting2.plotCols</a> and at <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.plotting.Plotting2-module.html#plotTiles"> Plotting2.plotTiles</a>.

## 7. Database Management
ComboCode provides basic database functionality to track and book-keep modeling output and data files. Any model that is calculated successfully is listed in the appropriate model database with its parameters and is assigned a model ID. Whenever ComboCode is ran, and requested models are found in the database, the model ID is returned instead of calculating the model anew. 

### Types of databases


### Database locking


### Cleaning your databases
The databases include a way to track which models are currently being calculated in **any** CC session (or, shell). This works through an "IN\_PROGRESS" entry in the keyword definition of a given cooling, mline, sphinx, or MCMax model, and regular synchronisation between the Database() instance in the CC session and the harddisk version of the database. 

However, if for some reason a CC session is terminated and results in an error, it may be possible "IN\_PROGRESS" entries remain in the database, while no model is currently being calculated. You can clean databases off these left-over "IN\_PROGRESS" entries by doing the following.  Only do this if you are positive **no other CC session is currently running**! Imagine a CC session is waiting for an mline model to finish (and no other CC session is running), open an ipython shell and do the following (for PATH_GASTRONOOM=MyModels):

    >>> from cc.tools.io import Database
    >>> db_fn = '/Users/user_name/GASTRoNOoM/MyModels/GASTRoNOoM_mline_models.db'
    >>> Database.cleanDatabase(db_fn)
    >>> exit()

This effectively removes all "IN\_PROGRESS" entries from the databases. It is possible that you run a CC session, which ends up waiting for another CC session to finish, while no other CC session is currently running. This means such a left-over "IN\_PROGRESS" entry is encountered. Open a separate ipython shell, and run the above script. Once finished, the CC session will continue (it will say the mline model failed, since it is no longer present in the database). You can re-run ComboCode if you want to re-try the model. You can run the exact same script for other databases, including cooling, sphinx and MCMax. 

Note that older versions of the databases may sometimes contain model\_ids for mline and sphinx (in their respective databases) that contain no molecules or transitions. These are not allowed anymore in the current version of ComboCode. Running this method also removes those empty model\_ids. This must be done in case you encounter one of the two following error messages: 

    KeyError: 'Empty molec id found in the database. This should not be possible.'
    KeyError: 'Empty trans id found in the database. This should not be possible.'

## 8. The EnergyBalance
The <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance-module.html">EnergyBalance module</a> is a standalone numerical implementation of the energy balance in stellar winds, written entirely in Python. The module consists of the differential equation dTdr and the class EnergyBalance, which is what you will use in practice. 

This section includes:
- an overview of how to set up the EnergyBalance (EB) class with input, how to run the EnergyBalance and solve the differential equation, and how to use and plot the results. 
- the class makes heavy use of the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.profilers-module.html">profilers module</a>, warranting an in-depth discussion of how the Profiler and its inheriting classes work. 
- the EnergyBalance can run iteratively with the radiative-transfer code ALI. The <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.codes.ALI-module.html">ALI module</a> provides a method to do this for you, and shows how you can set up your own script as well.

Note that the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.tests-module.html">tests module</a> contains a lot of examples regarding profilers, the EnergyBalance class, etc. with examples.

### Setting up and running the EnergyBalance
When calculating the energy balance using EB, you essentially take two steps, in addition to whatever plotting you want to do. This example uses the default settings (with exception of imax=10, the maximum amount of iterations allowed).

1) Create the EB object with the necessary input
    >>> from cc.modeling.physics import EnergyBalance as EB
    >>> eb = EB.EnergyBalance()
2) Iterate the temperature profile until convergence or until 10 iterations have been done
    >>> eb.iterT(imax=10,warn=0)
3) Plot, e.g., the temperature profiles of all iterations, and the cooling and heating rates of the final iteration
    >>> eb.plotT()
    >>> eb.plotRates()

If the T calculation did not yet converge in the 10 iterations calculated in the above example (imax can be defined as whatever you want), you can continue iterating by calling the iterT method again with a larger imax. EB leaves off from the last result, and simply continues. This means that you can also update the level populations at this stage for the next set of iterations. 
    
    >>> eb.updatePop(m='12C16O',fn='new_CO_populations.pop')
    >>> eb.iterT(imax=20)
    
The <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html">online documentation</a> gives an overview of all the class methods of EB. Most of these methods are called when needed by the class itself, and you need not bother with it, unless you want to extract specific information. The only methods you will use in practice are iterT, calcT (if you only wish one single iteration), updatePop, and the plot methods listed at the end. For detailed information on the different keywords for each class method, I refer to the online documentation.

<a name="EBinput">
#### 1. Defining input for the EnergyBalance
</a>
All the static model settings (numerical as well as physical) are provided when the EB class is called. There are three instances that provide input for the EB class when it is created. 
    
1) Three templates comes with ComboCode: standard (standard), gastronoom, and mcp. The standard template gives the most generalised settings possible, while gastronoom and mcp mimic the settings of the respective codes exactly. For MCP, note that not all its options are available in EB, e.g., for the drift velocity, only the power law option can currently be used. A template can be chosen by setting template='mcp'. The templates are available in the aux/inputEnergyBalance_\***.dat files, and give an example for the next option. The templates should not be changed! 

2) A filename can be given to the EnergyBalance as fn='filename.dat'. The structure of the file follows that of the templates, but does not have to include all the keywords in the template. You can define your keywords there at will. Note that the template (standard, if it has not been specified) will always be read. The options in the file that you specify will then overwrite the template settings. 

3) You can pass any of the keywords given in the templates explicitly to the EnergyBalance class. They will overwrite both the template as well as the file if it was specified. As an example, a=[0.005e-4,0.25e-4,100,1] sets the grain size grid with minimum and maximum sizes 0.005 and 0.25 micron, and 100 points distributed in logarithmic scale. 
    
The following example creates an EB based on the MCP template, adjusting for keywords in your own inputfile, and setting the grain size to a single value of 0.01 micron, and including only two heating terms: dust-gas collisions and thermal accommodation.

    >>> eb = EB.EnergyBalance(template='mcp',fn='my_own_inputfile.dat',a=0.01e-4,hterms=['dg','dt'])

Going through the standard inputfile template you will find several blocks of input. We consider here each block.

1) The first block gives the coordinate grids for the radius (r, in cm), the grain  size (a, in cm), and the wavelength (l, in cm). They are defined as a list, giving the input for the method <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.tools.numerical.Gridding-module.html">Gridding.makeGrid</a> in a fixed order: [min value, max value, number of points, linear (0) or logarithmic (1) distribution, floats (0) or integers (1)]. Only the first three values are mandatory. The last two values are 0 by default. Note that the grain size can be given as a single value (without brackets), in which case EB treats it as an average grain size much the same way as MCP. 

    r = [1e14,1e17,1000,1]
    a = [0.005e-4,0.25e-4,200,1]
    l = [0.1e-4,1000e-4,1000,1]

2) The second block defines the different relevant profiles for the gas expansion velocity, the initial kinetic gas temperature, the opacity, the gas and dust mass-loss rates, and the dust temperature. These profiles are defined through the use of the Profiler class, and are discussed in <a href="Manual.md#profilers">the next section</a>. Note that the velocity and the two mass-loss rate values can be given as a single value as well, in which case they are kept constant throughout the entire stellar wind. 

    v         = vbeta r0=1e14 v0=vs vinf=10e5 beta=1 vi_mode=beta rstar=3.32475e13 beta_inner=0.5
    Tinit     = Teps T0=1443. r0=1e14 epsilon=0.5 inner=1
    opac      = read_opacity species=AMCDHSPREI index=1 k=5 ext=3
    mdot      = 1e-7
    mdot_dust = 2e-10
    Td        = Tdust T0=2500. r0=1.3e13 s=1

3) The third block defines the stellar temperature and radius, and the type of stellar spectrum. Ltype also makes use of the Profiler class, but typically is set by a single function name. The function should be present in the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.profilers.Radiance-module.html">Radiance module</a>, where you will currently only find 'blackbody' as an option. Not to worry, profilers can always read files as well, so that is an option as well. 

    Tstar = 2500.
    rstar = 1e13
    Ltype = blackbody

4) The fractional atomic hydrogen and helium abundances are defined next. Based on these values the mean molecular weight is calculated. If you want the mean molecular weight to be different from what fH and fHe would result in, you can define mu=1.4 in the inputfile or in the EB.EnergyBalance() call, for instance.

    fH  = 0.0
    fHe = 0.0

5) The next block concerns dust properties, specifically the average specific density, the porosity, and the grain size distribution. For GSD currently only MRN is available. The function should be defined in the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.profilers.Grainsize-module.html">Grainsize module</a>. GSD is ignored if the variable **a** is a single value. 

    sd  = 1.8
    P   = 0.0
    GSD = MRN

6) The sixth block sets the heating and cooling terms: cterms can include line cooling (lc), and h2 cooling (h2); hterms can include dust-gas collisions (dg), thermal accommodation (dt), photoelectric heating (pe), and cosmic rays (cr). The adiabatic coefficient gamma can be set to a single value, or be determined by the code based on H2 properties (h2, makes use of the data in aux/h2\_physical\_properties.dat) depending on the temperature of each iteration (since gamma technically depends on T). The option gamma=h2_mcp is also available, which sets the adiabatic coefficient to two different values based on the T regime: 5/3 at T < 350K, 7/5 otherwise. 

Lastly, the method used for some of the heating and cooling terms can be specified. 
    
   a) Photoelectric heating (pe_method, see the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html#Hpe">Hpe method</a> for more details): Either bakes or draine. For bakes, G0 factor, and minimum grain size scaling, and the C and O abundances with respect to H have to be specified. See the gastronoom template for an example. For draine the Kpe factor has to be speciefied. See the mcp template for an example.

   b) Cosmic ray heating (cr_method, see the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html#Hcr">Hcr method</a> for more details): Either groenewegen or standard. The method is that of Goldsmith and Langer 1978 in the standard case. The same applies to the groenewegen case, but then all gas mass is places into H2 even if fH or fHe are non-zero. 

   c) H2 cooling (he_method, see the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html#Ch2">Ch2 method</a> for more details): Either groenewegen or decin. 

   d) Dust-gas collisional and thermal accommodation heating (heatmode, see the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html#Hdg">Hdg method</a> and the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.physics.EnergyBalance.EnergyBalance-class.html#Hdt">Hdt method</a> for more details): Either classical or gs2014. Currently, only classical is fully implemented, while gs2014 will only use the Gail & Sedlmayr 2014 equations for dust-gas collisions and thermal accommodation. The differential equation dTdr and all other heating and cooling terms are treated the same as in classical mode. 

    cterms = [lc,h2]
    hterms = [dg,dt,pe,cr]
    heatmode = classical
    gamma = h2
    pe_method = draine Kpe=1e-26
    cr_method = standard
    h2_method = decin

7) The next block deals with all the gas-dust drift parameters: w_mode is either standard (derives the drift velocity profile based on the balance between radiation pressure and drag force), or vbeta2D (derives the terminal drift the same way as the standard mode, but uses a power law for the inner wind, requiring a beta to be defined as well, see the mcp template). When w\_mode is standard, the thermal component can be chosen, see the w\_thermal variable in the <a href="http://robinlombaert.github.io/ComboCode/ComboCode.cc.modeling.profilers.Velocity-module.html#driftRPDF">driftRPDF method</a> (where epstein is the recommended mode, note that kwok is actually incorrect). Finally, the maximum drift velocity above which dust grains sputter and are destroyed can be given as w\_sputter, and the sticking coefficient is given as alpha.

    w_mode = standard
    w_thermal = epstein
    w_sputter = 20e5
    alpha = 0.0
    
8) Finally, the last block defines the radiative-transfer code from which level populations are used, the molecules included, and their respective collision rates and level populations. These are only used if lc is included in cterms in block 6. 

When defining the molecule key(s), give one line for each molecule you wish to include. The line gives the name of the molecule first (free to choose whatever name you wish, except for CO which should always be 12C16O), then the input for the abundance profiler (a single value for a constant abundance, or a function). This is in case lc is NOT included in cterms in block 6: the CO abundance must be provided, because it is required information for the photoelectric heating term in case of pe_method=bakes.

The inputfile may look something like this for the last block, if lc is NOT included in cterms: 

    #-- Line Cooling parameters not required. Still, give molecule with name and abundance info for 12C16O
    rtcode = gastronoom  # Not used if lc is not included
    molecule = 12C16O np.loadtxt fname=CO_abundance.dat usecols=[1,2] skiprows=9 unpack=1
    collis = []
    pop = []

When lc is included in the cterms in block 6, level populations must be given, and the abundance information from the radiative-transfer output is extracted instead. Then no abundance information has to be given on the molecule line. The order in which the molecule keys are given is important, as that associates the molecules with their respective collision rate files (given as a list in collis) and the level population files (given as a list in pop).

The inputfile may look something like this for the last block, if lc is included in cterms: 
    
    #-- Line Cooling parameters: Only included if lc in cterms. molecule contains name and abundance info
    rtcode = ali
    molecule = 12C16O 
    molecule = 1H1H16O
    collis = [co@rovib-yang.dat,o-h2o@faure2008.dat]
    pop = [my_CO_model.pop,my_H2O_model.pop]
    
In practice, you will very often define your molecules and molecular data in-line when calling the EnergyBalance. In the above examples, your call would look something like this: 
    
    >>> #-- No line cooling included
    >>> molecule = '12C16O np.loadtxt fname=CO_abundance.dat usecols=[1,2] skiprows=9 unpack=1'
    >>> eb = EB.EnergyBalance(template='mcp',molecule=molecule,cterms=[h2])
    >>> 
    >>> #-- Line cooling included, with two molecules (note that the filenames should include the full path): 
    >>> collis = ['co@rovib-yang.dat','o-h2o@faure2008.dat']
    >>> pop = ['my_CO_model.pop','my_H2O_model.pop']
    >>> eb = EB.EnergyBalance(template='mcp',molecule=['12C16O','1H1H16O'],cterms=[h2,lc],pop=pop,collis=collis)
    
<a name="EBiter">
#### 2. Iterating the EnergyBalance
</a>


Any iteration specifications (such as number of iterations, convergence criterion, etc.) are given to the iterT method. 

### Profilers

### Iteration with ALI


## 9. Statistical methods

### Measuring goodness of fit

### Trend analysis




## 10. Additional modules

### Line profile fitting

### Plotting line lists




