# Welcome to the ComboCode User Manual
## Introduction
### What is this manual?
This manual is meant as a guide to running ComboCode and its two numerical codes. This is not an all-inclusive, comprehensive manual for the ComboCode capabilities. However, in-depth, up-to-date documentation for the Python package is available on GitHub at <a href="https://IvS-KULeuven.github.io/ComboCode"> ComboCode Documentation</a>, as a collection of doc-strings. The source code is also available there. Together with the in-depth documentation and the cookbooks provided in this manual for extracting information and using the additional modules, you should be able to use the package to its full extent. For additional questions or remarks, please contact <a href="https://github.com/robinlombaert">R.~Lombaert</a>. 

For best results, it is recommended to use both radiative-transfer (RT) codes embedded in ComboCode: GASTRoNOoM (line RT) and MCMax (continuum RT), authored and maintained by L.~Decin and M.~Min, respectively. Both authors have to be contacted for use of the RT codes, also as part of ComboCode, and their contact details can be requested from <a href="https://github.com/robinlombaert">R.~Lombaert</a>.

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

### The ComboCode inputfile

### How do I run ComboCode?





## Data management

### Resolved Molecular emission (Radio)

### Unresolved molecular emission (Infrared)

### Spectral energy distribution




## Model management
### Combined dust and gas radiative transfer

### Reading and using model output

### Plotting model output

### Database management




## Statistical methods

### Measuring goodness-of-fit

### Trend analysis




## Additional modules

### Line profile fitting

### Plotting line lists




