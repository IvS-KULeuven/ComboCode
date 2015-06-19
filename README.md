# Welcome to the ComboCode README file
If you want to get started with ComboCode, read this file first. That's why it's the README file (u-hu).

## Introduction
### What is ComboCode?
ComboCode is a Python based package designed to work with radiative-transfer codes and the data they are meant to model. 
The radiative transfer is usually calculated for cool stellar winds of evolved stars, such as AGB stars.
The functionality includes:
* <b>Modeling</b>:
    - Currently works with GASTRoNOoM for gas radiative transfer, and with MCMax for dust radiative transfer
    - Allowing output of one code to be used as input for another
    - Databases for modeling output and easy parameter space searches
* <b>Data</b>: 
    - Management of data files associated with radio data, SEDs and spectroscopic data

### How do I run ComboCode?
<ol>
<li>First you have to install the package itself. See instructions below.</li>
<li>You will need the radiative-transfer codes to be able to run ComboCode to its fullest potential. For this you will
have to contact the authors of the codes.</li>
    * <b> GASTRoNOoM:</b> L. Decin (KU Leuven, Belgium) 
    * <b> MCMax:</b> M. Min (UvA, the Netherlands)
<li>Examples of how to run ComboCode will be added at a later time.</li>
</ol>



## Installation
A shell script for the installation of ComboCode will be provided soon as part of the repository. I will include: 
* Anaconda python distribution installation, including additional packages required to run ComboCode
* Downloading additional Data files not part of the repository but required to run ComboCode
* Creating a master branch in your local git repository that can be used to run the code

This is a work in progress!

## Documentation
Up-to-date documentation that goes with the package is available on GitHub at:
<a href="https://robinlombaert.github.io/ComboCode"> ComboCode Documentation</a>
