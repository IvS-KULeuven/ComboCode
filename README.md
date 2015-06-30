# Welcome to the ComboCode README file
If you want to get started with ComboCode, read this file first. That's why it's the README file (u-hu).

## Introduction & User Manual
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
<li>Examples of how to run ComboCode will be added at a later time.</li>
</ol>

### The codes currently in ComboCode
* <b> GASTRoNOoM:</b> L. Decin (KU Leuven, Belgium) 
* <b> MCMax:</b> M. Min (UvA, the Netherlands)



## Installation
First and foremost, if you have python installed, with all necessary dependencies (Anaconda distribution, pyfits, and the IvS repository), you can download the code to your hard disk right away: 
* Clone the git repository to create a local copy, located in ~/ComboCode/:
    - $ cd 
    - $ git clone https://github.com/IvS-KULeuven/ComboCode.git ComboCode

* Updating your own clone of ComboCode to the most recent version can be done with:
    - $ cd ~/ComboCode/
    - $ git pull

Note that this does not include the ~/ComboCode/Data/ files. You can copy those over from <a href="http://ster.kuleuven.be/~robinl/cc/Data/"> http://ster.kuleuven.be/~robinl/cc/Data/</a>. Afterwards, you can change those files at your whim. They do not come as part of the repository.

A shell script for the installation of ComboCode will be provided soon as part of the repository. It will include: 
* Anaconda python distribution installation, including additional packages required to run ComboCode
* Downloading additional Data files not part of the repository but required to run ComboCode
* Creating a master branch in your local git repository that can be used to run the code

This is a work in progress!

## Documentation
Up-to-date documentation that goes with the package is available on GitHub at:

<a href="https://robinlombaert.github.io/ComboCode"> ComboCode Documentation</a>

## Requirements
TBD

## Developer's Manual
If you want to make changes to ComboCode, you should fork the repository to your own github account so you can create your own branch on your machine that will not interfere with the main master branch. 

* Go to the main IvS-KULeuven/ComboCode.git page and fork the repository to your account. 

* Clone a copy of the code in your account to your local machine (fill in your github username)
    - $ cd 
    - $ git clone https://github.com/<YOUR_USERNAME>/ComboCode.git ComboCode

* Creating your a developer's branch can be done by:
    - $ cd ~/ComboCode/
    - $ git branch mybranch
    - $ git checkout mybranch

* Once you've added changes to the code, and you want to suggest those changes to be included in ComboCode, you can commit the changes and push to the origin of your account on GitHub:
    - $ git commit -a -m "Document your changes! Either in-line, or without the -m tag so a text editor opens for documentation."
    - $ git push origin  mybranch

* Then you can create a pull request to the IvS-KULeuven account, and wait for the pull request to be included in ComboCode, as descreibed <a href="https://help.github.com/articles/creating-a-pull-request/"> here</a>.


This section is a work in progress!



