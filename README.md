# Welcome to the ComboCode README file
## Introduction & User Manual
### What is ComboCode?
ComboCode is a Python based package designed to work with radiative-transfer codes and the data they are meant to model. 
The radiative transfer is usually calculated for cool stellar winds of evolved stars, such as AGB stars.
The functionality includes:
* <b>Modeling</b>:
    - Currently works with GASTRoNOoM for gas radiative transfer, and with MCMax for dust radiative transfer
    - Allowing output of one code to be used as input for another
    - Automatic line selection based on available data and line listing
    - Databases for modeling output and easy parameter space searches
    - Interaction with a supercomputer clusters built into the databases
    - Statistical analysis for samples and individual sources, and both resolved and unresolved emission lines
* <b>Data</b>: 
    - Management of data files associated with radio data, SEDs and spectroscopic data
    - Fitting routines for resolved emission lines

### How do I run ComboCode?
<ol>
<li>First you have to install the package itself. See instructions below.</li>
<li>You will need the radiative-transfer codes to be able to run ComboCode to its fullest potential. For this you will
have to contact the authors of the codes.</li>
<li>Examples of how to run ComboCode will be added at a later time.</li>
<li>ComboCode folder structure: 
<ol>
<li>cc -- Contains the Python modules.</li>
<li>aux -- Contains auxiliary files that ComboCode requires, and they come as part of the repository. These files are not to be chanegd by the user.</li>
<li>usr -- Contains files that ComboCode requires but are not part of the repository. The settings are user specific.</li>
<li>usr.dist -- The blueprint for usr/ that is part of the repository. This folder is to be copied to usr/ upon installation (see below).</li></ol>
</li>
</ol>

### The codes currently in ComboCode
* <b> GASTRoNOoM:</b> L. Decin (KU Leuven, Belgium) 
* <b> MCMax:</b> M. Min (UvA, the Netherlands)

## Requirements
Currently the code has been tested to run on Unix-based systems, more specifically Fedora and Mac OS X. In principle, any operating systems that fulfills the requirements listed below should be able to run ComboCode. The code runs ons machines with an internal memory of 8 GB, but less is likely fine as well. The memory requirements are primarily set by the numerical codes included in ComboCode.

First and foremost, you require a Python 2.7 (not Python 3!) distribution installed on your machine. I recommend Anaconda, which allows for very flexible package management. Specific packages required to be installed in your python distribution are: PyPDF2, h5py, ephem, astropy and lmfit. (e.g., in case of Anaconda, run "pip install ephem" in the shell after Anaconda installation)

Secondly, we are using git for the version control of the repository. Follow the instructions given <a href="https://help.github.com/articles/set-up-git/"> here</a> to set up git on your local machine, after you have created a user account at GitHub.

Thirdly, you require a recent installation of the gfortran compiler. 

Fourthly, you need the IvS repository, which can be installed as described in the readme.txt file located <a href="https://github.com/robinlombaert/IvSPythonRepository"> here</a>. With the IvS repository comes data that are used by the package. The bare minimum required to run ComboCode is available <a href="http://ster.kuleuven.be/~robinl/cc/ivsdata.tar.gz"> here</a>. Make sure to update the config.py file with the location of the unpacked ivsdata upon installation of the IvS repository. 

Lastly, ComboCode can be used to its fullest potential when working in tandem with the radiative-transfer codes GASTRoNOoM and MCMax. To use these codes, permission is required from the owners listed above. You can contact them directly, or through me (robinlombaert on GitHub). Once installed, make sure the executables of each code are linked in your Bin folder, and you will be able to run the codes through ComboCode.

## Installation
Once the requirements are sorted out, you can get to work with ComboCode. You can download the code to your hard disk right away (or if you intend to submit code for this repository: fork the repo, and clone that -- see Developer's Manual below): 
* Clone the git repository to create a local copy, located in ~/ComboCode/ (or whichever location your prefer, replace ~/ with the parent folder you want):
    - $ cd ~/
    - $ git clone https://github.com/IvS-KULeuven/ComboCode.git ComboCode

* Copy the contents of the usr.dist/ folder to the usr/ folder. The usr/ folder contains user-specific settings that are not updated with the repository. Any changes made to the structure of those files will come through usr.dist/ and won't affect the usr/ folder, allowing the user to save personal settings before updating usr/ files. Note especially usr/Path.dat which contains all the relevant folders for ComboCode.
    - $ cd ~/ComboCode/
    - $ mkdir usr/
    - $ cp usr.dist/* usr/.

* Add the ComboCode home folder you chose to the PYTHON\_PATH in your ~/.bash_profile.

* Updating your own clone of ComboCode to the most recent version can be done with:
    - $ cd ~/ComboCode/
    - $ git pull

* Lastly, when running ComboCode in conjunction with the radiative-transfer codes listed above, ComboCode will write the output to folders for GASTRoNOoM and MCMax separately. You can choose these locations by adding them to usr/Path.dat.

## Documentation
Up-to-date documentation that goes with the package is available on GitHub at:

<a href="https://IvS-KULeuven.github.io/ComboCode"> ComboCode Documentation</a>

## Developer's Manual
If you want to make changes to ComboCode, you should fork the repository to your own github account. The convention is to call your online version of the code *origin*. It is this version that your local machine will refer to when pulling and pushing changes. Once downloaded on a local machine, you can create your own developer's branch that will not interfere with the origin/master branch. In general, it is advised NEVER to work in the origin/master branch. Keep your origin/master branch up-to-date with the original repository (called upstream below), but don't meddle with it, and never merge your changes into the origin/master branch. 

### Setting up your developer's environment for ComboCode
* Go to the main IvS-KULeuven/ComboCode.git page and fork the repository to your account. 

* Clone a copy of the code in your account to your local machine (fill in your github user name, and change the ComboCode folder to whatever you want for your copy. Replace ~/ with the parent folder you want.)
    - $ cd ~/
    - $ git clone https://github.com/YOUR_USERNAME/ComboCode.git ComboCode

* Tell git what the original "upstream" version of the code is at IvS-KULeuven (i.e. the original online version of the code -- don't confuse with origin, i.e. *your* online version of the code):
    - $ cd ~/ComboCode/
    - $ git remote add --track master upstream https://github.com/IvS-KULeuven/ComboCode.git

* Copy the contents of the usr.dist/ folder to the usr/ folder. Note especially usr/Path.dat which contains all the relevant folders for ComboCode.
    - $ cd ~/ComboCode/
    - $ mkdir usr/
    - $ cp usr.dist/* usr/.

* Add the ComboCode home folder to the PYTHON\_PATH in your ~/.bash_profile.

* Tell git that any subfolders in the ComboCode home folder other than cc/, aux/ and usr.dist/ are to be ignored (including, e.g., usr/ or input/). For this, you can create (or update) the file ~/ComboCode/.git/info/exclude with paths/files to be excluded from any git tracking. An example of such a file for ComboCode is available <a href="http://ster.kuleuven.be/~robinl/cc/exclude"> here</a>. You can update this as you go if more folders or files end up in ~/ComboCode/ that you don't want tracked.

### Keeping your code up-to-date with the upstream version
* Now you can keep your master branch up-to-date with the upstream. As long as you never make any changes in your master branch, you can merge upstream/master with origin/master (make sure to be on the master branch):
    - $ git fetch upstream
    - $ git checkout master
    - $ git merge upstream/master

* Then update your GitHub repo with the changes:
    - $ git push origin master

### Making your own changes to the code
When developing code (such as hotfixing, or adding any other changes), always work in different branch, i.e. one that is not the master branch. 
* Creating a developer's branch can be done by (it is named dev here, but it can be named whatever you want):
    - $ cd ~/ComboCode/
    - $ git branch dev
    - $ git checkout dev

* Once you've added changes to the code, and you want to suggest those changes to be included in ComboCode, you can commit the changes and push to the origin of your account on GitHub:
    - $ git add changedFile.py
    - $ git commit -m "Document your changes! Either in-line, or without the -m tag so a text editor opens for documentation."

* If you don't want to add every file you've changed separately, you can commit all changed files that are being tracked by git with one command:
    - $ git commit -a -m "Document your changes! Either in-line, or without the -m tag so a text editor opens for documentation."

* Then you want to update your GitHub repo with the changes
    - $ git push origin  dev

* Then you have to go to your personal ComboCode GitHub page, and select the dev branch. Create a pull request to the upstream version of the code (i.e. the IvS-KULeuven account), and wait for the pull request to be merged with the upstream remote, as described <a href="https://help.github.com/articles/creating-a-pull-request/"> here</a>. Once the pull-request was accepted, you can proceed to update your master branch as described <a href="https://github.com/IvS-KULeuven/ComboCode#keeping-your-code-up-to-date-with-the-upstream-version">above</a>. 

* Finally, you want to make sure your developer's branch is updated with any changes to the upstream version by other users:
    - $ git merge master
