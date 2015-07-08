# Welcome to the ComboCode README file
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

## Requirements
Currently the code has been tested to run on Unix-based systems, more specifically Fedora and Mac OS X. In principle, any operating systems that fulfills the requirements listed below should be able to run ComboCode. The code runs ons machines with an internal memory of 8 GB, but less is likely fine as well. The memory requirements are primarily set by the numerical codes included in ComboCode.

First and foremost, you require a Python distribution installed on your machine. I recommend Anaconda, which allows for very flexible package management. If you have the Anaconda distribution installed, you only require the pyfits package to be installed in addition. 

Secondly, we are using got for the version control of the repository. Follow the instructions given <a href="https://help.github.com/articles/set-up-git/"> here</a> to set up git on your local machine, after you have created a user account for git.

Thirdly, you need the IvS repository, which can be installed as described <a href="https://github.com/JorisDeRidder/IvSPythonRepository"> here</a>. Make sure the IvS repository is included in your PYTHON_PATH. 

Lastly, ComboCode can be used to its fullest potential when working in tandem with the radiative-transfer codes GASTRoNOoM and MCMax. Contact the owners above for the respective codes, (and/)or contact me (robinlombaert on GitHub) to gain access to the source of those codes. 

## Installation
Once the requirements are sorted out, you can get to work with ComboCode. You can download the code to your hard disk right away (or fork the repo, and clone that -- see Developer's Manual below): 
* Clone the git repository to create a local copy, located in ~/ComboCode/:
    - $ cd 
    - $ git clone https://github.com/IvS-KULeuven/ComboCode.git ComboCode

* Updating your own clone of ComboCode to the most recent version can be done with:
    - $ cd ~/ComboCode/
    - $ git pull

Note that this does not include the ~/ComboCode/Data/ files. You can copy those over from <a href="http://ster.kuleuven.be/~robinl/cc/Data/"> http://ster.kuleuven.be/~robinl/cc/Data/</a>. Afterwards, you can change those files at your whim. They do not come as part of the repository.


In the future, a shell script for the installation of ComboCode will be provided. It will include: 
* Anaconda python distribution installation, including additional packages required to run ComboCode
* Installation of the IvS repository
* Downloading additional Data files not part of the repository but required to run ComboCode
* Creating a master branch in your local git repository that can be used to run the code

## Documentation
Up-to-date documentation that goes with the package is available on GitHub at:

<a href="https://IvS-KULeuven.github.io/ComboCode"> ComboCode Documentation</a>

## Developer's Manual
If you want to make changes to ComboCode, you should fork the repository to your own github account. The convention is to call your online version of the code *origin*. It is is this version that your local machine will refer to when pulling and pushing changes. Once downloaded on a local machine, you can create your own developer's branch that will not interfere with the origin/master branch. In general, it is advised NEVER to work in the origin/master branch. Keep your origin/master branch up-to-date with the original repository (called upstream below), but don't meddle with it, and never merge your changes into the origin/master branch. 

### Setting up your developer's environment for ComboCode
* Go to the main IvS-KULeuven/ComboCode.git page and fork the repository to your account. 

* Clone a copy of the code in your account to your local machine (fill in your github username)
    - $ cd 
    - $ git clone https://github.com/<YOUR_USERNAME>/ComboCode.git ComboCode

* Tell git what the original "upstream" version of the code is at IvS-KULeuven (i.e. the original online version of the code -- don't confuse with origin, or *your* online version of the code):
    - $ git remote add --track master upstream https://github.com/IvS-KULeuven/ComboCode.git

* Tell git that any subfolders in ~/ComboCode/ other than cc/ and aux/ are to be ignored. For this, you can create (or update) the file ~/ComboCode/.git/info/exclude with paths/files to be excluded from any git tracking. An example of such a file for ComboCode is available <a href="http://ster.kuleuven.be/~robinl/cc/exclude"> here</a>.

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
