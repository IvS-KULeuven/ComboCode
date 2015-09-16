# Overview: Conversion VIC3 to ThinKing
## Access
Connect to ThinKing through ssh at the address `username@login.hpc.kuleuven.be`
where the VSC username has to be filled in (e.g. vsc30226).

## Infrastructure
<ol>
<li>Improved login nodes: graphical interfaces possible</li>
<li>Two Ivy Bridge CPU's with 10 cores per node</li>
<li>Memory per node: 64 x 112 cores, 128 x 32 cores</li>
<li>Memory per core: ~ 3.2 GB, ~ 6.4 GB</li>
<li>If more memory is needed: Use Cerebro (sharing memory system on ThinKing)</li>
</ol>

## Code 
In principle all data that used to be on VIC3 should be available on ThinKing. 
However, ideally code is recompiled on the new server to improve performance
significantly. 
<ol>
<li>Tool chains are available, which are essentially libraries for packages to compile Python code, Fortran code, etc. </li>
<li>Currently available (April 2014): Intel (ifort), GNU  (gcc,gfortran,g77) </li>
<li>Python libraries and compilers at Python/2.7.6-goolf-2014a, which is the Python version number and the name of the tool chain. </li>
<li>Documentation through: module show toolchain/version </li>
<li>FFTW package should be available!</li>
<li>Access to other libraries is possible but not recommended.</li>
</ol>

## Queues
Do not use specified queues anymore! Use qsub or wsub (for worker queue system).

Names of the queues were changed from qreg, qlong, etc. to q1h, q24h, etc. to indicate queues appropriate for expected duration of calculation.

When queuing:
<ol>
<li>When queuing a job, specify, e.g., `#PBS -l nodes=1:ppn=20`</li>
<li>Bunch of jobs sent to a node. Perhaps pack per 20? You pay for the whole node in one go anyway!</li>
<li>Keep in mind memory limits. 3.2 GB per core should be OK.</li>
</ol>

When using the monitor tool: 
<ol>
<li>Load through: `module load monitor`</li>
<li>To check memory used, current jobs running, etc.</li>
<li>Works as before.
</ol>

More information available on the ThinKing documentation webpages!