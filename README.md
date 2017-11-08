TiGraMa
=========

Copyright 2014-2017, All rights reserved

Technical University of Denmark, Kongens Lyngby, Denmark

Code written by A. Cereser and S. Schmidt

--------------------------

![Alt Text](https://media.giphy.com/media/vFKqnCdLPNOKc/giphy.gif)

--------------------------

Welcome to TiGraMa, a software suite to analyze **Ti**me-of-flight **Gr**ain **Ma**pping data collected using time-of-flight 3D neutron diffraction (ToF 3DND)<sup>[1](#myfootnote1),[2](#myfootnote2)</sup>.

All software was developed analyzing software collected at J-PARC, beamline ID06 (SENJU).

TiGraMa can analyze data collected using the two available ToF 3DND methodologies:

  - Methodology 1, with data collected by an imaging detector mounted in transmission, with high spatial and temporal resolution

  - Methodology 2, with data collected both by an imaging detector mounted in transmission and my diffraction detector banks

_In the current version, only dataset collected using Methodology 1 can be analyzed. The support to Methodology 2 will be added in a later version._

Practical considerations
------------------------

All the TiGraMa scripts were developed during proof-of-principle data analyses. As such, the code is _not_ optimized for speed. To minimize computing time, it is recommended to run the TiGraMa scripts on a high performance computer.

Reconstruction steps
--------------------

0. Preprocess the collected frames. The aim of the steps is to reduce the noise and improve the signal.

   0.1. Perform ![overlap correction](http://iopscience.iop.org/article/10.1088/1748-0221/9/05/C05026/meta) using a dead-time correction algorithm.

   The code to perform the overlap correction was originally developed by A. Tremsin for the Windows platform, and has been ported to Unix by P.M. Larsen. To run it on a Mac, you need to

    * Install ![Homebrew](http://stacks.iop.org/1748-0221/9/i=05/a=C05026?key=crossref.88229b2d88c5ffd1bc62280555bdb4a1) or another package manager

    * On a terminal, launch `brew install cfitsio`

    * From the `Overlap_correction/linux_correction` folder, compile the code using `make clean && make`

    * Run `./fits_correction source_folder/example_file.fits destination_folder`

   To use the code from a Linux computer, it should be enough to change the path to `cfitsio` in the `Makefile` and in the `.h` files, and `clang++` to `g++`

   0.2. Calculate, from the open beam collected before and after the measurements, the open beam for the considered projection. Script: `Estimate_OB_no_median.m`

   0.3. Divide the collected dataset by the corresponding open beam. Script: `OB_correction.m`, which calls `OB_correction_function.m`

   0.4. Normalise using a rolling median. Code: `Correction_background.m`, which calls `Correction_background_function.m`

   Step 0.1 should run on a laptop (with data on an external drive) before uploading data to a server, where the more computationally consuming steps (0.2 to 0.4) run using Matlab.

1. Filter the frames using a signal enhancement filter, such as Murofi <sup>[1](#myfootnote1)</sup> (not included in TiGraMa)


References
----------

<a name="myfootnote1">1</a>: A. Cereser, M. Strobl, S. A. Hall, A. Steuwer, R. Kiyanagi, A. S. Tremsin, E. Bergbäck Knudsen, T. Shinohara, P. K. Willendrup, A. Bastos da Silva Fanta, S. Iyengar, P. M. Larsen, T. Hanashima, T. Moyoshi, P. M. Kadletz, P. Krooß, T. Niendorf, M. Sales, W. W. Schmahl & S. Schmidt (2017). “Time-of-Flight Three Dimensional Neutron Diffraction in Transmission Mode for Mapping Crystal Grain Structures.” Scientific Reports 7 (1): 9561, 9561. doi:[10.1038/s41598-017-09717-w](https://www.nature.com/articles/s41598-017-09717-w).

<a name="myfootnote2">2</a>: A. Cereser, S. A. Hall, A. Steuwer, M. Strobl & S. Schmidt (2016). [Time-of-flight 3D Neutron Diffraction for Multigrain Crystallography](http://findit.dtu.dk/en/catalog/2349663834). Department of Physics, Technical University of Denmark.

License
-------

This software is covered by the GNU General Public License.


[1]: http://stacks.iop.org/1748-0221/9/i=05/a=C05026?key=crossref.88229b2d88c5ffd1bc62280555bdb4a1
