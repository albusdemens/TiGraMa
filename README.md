TiGraMa
=========

Copyright 2014-2017, All rights reserved

Technical University of Denmark, Kongens Lyngby, Denmark

Code written by A. Cereser and S. Schmidt

--------------------------

![Alt Text](https://github.com/albusdemens/TiGraMa/blob/master/anim.gif)

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

The reconstruction steps are described in details in <sup>[2](#myfootnote2)</sup>. The Matlab scripts run on versions R2015a - R2017a. :exclamation: Before running a script, change the folder paths :exclamation:

**Grain shape reconstruction**

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

2. Isolate the extinction spots with the desired area. Script: `2_Blobs_detector.m`, which calls `2_Blobs_detector_function.m`

*The next steps are computationally intensive. You should probably transfer your data to a HPC machine*

3. Divide the blobs by projection using `3_Divide_blobs_proj.sh` or something similar

4. For each projection, calculate the similarity between different blobs.

   * Scripts: `4_Blobs_similarity.m`, which calls `4_Blobs_similarity_function.m`  and `4_Shape_comparison_function.m`

   * Output: files in `Results_blobs_comp/`. There are two types of output: the files in `Map_blobs/` list the properties of the blobs, and the files in `List_minima/` list the parameters used to measure how similar two blobs are

5. Using the angular parameters in `List_minima/` combine similar blobs using `5_Track_lambda.m`, which calls `5_Track_lambda_function.m`. Output: `Lambda_IM_before_HT.txt`, listing the properties of the combined images, and the combined images, stored in `Comb_blobs/`

Approximately, the number of extinction spots per projection should be the same

6. Separate the combined images in folders (one per projection) using `6_commands_sim_diff_lambda.sh`

7. Export the information relative to the images in `Comb_blobs/` using `7_Extract_data_Comb_blobs.m`. Output: `Data_blobs.txt`

8. Plot the distribution of the number of extinction spots per projection. Script: `8_Plot_blobs_distr.m`

9. Clean the set of combined images by considering the centre of mass of the combined extinction spots. This is done by considering a series of cylinders, and for each cylinder combining the spots within it that have similar area (A_min > 0.5*A_max). Additional operations are performed to avoid repetitions

   * Scripts: `9_IM_comb_2.m`, which calls `9_Funct_im_comb_2.m` and `9_Funct_clean_blobs.m`

   * Output: `Data_blobs_final.txt` and `Image_combination_before_HT.txt`

10. Considering a rolling interval, group the extinction spots using the Hough transform.

   * Scripts: `10_Hough_trans_pol_rolling_final.m`, `10_Funct_HT_pol_rolling_final.m`, `10_Funct_HT_sec_Murofi_final.m`, `10_Funct_Hough_split_final.m`, `10_Min_max_fitting_function.m`. Structure:

   <center><img src=Structure_scripts_HT.png height=180 /></center>

   * At first, look at the point distribution using `10_Funct_HT_sec_Murofi_final.m` and tune the parameters. Then comment the figure option and run for all intervals

   * Output: `CM_alpha_R.txt`

11. Clean the results of the Hough transform, grouping the centre of mass values. In this way, we determine the centre of mass of the grains and group the corresponding combined extinction spots.

   * Scripts: `Point_grouping_to_3D.m`, which calls `Plot_and_group_CM_Murofi.m`, `Clean_result_HT.m` and `Comb_im_before_reconstr.m`

   * Output: `3D_reconstr_input.txt`, `Cleaned_results_HT.txt` and `CM_clusters_final.txt`

12. Reconstruct the 3D shape of the grains. Script: `12_Voxels_tagging_final_P.m`

13. Visualize the reconstructed volume using [Paraview](https://www.paraview.org/). For a quick introduction to how to visualize volumes using Paraview, see Sec. 1.4.4 of the [Recon3D manual](https://github.com/albusdemens/Recon3D/blob/master/Manual_Recon3D.pdf).

Scripts: `13_Plot_3D_voxels.m` and  `13_Slides_to_vtk.m`. If you are using the _matryoshka doll approach_, use `13_Count_grains.m` instead.

**Grain orientation**

14. For all reconstructed grains, get (omega, lambda) values of the corresponding extinction spots. Scripts: `14_OL_curves_final.m`, which calls `14_OL_curves_fnct.m`, and `14_Unify_OL_values.m`. Final output: `OL_final_cleared.txt`

15. Calculate the centre of mass of the grains. Script: `15_HT_orientation_sections.m`, which calls `15_F_count_points_cleaned.m`, `15_Funct_acc_matrix.m` and `15_Funct_fit_HT.m`. Output: `Final_grains_CM.txt`

16. Print in a file the properties of the isolated extinction spots images. Script: `16_Isolated_images_properties.m`. Output: `Properties_isolated_blobs.txt`

17. Clean the (omega, lambda) list, keeping only the values relative to images containing the CM of the corresponding grain. Code: `OL_values_grain_CM.m`, output: `OL_final_OK.txt`

18. Calculate the grain orientation from the distribution of (omega, lambda) values. Script: `Save_omega_lambda_fit.m`, calling `Plot_omega_lambda_single_grain.m` and `indexToF.m`

Contributions
-------------
Feel free to download TiGraMa and use it to analyze your dataset. If you spot a bug, send an email to alberto.cereser@gmail.com or flag it here on GitHub.

How to quote the code
---------------------
TiGraMa is a scientific software package and if you use it for your research work you should reference it. Please refer to the article presenting ToF 3DND<sup>[1](#myfootnote1)</sup>.

References
----------

<a name="myfootnote1">1</a>: A. Cereser, M. Strobl, S. A. Hall, A. Steuwer, R. Kiyanagi, A. S. Tremsin, E. Bergbäck Knudsen, T. Shinohara, P. K. Willendrup, A. Bastos da Silva Fanta, S. Iyengar, P. M. Larsen, T. Hanashima, T. Moyoshi, P. M. Kadletz, P. Krooß, T. Niendorf, M. Sales, W. W. Schmahl & S. Schmidt (2017). “Time-of-Flight Three Dimensional Neutron Diffraction in Transmission Mode for Mapping Crystal Grain Structures.” Scientific Reports 7 (1): 9561, 9561. doi:[10.1038/s41598-017-09717-w](https://www.nature.com/articles/s41598-017-09717-w).

<a name="myfootnote2">2</a>: A. Cereser, S. A. Hall, A. Steuwer, M. Strobl & S. Schmidt (2016). [Time-of-flight 3D Neutron Diffraction for Multigrain Crystallography](http://findit.dtu.dk/en/catalog/2349663834). Department of Physics, Technical University of Denmark.

License
-------

This software is covered by the GNU General Public License.


[1]: http://stacks.iop.org/1748-0221/9/i=05/a=C05026?key=crossref.88229b2d88c5ffd1bc62280555bdb4a1
