# SarcOrgTextureAnalysis
Sarcomere Organization Texture Analysis - code from Sutclife et al Scientific Reports 2018, PMID 29352247

Sarcomere Organization Texture Analysis

v1.0: published version

For licensing, please see LICENSE.txt

If you use this software in a scientific publication, please cite:

Matthew D. Sutcliffe, Philip M. Tan, Antonio Fernandez-Perez, Young-Jae Nam, Nikhil V. Munshi and Jeffrey J. Saucerman. (2018). High content analysis identifies unique morphological features of reprogrammed cardiomyocytes. Scientific Reports 8, 1258.

You should first run SOTA with the provided example images (see below).

Summary of main files
---------------------
master.m: loads images and runs analyses, generating figures of analyzed results
eCM_segment_function.m: performs cell segmentation
morph_texture_function.m: performs sarcomere organization analysis on segmented images
ideal.m: generates idealized images and then calls function for sarcomere analysis

Requirements
------------
MATLAB
MATLAB Image Processing Toolbox
MATLAB Signal Processing Toolbox
MATLAB Statistics and Machine Learning Toolbox

Example images
--------------
We recommend first running analysis with idealized images. The cell images used for the manuscript have been archived on figshare at:
https://doi.org/10.6084/m9.figshare.11390022
Place the included Images folder within the directory for SOTA and run master.m.
