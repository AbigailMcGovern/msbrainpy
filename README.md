# msbrainpy
A Python package relating to mouse brains (i.e., the 'ms' and 'brain' parts) and Abigail's PhD (vaguely enough). Clearly taking the descriptions very seriously... 


## Contents
CURRENT:
* brainseries - functions for processing series of sparse brain slices. Specifically, in situ hybridisation images. 
* atlas - dealing with atlas-related data
  - wrangle - manipulating data frames containing atlas related data and managing data in the file system.
  - visualise - plotting tools specific to anatomical data.
  - query - query the Allen Brain Institute's RMA API to obtain atlas-related data. Supports obtaining gridded gene expression data, Structure Unionize gene expression data, and raw insitu images. Plan to access GenePaint insitu images in the future.

OLD - LIGHTSHEET:
* io - reading and writing image data
* base - a number of fundemental functions that didn't fit elsewhere
* chain - several classes and functions for the serial processing of data using combinations of functions and specific arguments
* teraPrep - writing and nameing files according to TeraStitcher input requirments (for stitching lightsheet output)
* quantify - quantification of features (so far nuclear detection) in stitched light sheet data
  - processing - smaller fundemental processing functions 
  - sample - applying processing pipelines to a directory of images sampled from a stitched lightsheet data set
  - volume - applying processing piplines to larger volumes
* map - mapping quantifications onto an atlas

