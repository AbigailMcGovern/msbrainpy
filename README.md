# msbrainpy
Tools for working with light sheet microscopy data from whole mouse brains and for Allen Mouse Brain Atlas aligned data
from a variety of sources. 

## Modules and Subpackages
* io - reading and writing imageing data
* base - a number of fundemental functions that didn't fit elsewhere
* chain - several classes and functions for the serial processing of data using combinations of functions and specific arguments
* teraPrep - writing and nameing files according to TeraStitcher input requirments (for stitching lightsheet output)
* quantify - quantification of features (so far nuclear detection) in stitched light sheet data
  - processing - smaller fundemental processing functions 
  - sample - applying processing pipelines to a directory of images sampled from a stitched lightsheet data set
  - volume - applying processing piplines to larger volumes
* map - mapping quantifications onto an atlas
* atlas - dealing with atlas-related data
  - wrangle - manipulating data frames containing atlas data 
  - visualise - plotting tools specific to anatomical data
  - query - query the Allen Brain Institute's RMA API to obtain atlas-related data
