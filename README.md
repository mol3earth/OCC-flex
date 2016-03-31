# OCC-flex
My Masters Research routine

This is a number of matlab files that will take a bathymetric grid file, filter and analyze it in a number of ways, 
take profiles through it over features that you can choose yourself, or have the program choose for you. 
Then with these profiles, you choose features of interest to a modeling routine I adapted from Schouten et al., 2010.
You can hand pick or have the program analyze and choose features automatedley (caveat emptor here). 
When these profiles have data of interest, there are scripts to do inversions and compare them
to a suite of flexural models. This suite is constructed according to variables you define. 

For each feature a display of variable fits by sigma is plotted.
As well as the best fitting model is displayed over the original features topography.  

Provided is a structure that already contains profiles, and relevant information so that the inversion routine can be tested. 
