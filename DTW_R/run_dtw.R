# script for running DTW based on the functions open_matrix_dtw.R and downsample.R,
# assuming that, after loading and un-listing, the input to downsample.R is a matrix of time bins by sampling points.
# e.g. with 1000 time bins, 2-second bin length, and 512Hz sampling rate, the input would be a 1000-by-1024 matrix.
	
	library(dtw)
	library(R.matlab)    # needed for loading .mat data into R

	source('~/downsample.R')   # assuming the function files are under the default directory
	source('~/open_matrix_dtw.R')

	# after converting the loaded .mat file from a list to a matrix (or multiple matrices), 
	# assuming that said matrix is named "input":

	resampled_50Hz <- downsample(input,10)	  # downsample to 50Hz, if the bin size is 2 seconds
	distance_matrix <- open_matrix_dtw(resampled_50Hz)	# may take a few hours if the number of bins is ~2000

	# convert the matrix into a distance object, if the "hclust" function in R is to be used
	Dist <- as.dist(distance_matrix)
