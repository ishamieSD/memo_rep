# create a distance matrix based on OBE-DTW (subsequence matching), with input matrix (channel-by-time)
#
# Created by Xi Jiang, Sept.1st, 2014
# Dependency: package "dtw" for dynamic time warping

open_matrix_dtw <- function (m){ # m is the input matrix, with every row representing a time segment

	library(dtw)
	dims <- dim(m)
	
	D <- NULL
	for (i in seq(1:dims[1])){
		query <- m[i,]
		D_row <- NULL

		for (j in seq(1:dims[1])){
			reference <- m[j,]
			alignment <- dtw(query,reference,step.pattern=asymmetric,open.begin=TRUE,open.end=TRUE)
			D_row <- c(D_row,alignment$distance)
		}
		D <- rbind(D,D_row)
		print(paste0("Current row: ", i))
	}
	return(D)
}
