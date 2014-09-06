# downsample a matrix
#
# Modified by Xi Jiang, Sept.1st, 2014

downsample <- function (m, N){ # m is the input matrix, with every row representing a time segment, and keep every N sample

	dims <- dim(m)
	downsample <- NULL
	for (i in seq(1:dims[1])){
		v <- m[i,]
		seed <- c(TRUE,rep(FALSE,N-1))
		cont <- rep(seed,ceiling(length(v)/N))[1:length(v)]
		downsample <- rbind(downsample,v[which(cont)])
	}
	return(downsample)
}
