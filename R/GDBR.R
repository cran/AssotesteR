GDBR <- 
function(y, X, distance="IBS", weights=NULL, perm=NULL)
{
    ## checking argumetns
    if (!is.vector(y) || mode(y) != "numeric")
	    stop("argument 'y' must be a numeric vector")
	if (any(is.na(y))) 
	    stop("Sorry =(   No missing data allowed in argument 'y' ")	
    if (!all(y %in% c(0, 1)))
        stop("Sorry =(   argument 'y' must contain only 0 and 1")
    if(!is.matrix(X) & !is.data.frame(X))
        stop("argument 'X' must be a matrix or data.frame")
    if (nrow(X) != length(y))
        stop("'X' and 'y' have different lengths")
    if (!is.matrix(X)) X = as.matrix(X)
    if (!(distance %in% c('IBS', 'wIBS')))
	    stop(paste("\n", "Sorry =(   I don't recognize distance:", distance, "\n",
		     "Choose one from: 'IBS' or 'wIBS'")) 	
	if (distance == "wIBS" && is.null(weights))
	    stop("when using distance 'wIBS', you must provide a vector of 'weights'")
	if (!is.null(weights))
	{
	    if(!is.vector(weights) || mode(weights) != "numeric" || any(weights < 0))
		    stop("argument 'weights' must be a numeric vector with non-negative values")
	    if(length(weights) != ncol(X)) 
		    stop("length of 'weights' does not match number of columns in 'X'")
	}
	if (!is.null(perm))
	{
		if (mode(perm) != "numeric" || length(perm) != 1
		    || perm < 0 || (perm %% 1) !=0) 
		{
			warning("argument 'perm' incorrectly defined. Value perm=100 is used")
			perm = 100
		}
	} else perm = 0

	n = nrow(X)
	## Genomic similarity matrix
	Sim = switch(distance, 
	    "IBS" = gdbr_IBS(X),
		"wIBS" = gdbr_wIBS(X, weights)
	)
    ## distance matrix
    D = 1 - Sim
    A = (-1/2) * D^2
	I = diag(1, n)
	## association matrix
	A.cen = scale(A, scale=FALSE)
	G = A.cen - rowMeans(A.cen)
	
    # get Fstat ((very computer-intensive!!!)
    gdbr.stat = .gdbr.fstat(y, G)
	
	## permutation procedure
	perm.pval = NA
	if (perm > 0)
	{
		I = diag(1, length(y))
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			y.perm = y[perm.sample] - mean(y)
			# get projection 'hat' matrix
			H = y.perm %*% solve((t(y.perm) %*% y.perm)) %*% t(y.perm)
			# calculate F statistic
			Fstat.num = sum(diag(H %*% G %*% H))
			Fstat.denom = sum(diag((I - H) %*% G %*% t(I - H)))
			x.perm[i] = Fstat.num / Fstat.denom
		}
		# p-value 
		perm.pval = sum(x.perm > gdbr.stat) / perm
	}

	## results
	name = "GBDR: Genomic Distance-Based Regression"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, distance)
	names(arg.spec) = c("controls", "cases", "variants", "n.perms", "distance")	
    res = list(gdbr.stat=gdbr.stat, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}
