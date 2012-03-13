SEQSUM <- 
function(y, X, perm=100)
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
    if (mode(perm) != "numeric" || length(perm) != 1
            || perm < 0 || (perm %% 1) !=0)
    {
        warning("argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
    }
    
	## vector for storing stastitic and signs
	seqsum.stat = rep(0, ncol(X))
	signs = rep(-1, ncol(X))

	## start with first variant
	X.new = X
	## positive effect
    seqsum.stat = .uni.score(y, rowSums(X.new, na.rm=TRUE))	
	## negative effect
	X.new[,1] = (-1) * X[,1]
    aux.neg = .uni.score(y, rowSums(X.new, na.rm=TRUE))
    ## who is the best
	if (seqsum.stat > aux.neg)
	{
	    X.new[,1] = X[,1]
	    signs[1] = 1	
	} else seqsum.stat = aux.neg

	## continue with the other variants
	for (j in 2:ncol(X))
	{
		## negative effect model
		X.new[,j] = (-1) * X[,j]
		aux.neg = .uni.score(y, rowSums(X.new, na.rm=TRUE))
		# who is the best?
		if (seqsum.stat > aux.neg)
		{
		    X.new[,j] = X[,j]
			signs[j] = 1
		} else seqsum.stat = aux.neg
	}

	## permutations
	perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			y.perm = y[perm.sample]
			## start with first variant
			X.new = X
			perm.seqsum = .uni.score(y.perm, rowSums(X.new, na.rm=TRUE))
			## negative effect
			X.new[,1] = (-1) * X[,1]
			aux.neg = .uni.score(y.perm, rowSums(X.new, na.rm=TRUE))
			if (perm.seqsum > aux.neg) 
			{
				X.new[,1] = X[,1]
			} else perm.seqsum = aux.neg
			## continue with the other variants
			for (j in 2:ncol(X))
			{
				X.new[,j] = (-1) * X[,j]
				aux.neg = .uni.score(y.perm, rowSums(X.new, na.rm=TRUE))
				if (perm.seqsum > aux.neg) 
				{
					X.new[,j] = X[,j]
				} else perm.seqsum = aux.neg
			}
			x.perm[i] = perm.seqsum
		}
		# p-value 
		perm.pval = sum(x.perm > seqsum.stat) / perm	
	}
	
	## results
	name = "SEQSUM: Sequential Sum Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
	names(arg.spec) = c("controls", "cases", "variants", "n.perms")	
    res = list(seqsum.stat=seqsum.stat, perm.pval=perm.pval, signs=signs, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}
