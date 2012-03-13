CMAT <- 
function(y, X, maf=NULL, weights=NULL, perm=100)
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
	if (!is.null(maf))
	{
	    if (mode(maf) != "numeric" || length(maf) != 1
    	    || maf <= 0  || maf >= 1)
    		stop("argument 'maf' must be a value between 0 and 1")
	}
    if (!is.null(weights))
	{ 
	    if (mode(weights) != "numeric" || !all(weights >= 0)) 
		    stop("argument 'weights' must contain non-negative numbers")
        if (length(weights) != ncol(X)) 
	        stop("length of 'weights' differs from number of columns in 'X'")
    } else {
	    weights = rep(1, ncol(X))
	}
    if (mode(perm) != "numeric" || length(perm) != 1
         || perm < 0 || (perm %% 1) !=0) 
    {
          warning("argument 'perm' incorrectly defined. Value perm=100 is used")
  	 	  perm = 100
    }
    
    # get rare variants if maf != NULL
    rare = "NULL"
    if (!is.null(maf))
    {
    	## get minor allele frequencies
    	MAFs = colMeans(X, na.rm=TRUE) / 2    
		## are there any rare variants?
		rare = sum(MAFs < maf)
		if (rare > 0) 
		{
			X.new = X[ , MAFs < maf]
			weights = weights[MAFs < maf]
			rare = sum(weights != 0)
			if (sum(weights) == 0)
			    stop(paste("Oops!, with maf =", maf, ", all weights are zero"))
		}		
		if (rare == 0)
			stop(paste("Sorry, no rare variants detected below maf=", maf, sep=""))	
    } else {
    	X.new = X
    	maf = "NULL"
    }

	## apply cmat.method
	cmat.stat = .cmat.method(y, X.new, weights)

	## permutations
    perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			x.perm[i] = .cmat.method(y[perm.sample], X.new, weights) 
		}
		## p-value 
		perm.pval = sum(x.perm > cmat.stat) / perm
	}

    ## results
	name = "CMAT: Cumulative Minor Allele Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, perm)
	names(arg.spec) = c("controls", "cases", "variants", "rarevar", "maf", "n.perms")	
    res = list(cmat.stat=cmat.stat, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
    return(res)    	  
}

