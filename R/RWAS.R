RWAS <-
function(y, X, perm=NULL)
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
#    if (!is.null(weights))
#	{ 
#	    if (mode(weights) != "numeric" || !all(weights >= 0)) 
#		    stop("argument 'weights' must contain non-negative numbers")
#       if (length(weights) != ncol(X)) 
#	        stop("length of 'weights' differs from number of columns in 'X'")
#    } else {
#	    weights = rep(1, ncol(X))
#	}
	if (!is.null(perm))
	{
	    if (mode(perm) != "numeric" || length(perm) != 1
		    || perm < 0 || (perm %% 1) !=0) 
		{
	        warning("Argument 'perm' incorrectly defined. Value perm=100 is used")
			perm = 100
		}
	} else perm=0
 
	## running rwas
    weights = rep(1, ncol(X))
    rwas.stat = .rwas.method(y, X, weights)	
    asym.pval = 1 - pnorm(rwas.stat)

	## permutations
	perm.pval = NA
	if (perm > 0)	     
	{
        x.perm = rep(0, perm)
        for (i in 1:perm)
        {
 	        perm.sample = sample(1:length(y))
  	        x.perm[i] = .rwas.method(y[perm.sample], X, weights) 
	    }
        # p-value 
        perm.pval = sum(x.perm > rwas.stat) / perm
	}
	
	## results
	name = "RWAS: Rare-Variant Weighted Aggregate Statistic"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
	names(arg.spec) = c("controls", "cases", "variants", "n.perms")	
	res = list(rwas.stat=rwas.stat, asym.pval=asym.pval, perm.pval=perm.pval, 
	      args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

