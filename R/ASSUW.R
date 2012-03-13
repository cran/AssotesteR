ASSUW <-
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
    ## get U and V
    getuv = .getUV(y, X)
    U = getuv$U
    V = getuv$V
    ## run score method
	stat.assuw = .assuw.method(U, V)
	assuw.stat1 = stat.assuw[1]  # stat normal
	p1.assuw = stat.assuw[2]     # pval normal

	## permutations
	perm.pval = NA
	if (perm > 0)
	{
        p1.perm = rep(0, perm)
		ymean = mean(y)
        for (i in 1:perm)
        {
 	        perm.sample = sample(1:length(y))
            # center phenotype y
		    y.perm = y[perm.sample] - ymean
            # get score vector
            U.perm = colSums(y.perm * X, na.rm=TRUE)
  	        perm.assuw = .assuw.method(U.perm, V)
		    p1.perm[i] = perm.assuw[2]
	    }
        # p-value 
        perm.pval = sum(p1.perm < p1.assuw) / perm  # normal
	}
	
	## results
	name = "ASSUW: Adaptive Weighted Sum of Squared Score"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
	names(arg.spec) = c("controls", "cases", "variants", "n.perms")	
    res = list(assuw.stat=assuw.stat1, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

