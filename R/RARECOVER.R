RARECOVER <- 
function(y, X, maf=0.05, dif=0.5, perm=100)
{
    ## checking arguments
    if (!is.vector(y) || mode(y) != "numeric")
        stop("argument 'y' must be a numeric vector")
    if (any(is.na(y)))
        stop("Sorry =(   no missing data allowed in 'y'")
    if (!all(y %in% c(0, 1)))
        stop("Sorry =(   argument 'y' must contain only 0 and 1")
    if(!is.matrix(X) & !is.data.frame(X))
        stop("argument 'X' must be a matrix or data.frame")
    if (nrow(X) != length(y))
        stop("'X' and 'y' have different lengths")
    if (!is.matrix(X)) X = as.matrix(X)
    if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>1)
        stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
    if (mode(dif) != "numeric" || length(dif) != 1
            || dif <= 0  || dif >= 1)
        stop("argument 'dif' must be a value between 0 and 1")
    if (mode(perm) != "numeric" || length(perm) != 1
            || perm < 0 || (perm %% 1) !=0)
    {
        warning("argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
    }

    ## get minor allele frequencies
    MAF = colMeans(X, na.rm=TRUE) / 2   
    ## how many variants < maf
    rare = sum(MAF < maf)
	if (rare == 0)
	    stop(paste("Ooops!  No rare variants detected below maf=", maf, sep=""))
    ## collapsing
    if (rare == 1) 
    {   
	    # if rare variants <= 1, then NO collapse is needed
        X.new = X
    } else {
        X.new = X[,MAF < maf]	   
    }
    ## change genotype 2 into 1
    X.new[X.new == 2] = 1
    rarecov = .rarecov.method(y, X.new, dif)
    rc.stat = rarecov$stat
    chisq.pval = rarecov$pval
    rc.sel = rarecov$sel
    names(rc.sel) = NULL

    ## permutations
    perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			rarecov.perm = .rarecov.method(y[perm.sample], X.new, dif)
			x.perm[i] = rarecov.perm$stat
		}
		## p-value
		perm.pval = sum(x.perm > rc.stat) / perm
	}

    ## results
    name = "RARECOVER Algorithm"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, length(rc.sel), perm)
	arg.spec = as.character(arg.spec)
    names(arg.spec) = c("controls", "cases", "variants", "rarevar", "maf", "varsel", "n.perms")
	sel.names = names(X)[rc.sel] 
	if (is.null(sel.names)) 
	    sel.names = paste("var", rarecov$sel, sep="")
	names(rc.sel) = sel.names
    res = list(rc.stat=rc.stat, perm.pval=perm.pval, 
               set=rc.sel, args=arg.spec, name=name)
    class(res) = "assoctest"
    return(res)
}










