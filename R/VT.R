VT <-
function(y, X, maf=0.05, perm=100)
{
    ## checking argumetns
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
    if (mode(maf) != "numeric" || length(maf) != 1
            || maf <= 0  || maf >= 1)
        stop("argument 'maf' must be a value between 0 and 1")
    if (mode(perm) != "numeric" || length(perm) != 1
            || perm < 0 || (perm %% 1) !=0)
    {
        warning("argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
    }
    
    ## running vt method
    mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
    h.maf = sort(unique(mafs))
    vt.stat = .vt.method(y, X, mafs, h.maf)
	
    ## permutations
	perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			x.perm[i] = .vt.method(y[perm.sample], X, mafs, h.maf)
		}
		## p-value
		perm.pval = sum(x.perm > vt.stat) / perm
	}
	
	## results
	name = "VT: Variable Threshold"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), maf, perm)
	names(arg.spec) = c("controls", "cases", "variants", "maf", "n.perms")	
    res = list(vt.stat=vt.stat, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

