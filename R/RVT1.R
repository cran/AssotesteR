RVT1 <-
function(y, X, maf=0.05, perm=100)
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
    if (mode(maf) != "numeric" || length(maf) != 1
            || maf <= 0  || maf >= 1)
        stop("argument 'maf' must be a value between 0 and 1")
    if (mode(perm) != "numeric" || length(perm) != 1
            || perm < 0 || (perm %% 1) !=0)
    {
        warning("argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
    }
	
    ## get minor allele frequencies
    MAFs = colSums(X, na.rm=TRUE) / (2*nrow(X))
	## are there any rare variants?
    rare = sum(MAFs < maf) 
	if (rare == 0)
	    stop(paste("\n", "Oops: No rare variants below maf=", 
		      maf, " were detected. Try a larger maf", sep=""))
	## get only rare variants
	X.new = X[ , MAFs < maf]	
	## convert to 0 and 1 (0: no rare copies;  1: rare copies)
	X.new[X.new == 2] = 1 
	## proportion of rare variants
	x.prop = rowMeans(X.new, na.rm=TRUE)
    # center phenotype y
    y.new = y - mean(y)
    # get score vector U
    U = sum(y.new * x.prop)
	## V
	xv = sum((x.prop - mean(x.prop))^2)
	V = mean(y) * (1 - mean(y)) * xv
	## get score
	score = sum(U^2 / V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
        score = 0
    asym.pval = 1 - pchisq(score, 1)
	
	## permutations
	perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			# center phenotype y
			y.perm = y[perm.sample] - mean(y[perm.sample])
			# get score vector
			U.perm = sum(y.perm * x.prop)
			x.perm[i] = sum(U.perm^2 / V)
		}
		# permuted p-value 
		perm.pval = sum(x.perm > score) / perm	
	}
	
	## results
	name = "RVT1: Rare Variant Test 1"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, perm)
	arg.spec = as.character(arg.spec)
	names(arg.spec) = c("controls", "cases", "variants", "rarevar", "maf", "n.perms")	
    res = list(rvt1.stat=score, asym.pval=asym.pval, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

