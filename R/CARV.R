CARV <- 
function(y, X, waf=FALSE, signs=FALSE, approach="hard", maf=0.05, perm=100)
{
    ## checking arguments
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
    if (!is.logical(waf))
         stop("argument 'waf' must be TRUE or FALSE")
    if (!is.logical(signs))
         stop("argument 'signs' must be TRUE or FALSE") 
    if (!(approach %in% c("stepup", "variable", "hard")))
	    stop("argument 'approach' must be one of 'stepup', 'variable' or 'hard'")
	if (mode(maf) != "numeric" || length(maf) != 1 || maf <= 0 || maf > 1)
        maf = 0.05
    if (mode(perm) != "numeric" || length(perm) != 1 || perm < 0 || (perm %% 1) !=0)
    {
        warning("argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
    }
	
	## how many variants
	M = ncol(X)

	## weights 'ak' to incorporate allele frequencies
    ak = rep(1, M)
	if (waf) 
		ak = .weights.wss(y, X)   # weights a la madsen & browning
	
	## signs 'sk' of the variant effect 
	sk = rep(1, M)
	if (signs)
	{
	    nAs = apply(X[y==1,], 2, function(x) sum(!is.na(x)))
	    nUs = apply(X[y==0,], 2, function(x) sum(!is.na(x)))	 
		maf.A = colSums(X[y==1,], na.rm=TRUE) / (2*nAs)
		maf.U = colSums(X[y==0,], na.rm=TRUE) / (2*nUs)
		# more prevalent in cases than controls
		sk[maf.U > maf.A] = -1 
	}

	## prepare data for approaches
	ymean = mean(y)
    y.cen = y - ymean
    X.cen = scale(X, scale=FALSE)
    w = ak * sk	
	
    ## hard approach
    if (approach == "hard") 
	{
	    # get minor allele frequencies
		MAFs = colMeans(X, na.rm=TRUE) / 2
    	# are there any rare variants?
	    if (sum(MAFs < maf) == 0)
	        stop(paste("\n", "Ooops: No rare variants below maf=", 
				 maf, " were detected. Try a larger maf", sep=""))
		carv.stat = .hard.approach(y.cen, X.cen, w, MAFs, maf)
	} 
    ## variable approach
	if (approach == "variable")
	{
	    # get minor allele frequencies
		MAFs = colMeans(X, na.rm=TRUE) / 2
        several.maf = sort(unique(MAFs))
        carv.stat = .variable.approach(y.cen, X.cen, w, MAFs, several.maf)
	}
	## step-up approach
	if (approach == "stepup")
	{ 
	    carv.stat = .stepup.approach(y.cen, X.cen, w)
	}

	## permutations
	perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			# center phenotype y
			y.perm = y[perm.sample] - ymean
			# get score vector
			x.perm[i] = switch(approach, 
				"hard" = .hard.approach(y.perm, X.cen, w, MAFs, maf),
				"variable" = .variable.approach(y.perm, X.cen, w, MAFs, several.maf),
				"stepup" = .stepup.approach(y.perm, X.cen, w)
			)
		}
		# p-value 
		perm.pval = sum(x.perm >= carv.stat) / perm	
	}

	## results
	if (waf) mywaf="TRUE" else mywaf="FALSE"
	if (signs) mysigns="TRUE" else mysigns="FALSE"
	name = "CARV: Comprehensive Approach to Analyzing Rare Genetic Variants"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, mywaf, mysigns, approach, maf)
	names(arg.spec) = c("controls", "cases", "variants", "n.perms", "waf", "signs", "approach", "maf")	
    res = list(carv.stat=carv.stat, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

