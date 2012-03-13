CMC <-
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
    if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>=1)
        stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
    if (mode(perm) != "numeric" || length(perm) != 1
		    || perm < 0 || (perm %% 1) !=0) 
	{
	    warning("argument 'perm' incorrectly defined. Value perm=100 is used")
	    perm = 100
	}	

	## number of individuals N
    N = nrow(X)
    ## get minor allele frequencies
    MAF = colMeans(X, na.rm=TRUE) / 2   
    ## how many variants < maf
    rare.maf = MAF < maf
    rare = sum(rare.maf)
    ## collapsing
    if (rare <= 1) 
    {   
	    # if rare variants <= 1, then NO collapse is needed
        X.new = X
    } else {
        # collapsing rare variants into one column
        X.collaps = rowSums(X[,rare.maf], na.rm=TRUE)
        X.collaps[X.collaps != 0] = 1
        # joining collapsed to common variants
        X.new = cbind(X[,!rare.maf], X.collaps)	   
    }
	## change values to -1, 0, 1
	X.new = X.new - 1
   	## number of new variants
	M = ncol(X.new)
    ## Hotellings T2 statistic
    cmc.stat = .cmc.method(y, X.new)

    ## Asymptotic p-values
    # under the null hypothesis T2 follows an F distribution 
    f.stat = cmc.stat * (N-M-1)/(M*(N-2))
    df1 = M          # degrees of freedom  
    df2 = N - M - 1  # degrees of freedom  
    asym.pval = 1 - pf(f.stat, df1, df2)

    ## under the alternative hyposthesis T2 follows a chi-square distr
    # pval = 1 - pchisq(cmc.stat, df=M)
	    
	## permutations
	perm.pval = NA
	if (perm > 0)
	{
		x.perm = rep(0, perm)
		for (i in 1:perm)
		{
			perm.sample = sample(1:length(y))
			x.perm[i] = .cmc.method(y[perm.sample], X.new) 
		}
		# p-value 
		perm.pval = sum(x.perm > cmc.stat) / perm
	}
		
	## results
	name = "CMC: Combined Multivariate and Collapsing Method"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare,  maf, perm)
	arg.spec = as.character(arg.spec)
	names(arg.spec) = c("controls", "cases", "variants", "rarevar", "maf", "perm")	
    res = list(cmc.stat=cmc.stat, asym.pval=asym.pval, perm.pval=perm.pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)    
}

