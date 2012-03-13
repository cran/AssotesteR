CALPHA <-
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
	if (is.null(perm)) {
	    perm = 0
	} else {
  	    if (mode(perm) != "numeric" || length(perm) != 1
            || perm < 0 || (perm %% 1) !=0)
        {
            warning("argument 'perm' incorrectly defined. Value perm=100 is used")
            perm = 100
        }
	}
 
    nA = sum(y)
    nU = sum(y==0)
    p0 = nA / (nA + nU)
   
    m = ncol(X)
    # copies of the i-th variant type
    n = apply(X, 2, function(x) sum(x>0, na.rm=TRUE))
    # copies of the i-th variant type in the cases
    g = apply(X[y==1,], 2, function(x) sum(x>0, na.rm=TRUE))
   
    # Test statistic 
    calpha.stat = .calpha.method(y, X)
    # Variance of Talpha
    Valpha = 0
    for (i in 1:m) {
        for (u in 0:n[i]) {
            Valpha = Valpha + (((u - n[i]*p0)^2 - n[i]*p0*(1-p0))^2)*dbinom(u, n[i], p0)
	    }
    }
    names(Valpha) = NULL
    # Z score
    Zscore = calpha.stat / sqrt(Valpha)
   
    # asumptotic p-vaue
    if (Valpha==0) asym.pval=1 else
        asym.pval = 1 - pchisq(calpha.stat^2 / Valpha, df=1)
   
    # permutations
    perm.pval = NA
    if (perm != 0)
    {
        x.perm = rep(0, perm)
        for (i in 1:perm)
        {
 	        perm.sample = sample(1:length(y))
  	        x.perm[i] = .calpha.method(y[perm.sample], X) 
	    }
        # p-value 
        perm.pval = sum(x.perm^2 > calpha.stat^2) / perm
    }	  
      
	## results
    name = "CALPHA: c-alpha Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
    names(arg.spec) = c("controls", "cases", "variants", "n.perms")
    res = list(calpha.stat=calpha.stat, asym.pval=asym.pval, perm.pval=perm.pval, 
	    args=arg.spec, name=name)
    class(res) = "assoctest"	
    return(res)
}

