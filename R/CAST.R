CAST <- 
function(y, X, maf=0.05, test="fisher")
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
    if (!(test %in% c("fisher", "chisq")))
        stop("argument 'test' must be 'fisher' or 'chisq'")

    # number of individuals 'N' and variants 'M'
    N = nrow(X)
    ## get minor allele frequencies
    MAFs = colSums(X, na.rm=TRUE) / (2 * N)   
	## are there any rare variants?
    rare = sum(MAFs < maf)
	if (rare == 0)
	    stop(paste("\n", "Oops: No rare variants below maf=", 
		      maf, " were detected. Try a larger maf", sep=""))
	# collapsing genotypes across variants
	X.collapsed = rowSums(X[ , MAFs < maf], na.rm=TRUE)	
	## convert to 0 and 1 (0: no rare copies; 1: one or more rare copies)
	X.collapsed[X.collapsed != 0] = 1 

    # testing proportions between cases and controls
    props = table(y, X.collapsed)
	stat = switch(test, 
	   "fisher" = unlist(fisher.test(props)[c(3,1)]),
	   "chisq" = unlist(chisq.test(props)[c(1,3)])
	)
	names(stat) = NULL
	cast.stat = stat[1]
	pval = stat[2]

	## results
	name = "CAST: Cohort Allelic Sums Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, test)
	names(arg.spec) = c("controls", "cases", "variants", "rarevar", "maf", "test")	
    res = list(cast.stat=cast.stat, asym.pval=pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)
}

