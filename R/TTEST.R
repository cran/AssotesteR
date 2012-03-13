TTEST <-
function(y, X)
{
    ## checking argumetns
    if (!is.vector(y) || mode(y) != "numeric")
	    stop("argument 'y' must be a numeric vector")
	if (any(is.na(y))) 
	    stop("Sorry =(   No missing data allowed in argument 'y' ")	
    if (!all(y %in% c(0, 1)))
	    stop("Sorry =(  argument 'y' must contain only 0 and 1")
	if(!is.matrix(X) & !is.data.frame(X))
	    stop("argument 'X' must be a matrix or data.frame")    
	if (nrow(X) != length(y)) 
	    stop("'X' and 'y' have different lengths")
	if (!is.matrix(X)) X = as.matrix(X)

    ## number of individuals N, cases nA, controls nU
    N = nrow(X)
	M = ncol(X)
	nA = sum(y)
	nU = N - nA
    ## matrix of genotypes in cases
    Xx = X[y==1,] - 1  
    ## matrix of genotypes in controls  
    Yy = X[y==0,] - 1 
    ## get means
    Xx.mean = colMeans(Xx, na.rm=TRUE)
    Yy.mean = colMeans(Yy, na.rm=TRUE)
    ## center matrices Xx and Yy
    Dx = sweep(Xx, 2, Xx.mean)
    Dy = sweep(Yy, 2, Yy.mean)

    ## pooled covariance matrix
    if (sum(complete.cases(X)) == N)  # no missing values
    {  
        COV = (t(Dx) %*% Dx + t(Dy) %*% Dy) / (N-2)
	} else {  # with missing values
        ## covariance matrix of cases
        tDx = t(Dx)
        Sx = matrix(0, M, M)
        for (i in 1:nrow(tDx))
        {
            for (j in i:M)
            {
                Sx[i,j] = sum(tDx[i,] * Dx[,j], na.rm=TRUE)
	        }
        }
	    sx.diag = diag(Sx)
        Sx = Sx + t(Sx)
        diag(Sx) = sx.diag
        ## covariance matrix of controls
        tDy = t(Dy)
        Sy = matrix(0, M, M)
        for (i in 1:nrow(tDy))
        {
            for (j in i:M)
            {
                Sy[i,j] = sum(tDy[i,] * Dy[,j], na.rm=TRUE)
	        }
        }
	    sy.diag = diag(Sy)
        Sy = Sy + t(Sy)
        diag(Sy) = sy.diag
		## pooled covariance matrix
        COV = (1/(N-2)) * (Sx + Sy)	
	}

    ## general inverse
    if (nrow(COV) == 1) # only one variant
    { 
        if (COV < 1e-8) COV = 1e-8
        COV.inv = 1 / COV
    } else {
        COV.eigen = eigen(COV)
        eig.vals = COV.eigen$values  
        inv.vals = ifelse(abs(eig.vals) <= 1e-8, 0, 1/eig.vals)
        EV = solve(COV.eigen$vectors)
        COV.inv = t(EV) %*% diag(inv.vals) %*% EV
    }	
	
    ## Hotellings T2 statistic
    stat = t(Xx.mean - Yy.mean) %*% COV.inv %*% (Xx.mean - Yy.mean) * nA * nU / N
    T2.stat = as.numeric(stat)

    ## Asymptotic p-value
    ## under the null hypothesis T2 follows an F distribution 
    f.stat = T2.stat * (N-M-1)/(M*(N-2))
    df1 = M          # degrees of freedom  
    df2 = N - M - 1  # degrees of freedom  
    pval = 1 - pf(f.stat, df1, df2)
    
	## results
	name = "TTEST: Hotellings T2 Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X))
	names(arg.spec) = c("controls", "cases", "variants")	
    res = list(T2.stat=T2.stat, asym.pval=pval, args=arg.spec, name=name)
	class(res) = "assoctest"
	return(res)    
}

