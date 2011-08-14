## DESCRIPTION
## Genetic optimization method using recombination
##
## Timothee Poisot
## Universite Montpellier 2
##
## tpoisot@um2.fr
##

Rsq = function(x,y,pred) # R squared
{
	m <- mean(y)
	SStot <- sum((y-m)^2)
	SSerr <- sum((y-pred)^2)
	return(1-(SSerr/SStot))
}

mutate = function(lop) # Mutation
{
	for(i in 1:length(lop))
	{
		lop[i] = lop[i] + rnorm(1,0,max(0.1,lop[i]/6))
	}
	return(lop)
}

recombin = function(g1,g2,n) # Recombination
{
	tmat = matrix(0,ncol=length(g1),nrow=n)
	tmat[1,] <- g1
	tmat[2,] <- g2
	for(i in 3:nrow(tmat))
	{
		for(j in 1:length(g1))
		{
			tmat[i,j] <- ifelse(runif(1,0,1)<0.5,g1[j],g2[j])
		}
	}
	return(tmat)
}

GenOptim = function(data,goal,problem,error,seed,genparms)
{
	FitHist <- NULL
	InitFlag <- TRUE
	# Population matrix
	PopMat <- matrix(0,ncol=length(seed),nrow=genparms$nind)
	colnames(PopMat) <- names(seed)
	RPopMat <- PopMat
	# Fitness vector
	fitness <- numeric(length=genparms$nind)
	# Begin real code
	with(as.list(genparms),{
		# Initialize population by mutating the seed
		if(InitFlag)
		{
			for(ind in 1:nind) PopMat[ind,] <- mutate(seed)
			InitFlag <- FALSE
		}
		# Now go on with the search alogorithm
		for(gen in 1:ngen)
		{
			# Step 0 : create a new matrix
			tPopMat <- NULL
			# Step 1 : get the fitness of each element
			fit <- apply(PopMat,1,function(x) errorfunc(x,goal,data,problem))
			FitHist[gen] <- min(fit)
			# Step 2 : get the genotypes that reproduce
			OkGen <- PopMat[(rank(fit,ties='random')<=dorep),]
			# Step 3 : make random pairings
			OkGen <- OkGen[sample(c(1:dorep),replace=FALSE),]
			for(pairing in 1:(dorep/2))
			{
				recom <- recombin(OkGen[((2*pairing)-1),],OkGen[((2*pairing)),],2*(nind/dorep))
				recom <- t(apply(recom,1,mutate))
				tPopMat <- rbind(tPopMat,recom)
			}
			# Step 4 : clone the new matrix
			PopMat <- tPopMat
			colnames(PopMat) <- names(seed)
		}
		# Return the best result
		fit <- apply(PopMat,1,function(x) errorfunc(x,goal,data,problem))
		output = list(set=as.list(PopMat[which.min(fit),]),pop=PopMat,fits=fit,fit.log=FitHist)
		return(output)
	})
}

y.rnd = function(x,pars) {
	with(as.list(pars),{
		return(c*x^(-g))
	})
}

errorfunc = function(x,tomatch,basedat,rfunc) sum(abs(tomatch-rfunc(basedat,x)))/length(tomatch)
genparms	<- list(ngen=300,nind=100,recomb=TRUE,dorep=30)

seed.r		<- c(c=1,g=2.5)


# Uncomment to run
#GenOptim(x,y,y.rnd,errorfunc,seed.r,genparms)