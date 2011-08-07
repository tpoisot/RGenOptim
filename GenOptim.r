x <- seq(from=-20,to=20,by=0.25)
target = function(x,pars) {
	with(as.list(pars),{
		return((sin(a*x+b)-cos(x+d))+c*x)
	})
}
opars <- list(a=0.4,b=-0.65,c=-0.4,d=8)
y <- target(x,opars)
yn <- target(x+rnorm(length(x),0,0.5),opars)

errorfunc = function(x,tomatch,basedat) sum(abs(tomatch-target(basedat,x)))/length(tomatch)

genparms	<- list(ngen=150,nind=200,recomb=TRUE,dorep=20)
seed		<- c(a=0,b=0,c=0,d=0)

mutate = function(lop)
{
	for(i in 1:length(lop))
	{
		lop[i] = lop[i] + rnorm(1,0,max(0.1,lop[i]/6))
	}
	return(lop)
}

recombin = function(g1,g2,n)
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
			fit <- apply(PopMat,1,function(x) errorfunc(x,goal,data))
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
		fit <- apply(PopMat,1,function(x) errorfunc(x,goal,data))
		output = list(set=as.list(PopMat[which.min(fit),]),pop=PopMat,fits=fit,fit.log=FitHist)
		return(output)
	})
}
		   
out <- GenOptim(x,yn,target,errorfunc,seed,genparms)

plot(x,yn,col='grey',pch=19)
lines(x,y)
lines(x,target(x,out$set),col='red')