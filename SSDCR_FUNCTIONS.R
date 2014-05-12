# Libraries

library(glasso)
library(matlab)

#===========================================================================
#===========================================================================
#===========================================================================

# Functions - for AIC, BIC, RISK and OTHER splits 

#===========================================================================
#===========================================================================
#===========================================================================
# This function prepares the data for the functions

datasetup = function(data,n.subj,T){

# data   = data set
# n.subj = number of subjects in the data set
# T      = number of time points in the data set

tempdata = as.matrix(data)
x        = c()

for (i in 1:n.subj){
	x = c(x,1:T)
}

row.names(tempdata) = x
data.1    = list()

for (i in 1:n.subj){
	data.1[[i]] = tempdata[(1+(i-1)*T):(T*i),]
}
return(data.1)
}

#===========================================================================
# The main function that calls all other functions

dcrfunction = function(x,datalist,n.subj,T,lower,upper,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type,perm_type,quantile_1,quantile_2){

# datalist    = data set to be analyzed
# n.subj      = number of subjects in the data set
# T           = number of time points in the data set
# lower       = starting point for the algorithm
# upper       = total number of time points in the experiment
# split.index = a list for output
# maxn        = the minimum distance between change points
# lambda.list = the list of regularization parameters
# n.rep       = number of replications in the bootstrap procedures
# block       = the block length for the stationary bootstrap in percentage terms
# refit.flag  = refit the glasso estimation procedure to reduce bias
# method_type = use BIC to find change points
# perm_type   = use the stationary bootstrap procedure for inference on the change points 
# quantile_1  = lower quartile for BIC reduction         
# quantile_2  = upper quartile for BIC reduction            

subwithsignsplits = c()

for (jlo in 1:n.subj){
	
	y = as.matrix(datalist[[jlo]])
	
	split.index = list()
	splits      = split.index.everything(x,y,lower,upper,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type)
	split.index = refit.split.perm(splits,perm_type,method_type,y,x,T,lambda.list,block,n.rep)

	
	if (length(split.index)>0){
	
		read.splits_actual    = convert.split(split.index$actual)
		read.splits_boot.perm = convert.split(split.index$boot.perm)
	
		write(read.splits_actual,file = sprintf("ActualSplitsSubject_%s",jlo),ncolumns = 4,sep = "\t",append=TRUE)
		write(read.splits_boot.perm,file = sprintf("BootSplitsSubject_%s",jlo),ncolumns = 4,sep = "\t",append=TRUE)

		boot_type                 = 1 # permutation

		results  = sign.splits.boot.edges(boot_type,method_type,y,T,lambda.list,split.index,quantile_1,quantile_2,block)
		
		n.edge   = length(results$corr.quants)
		nnn.edge = length(results$corr.quants[[1]])

		junk = matrix(NA,ncol=n.edge,nrow=nnn.edge)

		for (i in 1:n.edge){
			kkk = length(results$corr.quants[[i]])
			junk[1:kkk,i] = results$corr.quants[[i]][1:kkk]
		}

		write.table(junk,file = sprintf("EdgePermDistSubject_%s",jlo),sep = "\t",append=TRUE)

		# Save significant splits
		if (length(results$sign.splits) > 0){
	
			read.signif.splits.actual = convert.split(results$sign.splits)
			write(read.signif.splits.actual,file = sprintf("SignificantSplitsSubject_%s",jlo),ncolumns = 4,sep = "\t",append=TRUE)
			
			subwithsignsplits = c(subwithsignsplits,jlo)
		} else {
			print("No Significant Splits")
		}
		
		# Save partial correlation matrices
		n.corr.matrix = length(results$partial.corr.matrix)
		junk = c()
		for (i in 1:n.corr.matrix){
			junk = rbind(junk,results$partial.corr.matrix[[i]]) 
		}

		k.junk = dim(junk)[1]
		rownames(junk) = 1:k.junk
		write.table(junk,file = sprintf("PrecisionMatricesSubject_%s",jlo), sep = "\t",append=TRUE)

	} else {
		print("No Splits")

		# Make undirected graph for the entire data set
		results = sign.splits.boot.edges.iid(boot_type,method_type,y,T,lambda.list,block)
		
		n.edge   = length(results$corr.quants)
		nnn.edge = length(results$corr.quants[[1]])
		junk = matrix(NA,ncol=n.edge,nrow=nnn.edge)
		for (i in 1:n.edge){
			kkk = length(results$corr.quants[[i]])
			junk[1:kkk,i] = results$corr.quants[[i]][1:kkk]
		}
		write.table(junk,file = sprintf("EdgePermDistSubject_%s",jlo),sep = "\t",append=TRUE)
		# Save partial correlation matrices
		n.corr.matrix = length(results$partial.corr.matrix)
		junk = c()
		for (i in 1:n.corr.matrix){
			junk = rbind(junk,results$partial.corr.matrix[[i]]) 
		}

		k.junk = dim(junk)[1]
		rownames(junk) = 1:k.junk
		write.table(junk,file = sprintf("PrecisionMatricesSubject_%s",jlo), sep = "\t",append=TRUE)
		}
	}

write.table(subwithsignsplits,file = "SubjectswithSignificantSplits")

}

#===========================================================================
# This function finds splits the time course into partitions based on BIC 
# reductions it also calculates the undirected graph

split.index.everything = function(x,y,lower,upper,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type){

# x                 = time series of time
# y                 = data 
# lower             = lower limit of split 
# upper             = upper limit of split 
# split.index       = empty list that will contain the splits
# maxn              = smallest gap between splits
# lambda.list       = list of lambdas
# n.rep             = number of permutations or bootstrap replicates
# block             = block size for stationary or block bootstrap
# refit.flag        = TRUE if refit glasso to take out bias
# method_type       = 1 if AIC reduction is of interest for finding the splits
#                   = 2 if BIC reduction is of interest for finding the splits
#                   = 3 if RISK reduction is of interest for finding the splits

low.s       = lower+maxn
upp.s       = upper-maxn

if (low.s <= upp.s) {

	T              = dim(y)[1]
	ind.interval   = which(x<=upper & x>lower)
    	S.interval     = var(y[ind.interval,])
	current        = GraphLassoPath.ALL(S.interval, y[ind.interval,], method_type, refit.flag, lambda.list)
	current.result = current$result
	lambda.full    = current$lambda


	indices   = low.s:upp.s
	
	k           = 0
	result.diff = rep(0,length(indices))
	lambda_mat  = matrix(0,ncol=2,nrow = length(indices))

	for(i in indices){
		
		print(i)
		k = k+1
		left.ind      = which(x<=i & x>lower)
		right.ind     = which(x<=upper & x>i)
		
		temp_l      = GraphLassoPath.ALL(var(y[left.ind,]), y[left.ind,], method_type, refit.flag, lambda.list)
		leftresult  = temp_l$result
		temp_r      = GraphLassoPath.ALL(var(y[right.ind,]), y[right.ind,], method_type, refit.flag, lambda.list)
		rightresult = temp_r$result
	
		result.diff[k]  = current.result-leftresult-rightresult
		lambda_mat[k,1] = temp_l$lambda
		lambda_mat[k,2] = temp_r$lambda
		
	} # indices loop

	if(max(result.diff)<= 0){split.index = split.index}
	
	else {

		new.split = indices[which.max(result.diff)]

			print(paste("Split at =", new.split))
			print(paste("result.diff =",max(result.diff)))
			print(paste("lambda_1 =",lambda_mat[(new.split-low.s+1),1]))
			print(paste("lambda_2 =",lambda_mat[(new.split-low.s+1),2]))
			cc                 = c(new.split,max(result.diff),lambda_mat[(new.split-low.s+1),1],lambda_mat[(new.split-low.s+1),2])
			split.index$actual = rbind(split.index$actual,cc)

			split.index  = Recall(x,y,lower,new.split,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type)
			split.index  = Recall(x,y,new.split,upper,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type)
			
			} # else loop
			
		} # if low.s<upp.s loop

else{split.index=split.index}	

return(split.index)
		
}


#===========================================================================
# This function permutes the whole data set or data within subjects
# Permute whole rows here, that is, permute the time course of all ROIs
# together
# The first element of the permutation will be a 1*n matrix where n is the 
# number of ROIs

permute_split = function(data){

# data   = data to be read in

n1   = dim(data)[1]
k1   = dim(data)[2] # Number of regions

b.y = c()
index = sample(seq(1,n1,1))
for (i in 1:T){
	b.y    = rbind(b.y,data[index[i],])
	}
return(b.y)
}

#===========================================================================
# This function bootstraps the data set 
# The first element of the bootstrap distribution will be a 1*n matrix where 
# n is the number of ROIs

boot_split = function(data){

# data   = data to be read in

n1   = dim(data)[1]
k1   = dim(data)[2] # Number of regions

b.y = c()
index = sample(seq(1,n1,1),replace=T)
for (i in 1:n1){
	b.y    = rbind(b.y,data[index[i],])
	}
return(b.y)

}

#============================================================================
# This function block bootstraps the whole data set 
# The first block element of the resampling will be a block*n matrix where n 
# is the number of ROIs

b_boot_split = function(data,block){

# data   = data to be read in
# block  = length of block, given in % terms

n1   = dim(data)[1]
k1   = dim(data)[2] # Number of regions

block = round(n1*block)
b.y   = NULL
ytemp = rbind(data,data)  
 
for (i in 1:ceil(n1/block)){
	index  = sample(seq(1,n1,1))[1]
	b.y    = rbind(b.y,ytemp[index:(index+block-1),])
}
b.y   = b.y[1:n1,]
return(b.y)
}


#============================================================================
# This function stationary bootstraps the whole data set 
# The first block element of the resampling will be a block*n matrix where n 
# is the number of ROIs

stat_boot_split = function(data,block){

# data   = data to be read in
# block  = length of block, given in % terms

n1   = dim(data)[1]
k1   = dim(data)[2] # Number of subjects

block      = round(n1*block)
ytemp      = data
b.x        = NULL
b.y        = NULL
r.num      = runif(n1, min=0, max=1)
index      = sample(seq(1,n1,1))[1]
temp_index = c(seq(1,n1,1),seq(1,n1,1))
b.x        = c(b.x,index)

for (i in 1:(n1-1)){
	if (r.num[i] < 1 - (1/block)){
			
		b.x   = c(b.x,temp_index[index + 1])
		index = index + 1
			
		}
	else
	if (r.num[i] >= 1/block){
		index  = sample(seq(1,n1,1))[1]
		b.x    = c(b.x,index)
            }
      }
b.x   = b.x[1:n1]
b.y   = rbind(b.y,ytemp[b.x,])
return(b.y)
}

#===========================================================================
# This function computes the graphical lasso solution for several lambda 
# values

GraphLassoPath.ALL = function(S, y, method_type, refit.flag=FALSE, lambda.list=NULL)
{
# S           = variance of y
# y           = y vector, training set
# method_type = 1 if AIC reduction is of interest for finding the splits
#             = 2 if BIC reduction is of interest for finding the splits
#             = 3 if RISK reduction is of interest for finding the splits
# refit.flag  = TRUE to refit glasso
# lambda.list = list of regularization parameters

	out      = glassopath(S, rholist=lambda.list)
	# Adding full covariance for comparison 
      #out.zero = glassopath(S, rholist=0)	
      mu       = colMeans(y)
      n        = nrow(y)
      p        = ncol(y)
      S.out    = y-repmat(mu, n, 1)
      S.out    = as.matrix(t(S.out))%*%as.matrix(S.out)

      Omega        = list()
      lambda.len   = dim(out$wi)[3]
	# Bind the two lists
      #out$rholist  = c(out.zero$rholist,out$rholist)
     	#out$wi[,,1]    = out.zero$wi[,,1]
     	result         = rep(0,lambda.len)

     	for (i in  c(1:lambda.len)){

      	lambda = out$rholist[i]
           	if (refit.flag &&  any(out$wi[,,i] == 0)){
	      	Omega[[i]] = glasso(S, rho=lambda.list[1]/1000, zero = which(out$wi[,,i] == 0, arr.ind = TRUE), thr = 1.0e-2, maxit=1e2, approx=FALSE, penalize.diagonal=TRUE, start="cold",w.init=NULL, wi.init=NULL, trace=TRUE)$wi
           
		} else{Omega[[i]] =  out$wi[,,i]}
		
		if (method_type == 1){
			result[i]    = sum(diag(Omega[[i]]%*%S.out))-n*log(det(Omega[[i]]))+ 2*sum(Omega[[i]][upper.tri(Omega[[i]],diag=TRUE)]!=0)
		
		} else if (method_type == 2){
			result[i]    = sum(diag(Omega[[i]]%*%S.out))-n*log(det(Omega[[i]]))+ log(n)*sum(Omega[[i]][upper.tri(Omega[[i]],diag=TRUE)]!=0)
		
		} else if (method_type == 3){
			result[i]    = sum(diag(Omega[[i]]%*%S.out))-n*log(det(Omega[[i]])) 
		}

	      #n.edge[i] = (sum(sum(Omega[[i]]!=0))-p)/2
		if (result[i]=="NaN"){result[i]=2^64}
     }
 	
	ind = which(result==min(result))[1]

     #ind = which(BIC == min(BIC))
     #print( (sum(sum(Omega[[ind]]!=0))-p)/2)
     return(list(Omega=Omega[[ind]], result=result[ind], lambda=out$rholist[ind],ind=ind))     
}


#===========================================================================
# Recalculates the BIC reduction between each pair of splits
# Calculates the permutation distribution
# Computes the undirected graph between each pair of significant change points
# Finds the bootstrap distribution for the edges in the undirected graphs

refit.split.perm = function(splits,perm_type,method_type,datasubj,x,T,lambda.list,block,n.rep=1000){
					
# splits              = splitting times
# perm_type           = 0 no bootstrap replications
#	                = 1 for simple bootstrap
#                     = 2 for block bootstrap
#                     = 3 for stationary bootstrap
# method_type         = 1 if AIC reduction is of interest for finding the splits
#                     = 2 if BIC reduction is of interest for finding the splits
#                     = 3 if RISK reduction is of interest for finding the splits
# datasubj            = data for one subject
# T                   = experimental length
# lambda.list         = list of lambdas for glasso
# name_results        = original splits
# name_boot_results   = split from bootstrapping or permutations
# block               = size of block for stationary and block bootstrap
# n.rep               = number of replications for bootstrapping

results = list()

split.times = splits$actual[,1]

order.split = c()
ccc = order(split.times)
for (i in 1:length(split.times)){	
	order.split = c(order.split,split.times[ccc[i]])
	}
# Order them
order.split  = c(0,order.split)
order.split  = c(order.split,T)
	
for (jrl in 1:(length(order.split)-2)){

	print(jrl)
	set.seed(jrl*11)
	lower1         = order.split[jrl]
	upper1         = order.split[jrl+2]
	ind.interval   = which(x<=upper1 & x>lower1)
	S.interval     = var(datasubj[ind.interval,])
	current        = GraphLassoPath.ALL(S.interval, datasubj[ind.interval,], method_type, refit.flag, lambda.list)
	current.result = current$result

	left.ind      = which(x<=order.split[jrl+1] & x>lower1)
	right.ind     = which(x<=upper1 & x>order.split[jrl+1])
		
	temp_l      = GraphLassoPath.ALL(var(datasubj[left.ind,]), datasubj[left.ind,], method_type, refit.flag, lambda.list)
	leftresult  = temp_l$result
	temp_r      = GraphLassoPath.ALL(var(datasubj[right.ind,]), datasubj[right.ind,], method_type, refit.flag, lambda.list)
	rightresult = temp_r$result
	
	result.diff  = current.result-leftresult-rightresult
	lambda_mat_l = temp_l$lambda
	lambda_mat_r = temp_r$lambda
		
	dd             = c(order.split[jrl+1],result.diff,temp_l$lambda,temp_r$lambda) 
	results$actual = rbind(results$actual,dd)

	#-----------------------------------------------------------
      # Now do the permutation/bootstrapping for this split
			
	if (perm_type > 0){
		for (j in 1:n.rep){
					
			set.seed(11*j)

			if (perm_type == 1) {
				new.data = permute_split(datasubj[ind.interval,])

			} else if (perm_type == 2) {
				new.data = b_boot_split(datasubj[ind.interval,],block)

			} else new.data = stat_boot_split(datasubj[ind.interval,],block)
					
					x.new = matrix((lower1+1):upper1,ncol=1)
					# Link time points in new.data matrix
					row.names(new.data)     <- x.new

    					S.interval            = var(new.data)
					
					current.junk              = GraphLassoPath.ALL(S.interval, new.data, method_type, refit.flag, lambda.list)
					current.result.perm.boot  = current.junk$result
	
					left.ind      = which(x.new<=order.split[jrl+1] & x.new>lower1)
					l.ind         = length(left.ind)
					right.ind     = which(x.new<=upper1 & x.new>order.split[jrl+1])
					r.ind         = length(right.ind)

					temp.l       = GraphLassoPath.ALL(var(new.data[1:l.ind,]), new.data[1:l.ind,], method_type, refit.flag, lambda.list)
					left.result  = temp.l$result
					temp.r       = GraphLassoPath.ALL(var(new.data[((l.ind+1):(r.ind+l.ind)),]), new.data[((l.ind+1):(r.ind+l.ind)),], method_type, refit.flag, lambda.list)
					right.result = temp.r$result

					new.data.result.diff  = current.result.perm.boot-left.result-right.result
					ddd                   = c(order.split[jrl+1],new.data.result.diff,temp.l$lambda,temp.r$lambda) 
					results$boot.perm = rbind(results$boot.perm,ddd)

					} # j loop
			} # if type loop
	}

return(results)

}

#===========================================================================
# Finds significant splits using the permutation distribution
# Computes the undirected graph between each pair of significant change points
# Finds the bootstrap distribution for the edges in the undirected graphs

sign.splits.boot.edges = function(boot_type,method_type,data,T,lambda.list,split.index,quantile_1,quantile_2,block,n.rep=1000){

# boot_type           = 0 no bootstrap replications
#	                = 1 for simple bootstrap
#                     = 2 for block bootstrap
#                     = 3 for stationary bootstrap
# method_type         = 1 if AIC reduction is of interest for finding the splits
#                     = 2 if BIC reduction is of interest for finding the splits
#                     = 3 if RISK reduction is of interest for finding the splits
# data                = data
# T                   = number of experimental time points
# lambda.list         = list of lambdas for glasso
# split.index         = results from refitting
# quantile_1          = lower quantile for confidence bounds of bootstrapping or 
#                       permutations. 
# quantile_2          = upper quantile for confidence bounds of bootstrapping or 
#                       permutations. 
# block               = size of block for stationary and block bootstrap
# n.rep               = number of replications for bootstrapping

results = list()
output  = split.index$actual
output  = as.matrix(output)

junk = rep(0,T)
n    = dim(output)[1]
for (i in 1:n){
	junk[c(output[i,1])] = output[i,2]
	}

# Stationary bootstrap or permutations confidence bounds
stat = split.index$boot.perm
stat = as.matrix(stat)

split.times = output[,1]

garbage = matrix(0,ncol=2,nrow=T)

for (i in c(split.times)){
	index = which(stat[,1] == i)
	if (quantile(stat[index,2],prob = quantile_1)>0)
		{garbage[i,1] = quantile(stat[index,2],prob = quantile_1)}
	if (quantile(stat[index,2],prob = quantile_2)>0)
		{garbage[i,2] = quantile(stat[index,2],prob = quantile_2)}
	}

junk.output = c()
for (i in 1:dim(output)[1]){
	if (output[i,2] > garbage[split.times[i],2]){
		junk.output = rbind(junk.output,output[i,])
		}					
	}

	if (length(junk.output)>0){ # for examples where there are no significant splits

	new.output = c()
	ccc = order(junk.output[,1])
	for (i in 1:dim(junk.output)[1]){	
		new.output = rbind(new.output,junk.output[ccc[i],])
		}

	results$sign.splits = new.output
	
	# To check for significant correlations
	new.output  = rbind(c(0,0,0,0),new.output)
	new.output  = rbind(new.output,c(T,0,0,0))
	
	for (i in 1:(dim(new.output)[1]-1)){

		print(i)
		lower        = new.output[i,1]
		upper        = new.output[i+1,1]
		ind.interval = which(x<=upper & x>lower)
	
		ccc = final.edges.graph.boot(method_type,data[ind.interval,],ind.interval,lambda.list)
	
		results$partial.corr.matrix[[i]] = ccc

		# Loop to calculate significant partial correlations in undirected
		# graphs
		boot.entries = c()

		for (j in 1:n.rep){

			set.seed(11*j)
			if (boot_type == 1) {
				new.data = boot_split(data[ind.interval,])
	
				} else if (boot_type == 2) {
					new.data = b_boot_split(data[ind.interval,],block)
	
				} else if (boot_type == 3) {
					new.data = stat_boot_split(data[ind.interval,],block)
						
				}

				ivor    = final.edges.graph.boot(method_type,new.data,ind.interval,lambda.list)
				n_omega = dim(ivor)[1]

				for (u in 1:(n_omega-1)){
					for (o in (u+1):n_omega){
							boot.entries = c(boot.entries,ivor[u,o])
						} # o loop
					} # u loop	
				

			}# j loop
			
			results$corr.quants[[i]] = boot.entries
			
	} # i loop
	
	} else {    # if (length(junk.output)>0) loop
	

	results = sign.splits.boot.edges.iid(boot_type,method_type,data,T,lambda.list,block)

	}
return(results)
}

#===========================================================================
# Finds significant splits using the permutation distribution
# Computes the undirected graph between each pair of significant cahnge points
# Finds the bootstrap distribution for the edges in the undirected graphs
# THIS IS FOR IID CASE ONLY

sign.splits.boot.edges.iid = function(boot_type,method_type,data,T,lambda.list,block,n.rep=1000){

# boot_type           = 0 no bootstrap replications
#	                = 1 for simple bootstrap
#                     = 2 for block bootstrap
#                     = 3 for stationary bootstrap
# method_type         = 1 if AIC reduction is of interest for finding the splits
#                     = 2 if BIC reduction is of interest for finding the splits
#                     = 3 if RISK reduction is of interest for finding the splits
# data                = data (subjects stacked)
# T                   = experimental length
# lambda.list         = list of lambdas for glasso
# block               = size of block for stationary and block bootstrap
# n.rep               = number of replications for bootstrapping

results = list()

lower = 0
upper = T	

ind.interval = which(x<=upper & x>lower)

ccc   = final.edges.graph.boot(method_type,data[ind.interval,],ind.interval,lambda.list)

results$partial.corr.matrix[[1]] = ccc

# Loop to calculate bootstrap distribution of partial correlations in the
# undirected graph

boot.entries = c()

for (j in 1:n.rep){

	set.seed(11*j)
	if (boot_type == 1) {
		new.data = boot_split(data[ind.interval,])
	
	} else if (boot_type == 2) {
		new.data = b_boot_split(data[ind.interval,],block)
	
	} else if (boot_type == 3) {
		new.data = stat_boot_split(data[ind.interval,],block)
						
	}

	ivor    = final.edges.graph.boot(method_type,new.data,ind.interval,lambda.list)
	n_omega = dim(ivor)[1]

	for (u in 1:(n_omega-1)){
		for (o in (u+1):n_omega){
			boot.entries = c(boot.entries,ivor[u,o])
			} # o loop
		} # u loop	
				

	}# j loop
			
	results$corr.quants[[1]] = boot.entries
	
return(results)
}

#============================================================================
# This function takes the partitons and recalculates the glasso along the
# full path of lambda values
# This is carried out in order to create the undirected graphs with the
# correct data
# It outputs the partial correlation matrix as well as a bootstrap
# distribution for each partition

final.edges.graph.boot = function(method_type,data.ind,ind.interval,lambda.list,refit.flag=TRUE){

# method_type  = 1 if AIC reduction is of interest for finding the splits
#              = 2 if BIC reduction is of interest for finding the splits
#              = 3 if RISK reduction is of interest for finding the splits
# data         = data 
# ind.interval = x values or interval to calculate graph
# lambda.list  = list of lambda values to reduce over

T.ind  = dim(data.ind)[1]
k      = dim(data.ind)[2]

x_data = ind.interval
x_data = as.matrix(x_data)

new.data.ind = data.ind
new.data.ind = as.matrix(new.data.ind)
row.names(new.data.ind) = x_data

S.interval       = var(new.data.ind)
current          = GraphLassoPath.ALL(S.interval, new.data.ind, method_type, refit.flag, lambda.list)
lambda.full      = current$lambda
	
Omega   = glasso(S.interval,rho=lambda.full,zero=NULL)$wi

# Refit if there is shrinkage
junk = which(Omega==0)
if (length(junk)>0){
	Omega = glasso(S.interval,rho=lambda.list[1]/1000,zero=which(Omega==0,arr.ind=TRUE))$wi
}

cor.matrix = matrix(0,ncol=k,nrow=k)
for (i in 1:dim(Omega)[1]){
	for(j in 1:dim(Omega)[2]){
		cor.matrix[i,j] = Omega[i,j]/sqrt(Omega[i,i]*Omega[j,j])
	}
}		 

cor.matrix[lower.tri(cor.matrix, diag=TRUE)] = 0

return(cor.matrix)

}

#===========================================================================
# This function converts the split.index output into readable form
# This is for output from cluster

convert.split = function(split.index){

# split.index = output from main split.index function

split.index = matrix(split.index,ncol=4)

n    = dim(split.index)[1]
k    = dim(split.index)[2]

junk = matrix(0,ncol=k,nrow=n)

if  (dim(split.index)[1] > 3){
	if (dim(junk)[1] %% 4 == 0){
	for (i in 1:(n/4)){
		junk[(1+(i-1)*4):(4+(i-1)*4),1] = split.index[i,]
		junk[(1+(i-1)*4):(4+(i-1)*4),2] = split.index[(i+(n/4)),]
		junk[(1+(i-1)*4):(4+(i-1)*4),3] = split.index[(i+(n/2)),]
		junk[(1+(i-1)*4):(4+(i-1)*4),4] = split.index[(i+(3*n/4)),]
		}
	
	} else if (dim(junk)[1] %% 4 == 1){
		 if (dim(junk)[1] == 5){
	
			for (i in 1:(floor(n/4))){
				junk[(1+(i-1)*4):(4+(i-1)*4),1] = split.index[i,]
				junk[(2+(i-1)*4):(1+(i*4)),4]   = split.index[(i+ceil(3*n/4)),]

			}
			#for (i in 1:(floor(n/4)-1)){
			#	junk[(4+(i-1)*4):(3+(i*4)),2]   = split.index[(i+ceil(n/4)),]
			#	junk[(3+(i-1)*4):(2+(i*4)),3]   = split.index[(i+ceil(n/2)),]
			#	}
			junk[n:n,1]     = split.index[floor(n/4)+1,1]
			junk[1:3,2]     = split.index[floor(n/4)+1,2:4]
			junk[(n-1):n,2] = split.index[floor(n/2)+1,1:2]
			junk[1:2,3]     = split.index[floor(n/2)+1,3:4]
			junk[(n-2):n,3] = split.index[floor(3*n/4)+1,1:3]
			junk[1,4]       = split.index[floor(3*n/4)+1,4:4]

		} else if (dim(junk)[1] > 5){
						
			for (i in 1:(floor(n/4))){
				junk[(1+(i-1)*4):(4+(i-1)*4),1] = split.index[i,]
				junk[(2+(i-1)*4):(1+(i*4)),4]   = split.index[(i+ceil(3*n/4)),]

			}
			for (i in 1:(floor(n/4)-1)){
				junk[(4+(i-1)*4):(3+(i*4)),2]   = split.index[(i+ceil(n/4)),]
				junk[(3+(i-1)*4):(2+(i*4)),3]   = split.index[(i+ceil(n/2)),]
				}
			junk[n:n,1]     = split.index[floor(n/4)+1,1]
			junk[1:3,2]     = split.index[floor(n/4)+1,2:4]
			junk[(n-1):n,2] = split.index[floor(n/2)+1,1:2]
			junk[1:2,3]     = split.index[floor(n/2)+1,3:4]
			junk[(n-2):n,3] = split.index[floor(3*n/4)+1,1:3]
			junk[1,4]       = split.index[floor(3*n/4)+1,4:4]

	}
	} else if (dim(junk)[1] %% 4 == 2){
		for (i in 1:(n/4)){
			junk[(1+(i-1)*4):(4+(i-1)*4),1] = split.index[i,]
			junk[(3+(i-1)*4):(2+(i*4)),2] = split.index[(i+ceil(n/4)),]
			junk[(1+(i-1)*4):(4+(i-1)*4),3] = split.index[(i+ceil(n/2)),]
			junk[(3+(i-1)*4):(2+(i*4)),4] = split.index[(i+ceil(3*n/4)),]
			}
		junk[(n-1):n,1] = split.index[floor(n/4)+1,1:2]
		junk[1:2,2]     = split.index[floor(n/4)+1,3:4]
		junk[(n-1):n,3] = split.index[floor(3*n/4)+1,1:2]
		junk[1:2,4]     = split.index[floor(3*n/4)+1,3:4]
		
	}	else if (dim(junk)[1] %% 4 == 3){
		for (i in 1:(floor(n/4))){
			junk[(1+(i-1)*4):(4+(i-1)*4),1] = split.index[i,]
			junk[(2+(i-1)*4):(1+(i*4)),2] = split.index[(i+ceil(n/4)),]
			junk[(3+(i-1)*4):(2+(i*4)),3] = split.index[(i+ceil(n/2)),]
			junk[(4+(i-1)*4):(3+(i*4)),4] = split.index[(i+ceil(3*n/4)),]
			}
		junk[(n-2):n,1] = split.index[floor(n/4)+1,1:3]
		junk[1,2]       = split.index[floor(n/4)+1,4:4]
		junk[(n-1):n,2] = split.index[floor(n/2)+1,1:2]
		junk[1:2,3]     = split.index[floor(n/2)+1,3:4]
		junk[n:n,3]     = split.index[floor(3*n/4)+1,1:1]
		junk[1:3,4]     = split.index[floor(3*n/4)+1,2:4]
	
			}
} else {
	if (dim(junk)[1] == 1) {junk = split.index

} else if (dim(junk)[1] == 2){ 
		junk[1:2,1] = split.index[1,1:2]
	     	junk[1:2,2] = split.index[1,3:4]
		junk[1:2,3] = split.index[2,1:2]
		junk[1:2,4] = split.index[2,3:4]
	
} else if (dim(junk)[1] == 3){ 
		junk[1:3,1] = split.index[1,1:3]
		junk[1,2]   = split.index[1,4]
		junk[2:3,2] = split.index[2,1:2]
	      junk[1:2,3] = split.index[2,3:4]
		junk[3,3]   = split.index[3,1]
		junk[1:3,4] = split.index[3,2:4]

} else if (dim(junk)[1] == 4){ 
		junk[1:4,1] = split.index[1,1:4]
		junk[1:4,2] = split.index[2,1:4]
		junk[1:4,3] = split.index[3,1:4]
		junk[1:4,4] = split.index[4,1:4]
		}

}
	return(junk)
}

#============================================================================
#============================================================================
