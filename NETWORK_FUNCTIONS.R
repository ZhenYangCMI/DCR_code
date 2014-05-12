#===========================================================================
# Library
library(sna)

#===========================================================================
# This function plots the change points for each subject and the undirected
# graphs

dcrplots = function(changepointsplot,undirectedgraphs,subsignsplits,T,n.subj,percentile,lab_row,n.rep){

# changepointsplot = 1 if plotting the significant change points for each subject
# undirectedgraphs = 1 if plotting the undirected graphs between each pair of
#                    change points
# subsignsplits    = the subject's with significant change points
# T                = number of time points in the experiment
# n.subj           = number of subjects
# percentile       = cutoff for the inference on the edges in the graphs
# lab_row          = labels for the vertices in the undirected graphs 
# n.rep            = Number of replications in the bootstrap procedures

if (changepointsplot>0){

	png("ChangePointPlot")
	subsplits   = c()
	yaxis       = c()
	subjectlist = subsignsplits 
	for (i in c(subjectlist)){
	
		filename   = sprintf("SignificantSplitsSubject_%s",i)
		splitstemp = read.table(filename)
		subsplits  = c(subsplits,splitstemp[,1])
		nn         = length(splitstemp[,1])
		yaxis      = c(yaxis,rep(i,nn))
	}
	plot(subsplits,yaxis,type="p",pch=4,xlab = "splits",ylab = "subject",xlim=c(0,215),axes=F,
		frame.plot=T,lwd=4)
	axis(1, at = seq(1,T,10))
	axis(2, at = 1:n.subj)
	for (i in c(subjectlist)){
	abline(h=i,lwd=1,lty=2)
	}

}
# Undirected graphs
if (undirectedgraphs>0){

	for (jkj in 1:n.subj){
		# load upper triangular partial correlation matrices
		filename   = sprintf("PrecisionMatricesSubject_%s",jkj)
		cor.matrix = read.table(filename)
	
		nm = dim(cor.matrix)[2]
		# load perm edges
		filename1    = sprintf("EdgePermDistSubject_%s",jkj)
		perm.edge    = read.table(filename1)

		color      = c(rep("blue",nm))
		v_scaler   = 35000

		# Windows
		plot.network.perm.edge.boot(cor.matrix,lab_row,color,perm.edge,percentile,jkj,n.rep)

		rm(perm.edge)
		}
	}
}
#===========================================================================
# This function plots the undirected graph with inference on the edges
# provided by the bootstrap distribution

plot.network.perm.edge.boot = function(cor.matrix,lab_row,color,perm.edge,percentile,jkj,n.rep){

# cor.matrix = stacked partial correlations matrix of sign. change points
# n.ser      = number of series in model
# lab_row    = labeling of vertices
# color      = color of vertices
# perm.edge  = the permutation distribution for the edges
# percentile = quantile to cut the permutation distribution at
# jkj        = subject in the title of undirected graphs
# n.rep      = Number of replications in the bootstrap procedures

n.ser       = dim(cor.matrix)[2]
k.n         = dim(cor.matrix)[1]

num.cor.mat = k.n/n.ser

n.perm.edge = dim(perm.edge)[1]
n.edge      = n.ser*(n.ser-1)/2

for (k in 1:num.cor.mat){
	
	#windows()
	cor.junk         = cor.matrix[(1+(k-1)*n.ser):(k*n.ser),]

	b.y = c()
	for (j in 1:n.edge){
		junk = seq(from =j, to=n.perm.edge, by = n.edge )
		b.y  = cbind(b.y,perm.edge[junk,k])
	}

	junk   = c()
	for (i in 1:n.edge){
		index = which(b.y[,i]!=0)
		junk = c(junk,length(index)/n.rep)
	}

	# significant correlations
	sign.matrix = matrix(0,ncol=n.ser,nrow=n.ser)
	iter        = 0
	for (u in 1:(n.ser-1)){
		for (o in (u+1):n.ser){
			iter = iter+1
			if( junk[iter] > percentile ){
				sign.matrix[u,o] = cor.junk[u,o]
				} # if loop
			} # o loop
		} # u loop	
	
	# Red for negative correlation/black for positive correlation
	garbage = matrix(1,ncol=n.ser,nrow=n.ser)
		for (i in 1:n.ser){
			for (j in 1:n.ser){
				if (sign.matrix[i,j]>0){
					garbage[i,j] = 2
				}
			}
	}	 

	garbage[lower.tri(garbage, diag=TRUE)] = 0
	row.names(garbage) = lab_row

	cor.matrixUpperTri = abs(sign.matrix)
	cor.matrixUpperTri[lower.tri(cor.matrixUpperTri, diag=TRUE)] = 0
	row.names(sign.matrix) = lab_row
	
	# Added if only one network in graph	
	cor.matrixUpperTri[1,1] = 1
	garbage[1,1] = 0
	
	newfilename = sprintf("UndirectedGraphforSubject_%s_Partition_%s",jkj,k)
	png(newfilename)

	gplot(cor.matrixUpperTri, g = 1, gmode = "graph",label = lab_row, coord = NULL, jitter = FALSE, thresh = 0,
	diag = F, mode = "circle", displayisolates = TRUE, interactive = FALSE, xlab = NULL, 
	ylab = NULL ,xlim =NULL, ylim = NULL, pad =0.5, label.pad = 0.27, displaylabels =T, 
	boxed.labels = F, label.pos = 6, label.bg = "white", vertex.sides = 10,
	arrowhead.cex = 1, label.cex = 2, loop.cex = 1, vertex.cex=2, label.col = 1, vertex.col = color, 
	label.border = 1,vertex.rot=c(45), vertex.border = 1, edge.lty = 1, label.lty = NULL, vertex.lty = 1, edge.lwd =cor.matrixUpperTri*50,
	label.lwd = par("lwd"), edge.col = garbage,edge.len = 0.5, edge.curve = 0.01, edge.steps = 50, loop.steps = 20,  object.scale = 0.01,
	uselen = F, usecurve =F, suppress.axes = TRUE, vertices.last = F, new = TRUE, layout.par = NULL)
	#filename  = sprintf("Subject_%s_Partition_%s",jkj,k)
	#title(main=filename)
	dev.off()
	}
}

#===========================================================================
#=============================================================================
# Parametric Function for testing whether a series of precision matrices
# are equivalent without converting to covariance matrices and taking out
# the matrices that are not positive definite

parametric.test.no.convert = function(results.subj,nsample,which.index){

# results     = a list containing all undirected graphs for all subjects
# nsample     = data sample size, vector
# which.index = which undirected graphs are you testing

# s is the number of covariance matrices
s = length(results.subj)

list.precision = list()
for (i in 1:s){
	list.precision[[i]] = results.subj[[i]][[which.index]]
}

# A matrix is positive definite if it is symmetric and all its eigenvalues
# are positive
remove.mats = c()
for (i in 1:length(list.precision)){
	n.eigen = length(eigen(list.precision[[i]], only.values=1)$values)
	junk = (eigen(list.precision[[i]], only.values=1)$values)>0
	if (n.eigen != sum(junk)){
		print(i)
		remove.mats = c(remove.mats,i)
		}	
	}

new.list.precision = list()
if (length(remove.mats)>0){
	indices = (1:length(nsample))[-remove.mats]
	iter    = 0
	for (i in indices){
		iter = iter + 1
		new.list.precision[[iter]] = list.precision[[i]]
	}
	nsample = nsample[-remove.mats] 
} else {new.list.precision = list.precision}

# p is the number of rows or columns of covariance matrix
p = dim(new.list.precision[[1]])[1]

nnn = sum(nsample)

k = -0.5*nnn*p*log(2*pi) - 0.5*nnn*p

Omega_0 = matrix(0,ncol=p,nrow=p)
s.new   = s - length(remove.mats)
for (i in 1:p){
	for (j in 1:p){
		junk = c()
		garbage = c()
		for (k in 1:s.new){
			junk = c(junk,new.list.precision[[k]][i,j])
			garbage = c(garbage,(nsample[k]*junk[k])/nnn)
			}
		Omega_0[i,j] = sum(garbage)
	}
}

l_0 = k - 0.5*nnn*log(det(Omega_0))

naom = c()
for (i in 1:s.new){
	junk = 0.5*nsample[i]*log(det(new.list.precision[[i]]))
	naom = c(naom,junk)
}
l_3 = k - sum(naom)

T   = 2*(l_3 - l_0)
dof = 0.5*(s.new-1)*p*(p+1)

p.val = 1-pchisq(T, dof)

result = list()
result$l_0 = l_0
result$l_3 = l_3
result$T = T
result$p.val = p.val
return(result)
}

#===========================================================================
#===========================================================================