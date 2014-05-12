#========================================================================
#========================================================================
# Set working directory
setwd("/home/data/Projects/Zhen/DCR")

# Load Functions
# First, you need to load the following packages: 'glasso' and 'matlab'
library(glasso)
library(matlab)
source('SSDCR_FUNCTIONS.R')

# Read in data as a .txt file
# Subjects should be stacked on top of one another
# For example, if there are 215 time points, 23 subjects and 10 ROIs
# the .txt file should contain 4945 rows and 10 columns
data = read.table('ROIs_4_HR.txt')

# Convert the data into a list of subjects
data   = data
n.subj = 23
T      = 215
data   = datasetup(data,n.subj,T)

# Inputs for splitting function
# Time
x = matrix(1:T,ncol=1)
# Starting point for the algorithm
lower             = 0
# Total number of time points in the experiment
upper             = 215
# A list for output
split.index       = list()
# The minimum distance between change points
maxn              = 40
# The list of regularization parameters
lambda.list       = seq(from=0.1,to=1, length.out=10)
# Number of replications in the bootstrap procedures
n.rep             = 10
# The block length for the stationary bootstrap in percentage terms
block             = 0.05 
# Refit the glasso estimation procedure to reduce bias
refit.flag        = TRUE
# Use BIC to find change points
method_type       = 2 # BIC
# Use the stationary bootstrap procedure for inference on the change points 
perm_type         = 3 
# Lower quartile for BIC reduction   
quantile_1        = 0.025
# Upper quartile for BIC reduction 
quantile_2        = 0.975

dcrfunction(x,data,n.subj,T,lower,upper,split.index,maxn,lambda.list,n.rep,block,refit.flag,method_type,perm_type,quantile_1,quantile_2)

# For each subject:
# 1. The candidate change points are saved as 'ActualSplitsSubject_'
# 2. The bootstrap distribution for inference on the change points are saved as
# 'BootSplitsSubject_'
# 3. The significant change points are saved as 'SignificantSplitSubject_'
# 4. The precision matrices (undirected graphs) are saved as 
# 'PrecisionMatricesSubject_'
# 5. The bootstrap distribution for the edges in the undirected graphs are
# saved as 'EdgePermDistSubject_'

#===========================================================================
# Plot significant change points against subjects

# You need to load the package 'sna'
# Load functions
source('NETWORK_FUNCTIONS.R')

# To plot the significant change points against the subjects (no plot, put 0)
changepointsplot = 1
# To plot the undirected graphs between each pair of significant change points
# for each subject (no plot, put 0)
undirectedgraphs = 1
# A list of subjects with signigicant change points (this is an output from
# the dcrfunction above)
subsignsplits    = read.table("SubjectswithSignificantSplits")$x
# The number of time points in the experiment
T                = 215
# The number of subjects
n.subj           = 23
# Cutoff for inference on the edges in the graphs
percentile       = 0.75
# Labels for the vertices in the undirected graphs 
lab_row          = c(1,2,3,4,"HR")

dcrplots(changepointsplot,undirectedgraphs,subsignsplits,T,n.subj,percentile,lab_row,n.rep)
graphics.off()

# The output is found in your directory. It includes:
# 1. A plot of significant change points against time for each subject (saved
# as "ChangePointPlot.png").
# 2. For each subject, a undirected graph between each pair of sigificant
# change points. They are saved under "UndirectedGraphforSubject__Partition_.png"

#============================================================================
#============================================================================
