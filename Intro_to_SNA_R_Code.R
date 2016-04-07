## Introduction to SNA with R
## Breanne Chryst and 
## CSSSI StatLab
## Feb 19, 2016

################################################################################
## R Basics ##
################################################################################
## R can be used as a calculator, it works as expected:
2+3
exp(2)
5^(2)

## Assigning a variable
x <- 5 # 5 has now been assigned to the variable x
x
x^2

## Creating a vector:
y <- c(3,7,5,1,2,3,2,5,5) # "c()" concatenates, creating a vector

## Extracting values of a vector:
y[2]

3:5 # the whole numbers from 3 to 5

y[3:5]

## "matrix()" creates a matrix from the values entered:
z <- matrix(y, nrow=3) # This is filled by column
z 

z <- matrix(y, nrow=3, byrow=T) 
# By changing the "byrow" option, we can fill the matrix by row
z 

## Extracting values from matrices:
z[2,] # Row
z[,3] # Column
z[2,3] # Value

## Create Dataframes:
dat <- as.data.frame(z)
names(dat) <- c("cat", "giraffe","bowlingball")
dat

## R has base functions:
mean(y)
length(y)
sd(y)
var(y)
prod(y) # Takes the product of each element in the vector
apply(z, 2, mean) # Very useful in avoiding for loops, also has useful cousins sapply and lapply

################################################################################
## Brief Introduction to Statistics with R ##
################################################################################
## Getting help
#help.start()	# Opens html help in web browser (if installed)
#help(help)	# find help on how to use help
#?help		# Same as above
#help.search("help")	# Find all functions that include the word 'help'

## Reading in your data
getwd()		# What directory are we in?
#setwd("~/Desktop")	
## Set working directory to the directory where we put the data	
dat <- read.table("http://www.stat.yale.edu/~blc3/IntroR2015/remote_weight.txt", header=T, sep="", row.names=NULL, as.is = TRUE)	
# Read data including headers, data separated by spaces, no row names
ls()			# List all variables stored in memory
head(dat) # Shows the first 6 rows of the data
head(dat, 10) # the first 10 rows of the data
tail(dat) # last 6 rows of the data

## Extracting data from the data frame
dim(dat)	        # Find out how many rows and columns in the data set
names(dat)	        # List all variable names in the dataset 
str(dat)           # Look at the structure of your data
# dat		        # See the data frame on the screen
dat[1:5,]	        # See the first 5 rows
dat[,"weight"]	    # See only the weight column
dat[,3]	        # Same as above
dat$weight	        # Yet another way
dat[1:5,"weight"]	# See only the first 10 values of the weight col.
#dat[,-1]	        # See all but the first column of data
dat.o <- dat	    # Copy the data frame to a data.frame named data.0
ls()			# Now we have 5 variables: 'x', 'y', 'z', 'data' and 'data.o'

## Getting familiar with the data
summary(dat)		# Generate summary statistics of data
apply(dat, 2, sd)		# Calculate standard deviations of all variables
var(dat)		# Variance on diagonal, covariance off diagonal
pairs(dat)	 	# A general view of data through scatter plots
pairs(dat[,-1]) 	# See scatterplots for all pairs of variables except the first ('id') in the data frame
plot(dat$weight, dat$remote)	# Scatterplot of 'weight' vs. 'remote'

# Changing data type
class(dat$gender)  		# What kind of variable is 'gender'?
dat$gender <- factor(dat$gender)	# Converts 'gender' from type integer to factor 
class(dat$gender)		# Verify that gender is now indeed of type factor
dat$gender			# See all data in column 'gender'; note "Levels: 0 1" at the bottom

# Basic Graphics
hist(dat$remote)			# Histogram of 'remote'
hist(dat$weight)			# Histogram of 'weight'
boxplot(dat$remote, dat$weight)  	# Boxplot of 'remote' and 'weight' 
boxplot(remote ~ gender, data = dat)  	# Boxplot of 'remote' conditioned on 'gender'
boxplot(weight ~ gender, data = dat, main = "Weight by Gender", 
        xlab = "Gender", ylab = "Weight", 
        col = c("forestgreen", "tomato"))	# Boxplot of 'weight' conditioned on 'gender'

################################################################################
## Introduction to SNA in R ##
################################################################################
# Installing the packages to be used in the analysis
#install.packages("igraph")
#install.packages("igraphdata")
# loading package igraph
library(igraph) # functions for igraph: http://igraph.org/r/doc/
# loading package igraphdata
library(igraphdata)

#A Simple Example:
g <- graph(c(1,2, 2,3, 3,4, 4,5, 3,1, 4,7, 2,5, 
              6,1, 3,6, 7,8, 7,9, 
              8,9, 8,10), directed=F)
g           # summary information
plot(g) # network picture

################################################################################
## Reading in data and creating plots ##
################################################################################
# Read in the edgelist I made up:
dat <- read.csv("http://www.stat.yale.edu/~blc3/SNA2016/Partners.csv")
# Make a igraph object from the edgelist
partners <-  graph.data.frame(dat, directed=F)

# Add edge colors determined by the variable "Type"
E(partners)$color <- c("black", "red", "green")[as.numeric(dat$Type)]

# Plotting my social network
plot(partners, vertex.size=30, edge.color = E(partners)$color, 
     vertex.label.color = "white",
     vertex.color="darkred", edge.width=3 )

# Read in the adjacency version of the data
dat1 <- read.csv("http://www.stat.yale.edu/~blc3/SNA2016/AdjMat.csv", header = T)
# Create an igraph object from the matrix
g1 <- graph_from_adjacency_matrix(as.matrix(dat1[,-1]), mode = "undirected")
# Plot the network
plot(g1, vertex.size=25,vertex.label.color = "white",
     vertex.color="darkgreen")

# Read in the Enron email data and plot
data(enron)
#enron <- delete_vertices(enron, c(72, 118))

#plot(enron, vertex.size=5, vertex.label=NA,
#     edge.arrow.size=.5)

# Read in the kite data and plot
data(kite)
plot(kite)

# Read in an existing Media data, first the nodes, then the edges
media.node <- read.csv("http://www.stat.yale.edu/~blc3/SNA2016/Media-NODES.csv")
head(media.node)

media.edge <- read.csv("http://www.stat.yale.edu/~blc3/SNA2016/Media-EDGES.csv")
head(media.edge)

# Make a graph out of the node and edge files
media <- graph.data.frame(media.edge, media.node, directed=T)

# The network object in R
media

# Plot of the network
plot(media, edge.arrow.size=.2, edge.color="goldenrod",
     vertex.color="goldenrod", vertex.frame.color="white",
     vertex.label=V(media)$media) 

# extracting adjacency matrix

# full adjacency matrix
get.adjacency(partners, sparse=F)
get.adjacency(partners, sparse=T)

# only upper triangle ('g' is undirected)
get.adjacency(partners, type="upper", sparse=FALSE)

# extracting edgelist
get.edgelist(partners)
get.edgelist(kite)

# Setting edge attributes:
E(partners)$weight <- runif(25,1,5)

# Setting vertex attributes:
V(partners)$gender <- c("M", "M", "M", "F", "M", "M", "M", "M", "M", "F", "M")
V(partners)$color= ifelse(V(partners)$gender == "M", "tomato", "gold")

# Network object with new attributes
partners

################################################################################
## Introduction to Social Network Visualization in R ##
################################################################################
# see ?igraph.plotting for detailed explanation of all
# the options

## Layouts
# some available layouts
# Default is Fruchterman-Reingold
# circle layout
plot(g, layout=layout.circle)
# Kamada-Kawai
plot(g, layout=layout.kamada.kawai,
     vertex.label=NA, vertex.size=5, edge.arrow.size=0.5)
# Multidimensional Scaling
plot(g, layout=layout.mds)

# Simulate a random network for visualiation:
net.bg <- barabasi.game(80) 

# Setting node and edge attributes for nice graphs
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- sample(c("orange", "blue"),80, replace = T)
V(net.bg)$label <- "" 
V(net.bg)$size <- 10
E(net.bg)$arrow.mode <- 0

# Base plot
plot(net.bg)

# All the possible layout options in R
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE) 

# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|spring|grid.3d|svd|fruchterman.reingold.grid", layouts)]

# Setting the graph window to a 2 by 2
par(mfrow=c(2,2))

# For loop to plot our simulated network under several layouts
for (layout in layouts) {
  l <- do.call(layout, list(net.bg)) 
  plot(net.bg, edge.arrow.mode=0, layout=l, main=layout) }

# Resetting graph window to contain just one plot
par(mfrow=c(1,1))

### Network plotting examples
# Curved edges
# 'edge.curved' can be between 0 and 1
plot(g, vertex.label.color="white", edge.curved=0.2,
     edge.arrow.size=0.5, vertex.size=10)

# vertex frames color
plot(g, vertex.label.color="black",
     vertex.frame.color="gray", vertex.color="white",
     edge.curved=0.2, edge.arrow.size=.5)


# Vertex and edge attributes with names including
# "color", "label", "size", etc. are used by the plotting
# function. Instead of specifying an argument to 'plot'
# we can set an appropriate attribute.
#
# Example: Highlight ties involving at least one Female

# default edge color
E(partners)$color <- "gray"

# verify
plot(partners, edge.arrow.size=.5, 
     vertex.label.color="white")

# set attribute 'color' for edges sent by females
female.ids <- V(partners)[gender=="F"]
female.ids
E(partners)[from(female.ids)]$color <- "blue"

plot(partners, edge.arrow.size=.5, 
     vertex.label.color="white", edge.curved=0.2)


## Media graph
# Generate colors base on media type:
colors <- c("gray50", "tomato", "gold")
V(media)$color <- colors[V(media)$media.type]

# Use the audience size value for the node size:
V(media)$size <- V(media)$audience.size*0.6

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(media)$label <- NA

#change arrow size and edge color:
E(media)$arrow.size <- .2
E(media)$edge.color <- "gray80"
E(media)$width <- E(media)$weight
plot(media) 
legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1)

################################################################################
## Introduction to Social Network Statistics in R ##
################################################################################
# Media network
vcount(media) # number of nodes
ecount(media) # number of edges
graph.density(media)  # density (very sparse)

# Enron network
vcount(enron)
ecount(enron)
graph.density(enron, loops=TRUE)

# Kite network
vcount(kite)
ecount(kite)
graph.density(kite, loops=TRUE)


### Degrees and their distribution ###
# vector of degrees (directed graph)
head(degree(enron))  # total degree
head(degree(enron, mode="in"))   # in-degree
head(degree(enron, mode="out"))  # out-degree

summary(degree(enron))

# vector of degrees (undirected graph)
degree(partners)  # total degree
degree(partners, mode="in")   # in-degree
degree(partners, mode="out")  # out-degree
# for undirected graphs the in and out-degrees are the same

summary(degree(partners))

# degree distribution
degtab <- table(degree(partners))
degtab

plot(degtab, main = "Degree Distribution for Partner Network", xlab = "Degree", ylab = "Frequency") # plot degree distribution

barplot(table(degree(kite)), main = "Degree Distribution for Kite Network", xlab = "Degree", ylab = "Frequency", col = "aquamarine")

## Subgraphs and components ###
## Components (clusters)
k <- clusters(enron)
k
# 'k' is a list with:
# membership  =  vector assigning nodes to components
# component size = number of nodes in each component
# number of components  

# "strong" = relationships have to go both ways
k2 <- clusters(enron, "strong")
k2

# Create a new graph for illuatration
k <- g %>%
  add_vertices(4, color = "red") %>%
  add_edges(c(13,14, 11,12, 11,13, 11,14))
plot(k)
k3 <- clusters(k, "strong")

### Paths ###

# Matrix of shortest paths between nodes
shortest.paths(kite)

average.path.length(kite)

### Subgraphs ###

# Create subgraph containing all nodes in the largest
# strongly connected component.
# Using 'induced.subgraph'.

which.max(k3$csize)
i <- which(k3$membership == which.min(k3$csize))
# largest component of the Enron data
k.lc <- induced.subgraph(k, V(k)[i])

plot(k.lc)


### Network diameter ###
# diameter: longest shortest path
# by default directed
diameter(g)
diameter(enron, directed=FALSE)

# get vertex ids of nodes on the longest shortest path
l <- get.diameter(media)
l

# color the edges adjacent to these vertices (color the
# shortest path itself)
E(media)$color <- "gray"
E(media, path=l )$color <- "red"
# color the vertices on the path
V(media)$color <- "lightblue"
V(media)[l]$color <- "red"

plot(media, vertex.size=5, vertex.label=NA,
     edge.width=2, edge.arrow.size=0.5,
     edge.curved=0.5)


### Centrality ###
# vector of betweenness centrality scores
b <- betweenness(g)
b
which.max(betweenness(kite))

# eigenvector centralities
ec <- evcent(g)
str(ec)
ec <- ec$vector
ec
which.max(evcent(kite)$vector)

# Dotchart is a barplot alternative
dotchart(ec)

# closeness
closeness(g)
head(closeness(enron))
which.max(closeness(kite))

### Homophily or assortivity ###
# by degree
assortativity_degree(enron)

# by media type
assortativity_nominal(media, V(media)$media.type, directed=F)

# by audience size
assortativity(media, V(media)$audience.size, directed=F)


### Transitivity ###
transitivity(g)
transitivity(enron)
transitivity(kite)

