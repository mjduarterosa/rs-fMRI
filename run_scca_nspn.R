# Clean workspace
rm(list=ls())

# Set working directory
setwd('/Users/maria/Documents/NSPN/code/')

# Load data
vars <- as.matrix(read.csv('/Users/maria/Documents/NSPN/code/vars.csv', header=FALSE))
conn <- as.matrix(read.csv('/Users/maria/Documents/NSPN/code/nets.csv', header=FALSE))

# Load packages
library(PMA)
library(sampling)

perm.out <- CCA.permute(x=conn,z=vars,typex="standard", typez="standard",standardize=TRUE, nperms=100)
out <- CCA(conn,vars, typex="standard", typez="standard",penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init,standardize=TRUE) # Save time by inputting  lambda and v

