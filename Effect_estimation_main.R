library(Hmisc)
library(MCMCpack)
library(cubature)
library(mvtnorm)
library(foreach)
library(doParallel)
cl<-makeCluster(40)
registerDoParallel(cl)
library(doRNG)
source("misc.R")  # misc.R has the functions for computations

# subpopulation:   a data framework with 9 columns with each row representing a cell.
#                  The columns are number of non-event in control arm, number of event in control arm,
#                  total number of subjects in control arm, number of non-event in treatment arm, number of event in treatment arm,
#                  total number of subjects in treatment arm, total number of non-event, total number of event, total number of 
#                  subject.
# n: the total number of subjects from which the subpopulation is selected
# prior.population: the population used for the prior estimation with the same data structure as subpopulation
# n.sam: number of samples for importance sampling
# n.iter: number of computation threads


out<-est.main(subpop=subpopulation, prior=prior.population)



