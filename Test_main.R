############################################################
###Title: Test Function for Extreme and Average Value Tests
###Specification: Survival Endpoints
############################################################


# For parallel computing option = "op1", use the following packages
#library('doRNG') 
#library('doParallel')

# load the supporting functions #
source("SupportingFunctions.R")

# Usage of fn.test()
fn.test(dataset,c1,c2,type,p,k,n2,a=2,option="op2",cl=2)

# Arguments
##
# dataset		a dataframe that has the following columns - 
## 					Column 1: SubjectID, 
## 					Column 2: Treatment Indicator, Trt=1 treatment,Trt=0 control
##								
## 					Binary & Continuous Endpoints:
##					Column 3: Binary outcome or Continuous outcome
##					
##					Survival Endpoints:
##					Column 3:  Time to event 
##					Column 4:  Event status 
##								EVENT=1 event occured, EVENT=0 event censored 
## 
##					The rest of the columns (from Column #c1 to Column #c2):binary covariates to define cells
##
##
# c1 			the starting column number of the binary covariates in the dataset
##
##
# c2            the ending column number of the binary covariates in the dataset
##
##
# type 			the type of endpoint and must be one of "cts", "bi" or "sv". type="cts" is for continuous
## 				endpoint, type="bi" is for binary endpoint, type="sv" for survival endpoint
##                              
##                              
# p 			the probability each cell will be selected as in the paper (0<p<1)
## 
##
# k 			a positive integer represents the number of subpopulations to be drawn as in the paper
## 
##
# n2 			the number of simulation ran to estimate the null distributions
## 
##
# a 			the required smallest number of cells included in one sub-population, can be any interger greater than
##				0, the default number is 2 (i.e. to form one sub-population, at least 2 cells are required).
##   
# option		indicates the option for step 3 in the fuction, and must either "op1" or "op2". 
## 				option = "op1" is to use Parallel Computing, option = "op2" is to use "for" loop. Detault is option = "op2".
##
# cl			the number of worker clusters created and used if option = "op1" is specified, otherwise, cl will be ignored.
##				The default cl = 2. 



# Examples #
## All the example datasets were generated with both beneficial group and harmed group ##

########################
## Continuous Outcome ##
########################
ctsexampledata=read.csv("ctsexampledata.csv")
ctsexampledata[1:5,]
# Two-sample t test to test overall difference btw Trt and Control
t.test(outcome~Trt,data=ctsexampledata)

set.seed(123)
fn.test(ctsexampledata,4,8,type='cts',0.1,500,100)


####################
## Binary Outcome ##
####################
biexampledata=read.csv("biexampledata.csv")
biexampledata[1:5,]
# Fisher Exact test to test overall difference btw Trt and Control
yt=biexampledata[biexampledata$Trt==1,3] # outcome in the treatment group
yc=biexampledata[biexampledata$Trt==0,3] # outcome in the control group
nt=length(yt)
nc=length(yc)
y.t=sum(yt)
y.c=sum(yc)
n.t=nt-y.t
n.c=nc-y.c
testdata=matrix(c(y.t,y.c,n.t,n.c),nrow=2,dimnames =
       list(c("Trt", "Control"),
            c("Yes", "No")))
fisher.test(testdata)

set.seed(123)
fn.test(biexampledata,4,8,type='bi',0.1,500,250)


######################
## Survival Outcome ##
######################
library('survival')
survivalexampledata=read.csv("survivalexampledata.csv")
survivalexampledata[1:5,]
# Overall Coxph Model to test difference btw Trt and Control
summary(coxph(Surv(TIME,EVENT)~Trt,data=survivalexampledata))

set.seed(123)
fn.test(survivalexampledata,5,9,type='sv',0.1,500,100)

