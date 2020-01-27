######## Step 0 #########################################################
## Read in dataset0 - needs to be in data frame format
## The data structure:
## Column 1: SubjectID, 
## Column 2: Treatment Indicator [called Trt in the function: Trt=1 in treatment,Trt=0 control]
## Column 3: outcome *
## * If it is a survival outcome, the Column 3 is the time to event [called TIME in the function]
## and column 4 is censoring status 1=event occured, 0=censored [called EVENT in the function]
## The rest of the columns (from Column c1 to Column c2) are binary covariate to define cells
## c1 and c2 are the starting and ending Column numbers of the binary covariates in the dataset
## type is the type of endpoint - type="cts" for continuous endpoint
##                              - type="bi" for binary endpoint
##                              - type="sv" for survival endpoint
## p is the probability each cell will be selected as in the paper
## k is the number of subpopulations to be drawn as in the paper
## n2 is the number of simulation ran to estimate the null distributions
## a is the required smallest number of cells included in one sub-population, the default number is 2
##   i.e. to form one sub-population, at least 2 cells are required
## Two options are given in step 3, choose one of them to use, detault is option = "op2"
## option = "op1" is Parallel Computing, cl can be specified to be the number of worker clusters
## created and used, the default cl = 2
## option = "op2" is simplw for loop

############# Z-statistics Functions ################

# z-statistics function for continuous endpoint
fnz.cts=function(ncell,dataset,p,a){
	flip0=rbinom(ncell,1,p)
	while (sum(flip0)<a) {flip0=rbinom(ncell,1,p)}
	Icell.c=which(flip0==1)
	dataset.sam=dataset[which(dataset$cellname %in% Icell.c==T),]
	yt=dataset.sam[dataset.sam$Trt==1,3] # outcome in the treatment group
	yc=dataset.sam[dataset.sam$Trt==0,3] # outcome in the control group
	t.test(yt,yc)$statistic

}

# z-statistics function for binary endpoint
fnz.bi=function(ncell,dataset,p,a){
	flip0=rbinom(ncell,1,p)
	while (sum(flip0)<a)  {flip0=rbinom(ncell,1,p)}
	Icell.c=which(flip0==1)
	dataset.sam=dataset[which(dataset$cellname %in% Icell.c==T),]
	yt=dataset.sam[dataset.sam$Trt==1,3] # outcome in the treatment group
	yc=dataset.sam[dataset.sam$Trt==0,3] # outcome in the control group
	nt=length(yt)
	nc=length(yc)
	pt=sum(yt)/nt
	pc=sum(yc)/nc
	pall=sum(c(yt,yc))/(nt+nc)
    se.bi=sqrt(pall*(1-pall)*(1/nt+1/nc))
	(pt-pc)/se.bi

}

# z-statistics function for survival endpoint
fnz.sv=function(ncell,dataset,p,a){
	flip0=rbinom(ncell,1,p)
	while (sum(flip0)<a)  {flip0=rbinom(ncell,1,p)}
	Icell.c=which(flip0==1)
	dataset.sam=dataset[which(dataset$cellname %in% Icell.c==T),]
	-coef(summary(coxph(Surv(TIME,EVENT)~Trt,data=dataset.sam)))[4]

}


fn.test=function(dataset0,c1,c2,type,p,k,n2,a=2,option="op2",cl=2){
	
###### rename column names #####
if(type=="sv"){
	names(dataset0)[1:4]=c("ID","Trt","TIME","EVENT")
} else {
	names(dataset0)[1:3]=c("ID","Trt","outcome")
}



###### Step 1       #######
###### Define cells #######

# Order the dataset by using covariates (Column c1 to c2)
varc=dataset0[,c1:c2] # varc is a dataset only consist of the binary vocariates#
Mdataset=dataset0[do.call(order,varc),]
l.dataset=dim(dataset0)[1]

fn.cell=function(l1){
	a=Mdataset[l1,c1:c2]
	b=Mdataset[(l1+1),c1:c2]
	ifelse(all(a==b)=='TRUE',0,1)
}

I1=apply(as.matrix(1:(l.dataset-1)),1,fn.cell)
I2=c(0,which(I1==1)) ## the last observation in each cell
ncell=length(I2)  #Number of cells
I2.1=I2+1
I2.2=c(I2[-1],l.dataset)
I3=I2.2-I2 # Number of subjects in each cell
cellname=rep(1:ncell,I3)
dataset=data.frame(cbind(Mdataset,cellname))


###### Step 2       #######
###### Calculated observed S(D) and T(D) as shown in the paper #######

# define fn.z first from one of the three z-statistics functions
if(type=="cts"){
	fn.z=fnz.cts
} else if (type=="bi"){
	fn.z=fnz.bi
} else if (type=="sv"){
	fn.z=fnz.sv
} else {
	print("error")
}
	
result1=replicate(k,fn.z(ncell,dataset,p,a))

## The observed S(D) and T(D) ##
## Extreme Value Test
S_D.1=max(result1,na.rm = T)
T_D.1=min(result1,na.rm = T)

## Average Value Test
S_D.2=mean(ifelse(result1>0,result1,0),na.rm = T)
T_D.2=mean(ifelse(result1<0,result1,0),na.rm = T)


###### Step 3       #######
###### Estimate the Distribution of S(D*) and T(D*) under the null hypothesis, as shown in the paper #######

# A function which runs z-statistics function using k selected sub-populations #
# assign fn.z to be one of the functions: fnz.cts, fnz.bi, fnz.sv
fk=function(k){
	# Shuffle the treatment assignment in each cell dataset0$Trt is the Trt variable
	dataset$Trt=unlist(tapply(dataset0$Trt,cellname,sample))
	
	res=replicate(k,fn.z(ncell,dataset,p,a))
	
	c(max(res,na.rm = T),min(res,na.rm = T), mean(ifelse(res>0,res,0),na.rm = T), mean(ifelse(res<0,res,0),na.rm = T))
}

## Under null hypothesis, the distribution of S(D*) and T(D*) ##

if(option=="op1"){
	#Option 1: Parallel Computing
	#cl=40 # the number of worker clusters will be created and used, default is 2
	registerDoParallel(cl) 
	Tstat=foreach(i=1:n2,.combine='rbind', .packages='survival') %dorng% fk(k)

} else if (option=="op2"){
	#Option 2: for loop
	Tstat=NULL

	for (i in 1:n2){
		Tstat0=fk(k)
		Tstat=rbind(Tstat,Tstat0)
	}
} else {
	print("error")
}


###### Step 4       #######
###### Calculate under the null distribution ##########
###### the probability of any value exceeding the observed S(D) and T(D) #######

## Extreme Value Test
I.S.1=round(length(which(Tstat[,1]>S_D.1))/n2,digits=4)
I.T.1=round(length(which(Tstat[,2]<T_D.1))/n2,digits=4)
## Average Value Test
I.S.2=round(length(which(Tstat[,3]>S_D.2))/n2,digits=4)
I.T.2=round(length(which(Tstat[,4]<T_D.2))/n2,digits=4)

TestResults=matrix(c(I.S.1, I.T.1, I.S.2, I.T.2),nrow=2)
colnames(TestResults)=c("Extreme Value Test","Average Value Test")
rownames(TestResults)=c("Benefit","Harm")
TestResults
}


