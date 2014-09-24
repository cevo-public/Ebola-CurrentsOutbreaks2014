library(ape)
library(TreePar)
setwd("/Users/tstadler/Documents/Data/Uni/Research/R/ebola")
source("ebolaFunctions.R")

############################################
############################################
# read data trees

setwd("/Users/tstadler/Documents/Data/Uni/Research/Datasets0613/Ebola/SciencePaperDataRambaut/only_SL_1Aug/Informed_prior")
trees<-read.nexus("2014_SL.HKY_strict.exp.trees")
trees<-trees[1001:10001]
trees<-trees[100*(1:90)]   


############################################
############################################
# SKYLINE analysis

# estimate death rate
deathfix<-0
# fix sampling probability to 0.7
sprobc <- 0.7
# 3 Re intervals (one prior to most ancestral sample, then 50/50 length) 
numbRe<-3  # alternative numbRe<-1
#calculate likelihood for first 2 trees (we did it for all 90)
numbTree<-2 # full analysis with numbTree<-90

likdif<-vector()
estimates<-c()
for (index in 1:numbTree){
test<-trees[[index]]
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
test<-drop.tip(test,droptip)
rootheight<-max(getx(test,sersampling=1)[,1])
x<-getx(test,sersampling=1)
times<-x[,1]
ttype<-x[,2]
if (numbRe ==1 ){
out<-optim(c(2,1,0.7),LikShiftsSTTebolaConst,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)} else {
out<-optim(c(rep(2,numbRe),1,0.7),LikShiftsSTTebola,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)}
estimates<-rbind(estimates,parepi(out,sprobc=sprobc,deathfix=deathfix))
}
estimatesraw<-estimates

estimatesmedian<-vector()
estimatesHPD<-vector()
estimatesmean<-vector()
estimatesvar<-vector()
for (i in 1:length(estimates[1,])){
	estimatesmedian<-c(estimatesmedian,median(estimates[,i]))
	estimatesHPD<-cbind(estimatesHPD,HPD(estimates[,i]))
	estimatesmean<-c(estimatesmean,mean(estimates[,i]))
	estimatesvar<-c(estimatesvar,var(estimates[,i]))
}
estimates<-round(rbind(estimatesmedian,estimatesHPD,estimatesmean,estimatesvar,estimatesraw),2)
print(estimates)
#                [,1] [,2] [,3] [,4] [,5]
#estimatesmedian 1.82 1.60 1.15 5.63  0.7
#(lower95)       1.79 1.33 0.68 4.29  0.7
#(upper95)       1.79 1.33 0.68 4.29  0.7
#estimatesmean   1.82 1.60 1.15 5.63  0.7
#estimatesvar    0.00 0.15 0.45 3.61  0.0
#(Tree1)         1.79 1.33 0.68 4.29  0.7
#(Tree2)         1.85 1.87 1.63 6.97  0.7


############################################
############################################
# EI model analysis
# this analysis takes longish (we did each tree on a seperate node on the cluster)

# here we analyse tree 1 with fixed sampling proportion of .7
i<-1
sprob<-0.7
phylo<-trees[[i]]
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
phylo<-drop.tip(phylo,droptip)

par<-c(1/5,2*1/5,1/5)
phylo$edge.length <- phylo$edge.length*365
out<-optim(c(par),LikTypesSTTebolaS,phylo=phylo,sprob=sprob)$par
epi<-c(out[2]/out[3],1/out[1],1/out[3],sprob)
# epi for tree 1 is 1.338087 2.215416 1.376223 0.700000


############################################
############################################
# superspreader model analysis
# this analysis takes longish (we did each tree on a seperate node on the cluster)

# here we analyse tree 1 with fixed sampling proportion of .7
i<-1
sprob<-0.7
phylo<-trees[[i]]
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
phylo<-drop.tip(phylo,droptip)
phylo$edge.length <- phylo$edge.length*365

test<-addroot(phylo,(1/24))
x<-getx(test,sersampling=1)
times<-x[,1]
ttype<-x[,2]
out2<-optim(c(2,1,0.7),LikShiftsSTTebolaConst,times=times,ttype=ttype,sprobc=sprob,deathfix=0,root=0,survival=0)
out<-optim(c(out2$par[1]/2,out2$par[1]/2,out2$par[1]/2,out2$par[2]),LikTypesSTTebolaSSpread,phylo=phylo,sprob=sprob)

# for tree 1:
# out2$par 0.3426510 0.2595218 2.6061846
# out2$value 318.2714
# out$par 0.18973391 0.75521212 0.03998312 0.23648649
# out$value 315.7676