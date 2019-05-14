### mean_mode

P<-read.table("phenotypes_warm_2000.txt",sep="\t",header=T)
G<-read.table("genotype_vcf012_warm_wo_extraNcol.txt",sep="",header=T,comment.char="",check.names=F,as.is=T)

Strain<-names(G)[-c(1:2)] 
P_filtered<-P[P$Strain %in% Strain,]

Plong<-aggregate(P_filtered[3:7],by=list(Strain=P_filtered$Strain,treatment=P_filtered$treatment),FUN=mean,na.rm=T)

Pwide<-reshape(Plong,idvar="Strain",direction="wide",timevar="treatment",v.names=c("mean_max","mean_mode","mean_shape"))

GT<-data.frame(t(G[,-c(1:2)])) # bara genotyper och transponerade
names(GT)<-paste(G[,1],G[,2],sep="_") # skapa SNP-namn genom att slå ihop contig och position
GT[GT==-1]<-NA
nsnp<-ncol(GT)

rare<-rep(FALSE,nsnp)
for (i in 1:nsnp) {
	tb<-table(GT[,i]) # counts of 0, 1, 2 for sample i
	if (length(tb)<=1) rare[i]<-TRUE
	else if (length(tb)==2 & min(tb)<2) rare[i]<-TRUE
}

GT_rare<-GT[,!rare]
nsnp<-ncol(GT_rare)
Strain<-names(G)[-c(1:2)]
snpdata<-cbind(Strain,GT_rare) # join sample name column with GT_rare
M<-merge(Pwide,snpdata)

snps<-colnames(M)[7:19356]
nstrains<-length(Strain)

### Linear model

p_mean_mode_lm<-rep(NA,nsnp)
p_mean_mode<-rep(NA,nsnp)
p_mean_mode_slope<-rep(NA,nsnp)
for (i in 1:nsnp){
	p_mean_mode_lm<-summary(lm(M[,5]~M[,6+i]))
	p_mean_mode[i]<-p_mean_mode_lm$coefficients[2,4]
	p_mean_mode_slope[i]<-p_mean_mode_lm$coefficients[2,1]
}
 
ntop<-20 # use the top ntop results of the regressions only

p_mean_mode_sorted<-sort(p_mean_mode)[1:ntop]
p_mean_mode_sorted_ind<-sort(p_mean_mode,index.return=TRUE)$ix[1:ntop]

### Permutation test

pind<-rep(1:nstrains) # permuted indices
nrep<-1000000 # number of permutations

p_mean_mode_pperm_top<-rep(0,ntop)
for(i in 1:ntop){
   myInd_mean_mode<-p_mean_mode_sorted_ind[i]
   mySlope_mean_mode<-p_mean_mode_slope[myInd_mean_mode] # top ntop slopes
   permSlope_mean_mode<-rep(NA,nrep)
   for(k in 1:nrep){
      curPerm<-sample(pind,nstrains,replace=FALSE) # shuffle the sample index
      permSlope_mean_mode[k]<-summary(lm(M[,5]~M[curPerm,6+myInd_mean_mode]))$coefficients[2,1]
   }
   pos_mean_mode<-match(mySlope_mean_mode,sort(append(permSlope_mean_mode,mySlope_mean_mode))) # find position among perm sample
   p_mean_mode_pperm_top[i]<-pos_mean_mode/(nrep+1)
}


results_mean_mode<-cbind(p_mean_mode_pperm_top,snps[p_mean_mode_sorted_ind])
write("mean_mode adjusted\n",file="perm_results_mean_mode.txt")
write.table(results_mean_mode,file="perm_results_mean_mode",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

