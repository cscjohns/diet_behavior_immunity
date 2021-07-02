################################################################################
#
# Correlation by Individual Level Product
#
################################################################################

# Format R session -------------------------------------------------------------
# Load packages
library('RColorBrewer')
library('ggplot2')
library('EMMREML')
library('qvalue')
source("data/fwdcilpcodetopost/joaquin_permFDR_function.R")

# Define function
makeTransparent <- function(someColor, alpha=25){
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3],
        alpha = alpha, maxColorValue = 255)})
}

# Load data
load("data/fwdcilpcodetopost/cyno_project_CILP.Rdata")


# Prep data to test for differential correlation between MED and WEST-----------

# scale data within the classes to be compared
tmp <- matrix(nrow = dim(r_matrix)[1], ncol = dim(r_matrix)[2])
for(i in c(1:dim(tmp)[1])){
  tmp[i, which(sample_info$diet == 'WD')] <-
    scale((r_matrix[i, which(sample_info$diet == 'WD')]),
          center = TRUE, scale = TRUE)
  tmp[i, which(sample_info$diet == 'MD')] <-
    scale((r_matrix[i, which(sample_info$diet == 'MD')] ),
          center = TRUE, scale = TRUE)
} 

e_all_norm <- tmp

# CILP--------------------------------------------------------------------------
# Test for differential correlation in LMM framework (focus on top DE genes)

# Create design and Z matrices
design <- model.matrix(~ sample_info$diet)
Z_matrix <- diag(dim(sample_info)[1])

sig_sort <- model_results[order(model_results$diet_pval), ]
sig_to_use1 <- sig_sort$stable_ID[1:140]

# get all possible pairwise combinations of genes
pairs1<-as.data.frame(t(combn( which(model_results$stable_ID %in% sig_to_use1), 2)))

# implement CILP approach using LMM framework
results1 <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs1)[1], function(i) {
  tmp1<-(e_all_norm[pairs1$V1[i],])
  tmp2<-(e_all_norm[pairs1$V2[i],])
  return( c( cor.test( e_all_norm[pairs1$V1[i],] , e_all_norm[pairs1$V2[i],])$estimate, cor.test( e_all_norm[pairs1$V1[i],which(sample_info$diet=='WD')] , e_all_norm[pairs1$V2[i],which(sample_info$diet=='WD')])$estimate, cor.test( e_all_norm[pairs1$V1[i],which(sample_info$diet=='MD')] , e_all_norm[pairs1$V2[i],which(sample_info$diet=='MD')])$estimate, emmreml(y=as.matrix(tmp1*tmp2),X=design,Z=as.matrix(Z_matrix),K=as.matrix(kinship),varbetahat=T,varuhat=T,PEVuhat=T,test=T)$pvalbeta[2,8] ))
}))))

names(results1)<-c('cor','cor_WD','cor_MD','p_value')
results1$gene1<-model_results$stable_ID[pairs1$V1]
results1$gene2<-model_results$stable_ID[pairs1$V2]

# results with non empirical q-value
results1$qvalue<-qvalue(results1$p_value)$qvalues

############
# test for differential correlation in LMM framework (focus on top DE genes)
# permute the sample labels, to get an empirical null
############

# 10 permutations shown here, but this number can be increased (n=100 was used in the main text, but will take awhile to run)
perm_results<-c()

for (i in 1:10){
  sample_info$diet_perm_id<-sample(1:length(sample_info$diet))
  design_perm<-model.matrix(~sample_info$diet[sample_info$diet_perm_id])
  kinship_perm<-kinship[sample_info$diet_perm_id,sample_info$diet_perm_id]
  
  results2 <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs1)[1], function(i) {
    tmp1<-(e_all_norm[pairs1$V1[i],])
    tmp2<-(e_all_norm[pairs1$V2[i],])
    return( c(  emmreml(y=as.matrix(tmp1*tmp2),X=design_perm,Z=as.matrix(Z_matrix),K=as.matrix(kinship),varbetahat=T,varuhat=T,PEVuhat=T,test=T)$pvalbeta[2,8] ))
  }))))
  names(results2)<-c('p_value')
  
  perm_results<-c(perm_results,results2$p_value)
  print(i) }

results2=perm.fdr(results1,matrix(perm_results,ncol=10),'p_value',plot=F,details=T)

############
# plot differential correlations between MED and WEST (permutation results as used in the main text)
############

# read in results
results2=read.delim('data/fwdcilpcodetopost/empirical_FDR_results_LMM.txt')
sig_results<-subset(results2,q_ST_perm<0.2)
not_sig_results<-subset(results2,q_ST_perm>0.2)

# plot changes in correlation in the whole dataset
par(mfrow=c(1,1))
plot(not_sig_results$cor_MD,not_sig_results$cor_WD,pch=20,col=makeTransparent("grey"),bty='n',xlab='Correlation coefficient: MD',ylab='Correlation coefficient: WD',cex.lab=1.25,cex.axis=1.25)
points(sig_results$cor_WD , sig_results$cor_MD,col='steelblue',pch=20)
abline(v=0,lty=2,lwd=2);abline(h=0,lty=2,lwd=2)
legend('topleft',col=c('lightgrey','steelblue'),c('NS',"FDR<0.2"),pch=c(20,20),bty='n')

# plot example with RF00683
tmp<-subset(sig_results,gene1=='ENSMFAG00000030036' | gene2=='ENSMFAG00000030036')
ex<-tmp[which( tmp$p_value==min(tmp$p_value)),]

par(mfrow=c(2,2))
plot( e_all_norm[ which(model_results$stable_ID==ex$gene1),which(sample_info$diet=='WD')] , e_all_norm[ which(model_results$stable_ID==ex$gene2),which(sample_info$diet=='WD')]  ,bty='n',xlab='Normalized expression of RF00283',ylab='Normalized expression of KLF11',bty='n',cex.lab=1.25,cex.axis=1.25,pch=20,col='darkorange',main='Western diet')
abline(lm ( e_all_norm[ which(model_results$stable_ID==ex$gene2),which(sample_info$diet=='WD')] ~ e_all_norm[ which(model_results$stable_ID==ex$gene1),which(sample_info$diet=='WD')]  ),lty=2,lwd=2,col='darkorange')

plot( e_all_norm[ which(model_results$stable_ID==ex$gene1),which(sample_info$diet=='MD')] , e_all_norm[ which(model_results$stable_ID==ex$gene2),which(sample_info$diet=='MD')]  ,bty='n',xlab='Normalized expression of RF00283',ylab='Normalized expression of KLF11',bty='n',cex.lab=1.25,cex.axis=1.25,pch=20,col='steelblue',main='Mediterranean diet')
abline(lm ( e_all_norm[ which(model_results$stable_ID==ex$gene2),which(sample_info$diet=='MD')] ~ e_all_norm[ which(model_results$stable_ID==ex$gene1),which(sample_info$diet=='MD')]  ),lty=2,lwd=2,col='steelblue')

############
# look at # of differential corrrelations a gene is involved in to identify hub genes
############

# number of differential corrrelations each gene in involved in
tmp<-as.data.frame(table(c(as.character(sig_results$gene1),as.character(sig_results$gene2))))
tmp<-merge(tmp,model_results[,1:2],by.x='Var1',by.y='stable_ID',all.x=T)
tmp<-tmp[order(tmp$Freq),]

# list of human TFs from TRRUST (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753191/)
tfs=read.delim('human_TFs.txt',header=F)
tfs$is_tf<-'Yes'
tmp_tf=merge(tmp,tfs,by.x='external_gene_name',by.y='V1',all.x=T)

# generate null distribution of number of differential corrrelations each gene in involved in
sig_dist<-c()
for (i in 1:100){
  sig_tmp<-pairs1[sample(1:dim(pairs1)[1],445),]
  sig_tmp2<-as.data.frame(table(c(sig_tmp$V1,sig_tmp$V2)))
  sig_dist<-c(sig_dist,sig_tmp2$Freq)}

df<-as.data.frame(c(tmp$Freq,sig_dist))
names(df)<-'sig_correlations'
df$Distribution<-'Null'
df$Distribution[1:dim(tmp)[1]]<-'Observed'

par(mfrow=c(1,1))
g1<-ggplot(df, aes(sig_correlations, fill = Distribution)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')+theme_bw(15)+ylab('Density')+xlab('Number of differential correlations per gene') + theme(legend.position = c(0.8, 0.8))+  scale_fill_brewer(palette="Dark2")
g1+ annotate(geom="text", x=41, y=0.025, label="RF00683",size=4)