#######TSRP(Transcriptome-based Stemness Regulator Prediction)##
#######Author:Jihong Yang#######################################
#######Email: yangjiho123@gmail.com#############################


#####step1-1
options(stringsAsFactors = FALSE)
DR<-read.table("Whole_data.txt",sep = "\t", header = TRUE) #
s <- read.table("stemness21markers.txt",sep = "\t", header = FALSE) #
DR1 <- DR[,2:dim(DR)[2]]
cormat <- round(cor(DR1),2)

cindex1 = colnames(cormat)
temp = matrix(0,1,length(s$V1))
for (i in 1:length(s$V1)){
  temp[1,i]  = which(cindex1==s$V1[i])
}
diag(cormat) = 0
cormat1 <- cormat[,temp]

write.table(cormat1, file="output/cormat1.txt", 
            row.names=TRUE, col.names=TRUE, sep = "\t")
            
#####step1-2
cluster1 = c (10,13,15,18,20)
cluster2 = c(2,3,4,14,17,19,21)
cluster3 = c(1,5,6,7,8,9,11,12,16)
colnames(cormat1)[cluster1]
colnames(cormat1)[cluster2]
colnames(cormat1)[cluster3]

index_21 = match(s$V,rownames(cormat1))
cormat2 = cormat1[-index_21,]
#match(s$V,rownames(cormat2))
#dim(cormat1)
#dim(cormat2)
c1 = as.data.frame(rowSums(cormat2[,cluster1]))/length(cluster1)
c2 = as.data.frame(rowSums(cormat2[,cluster2]))/length(cluster2)
c3 = as.data.frame(rowSums(cormat2[,cluster3]))/length(cluster3)

a = 1.697
b=1.255
c=-0.378
StemScore = a*c1+b*c2+c*c3
colnames(StemScore) = "StemScore1"
pos_index = which(StemScore$StemScore1>0)
neg_index = which(StemScore$StemScore1<0)
#length(pos_index)+length(neg_index)

pos = StemScore[pos_index,]
pos_name = rownames(StemScore)[pos_index]
neg = StemScore[neg_index,]
neg_name = rownames(StemScore)[neg_index]

z_pos = as.data.frame((pos-mean(pos))/sd(pos))
colnames(z_pos) = "pscore"
z_neg = as.data.frame((neg-mean(neg))/sd(neg))
colnames(z_neg) = "nscore"

index_cluster1 = which(z_pos$pscore> 1.65)
index_cluster2_1 = which(z_pos$pscore < -1.65)
index_cluster2_2 = which(z_neg$nscore> 1.65)
index_cluster3 = which(z_neg$nscore < -1.65)

cluster1_like_genes = pos_name[index_cluster1]
cluster2_like_genes = union(pos_name[index_cluster2_1],neg_name[index_cluster2_2])
cluster3_like_genes = neg_name[index_cluster3]
write.table(cluster1_like_genes, file="output/cluster1_like_genes.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")
write.table(cluster2_like_genes, file="output/cluster2_like_genes.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")
write.table(cluster3_like_genes, file="output/cluster3_like_genes.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")

#####step2-1
# code for overlap analysis
options(stringsAsFactors = FALSE)
Sig_up <- read.table("step2_up.txt",sep = "\t", header = FALSE) #
Sig_dw <- read.table("step2_dn.txt",sep = "\t", header = FALSE) #

coregenes <- read.table("step1_21coregenes.txt",sep = "\t", header = FALSE) #
#allgenes <- read.table("step1_allstemgenes.txt",sep = "\t", header = FALSE) #
length = dim(Sig_dw)[1]
M_up <- matrix(0,length,4)
M_dw <- matrix(0,length,4)
M_up_core <- M_up
M_dw_core <- M_dw
Ncore <- dim(coregenes)[1]
#Nall <- dim(allgenes)[1]
N_background = 8706
######for core genes######
for (i in 1:length){
  temp_up <- Sig_up[i,1]
  tempup1 <- as.data.frame(strsplit(temp_up, ";")[[1]])
  temp_dw <- Sig_dw[i,1]
  tempdw1 <- as.data.frame(strsplit(temp_dw, ";")[[1]])
  temp_ov_core_up <- intersect(tempup1$`strsplit(temp_up, ";")[[1]]`,coregenes$V1)
  temp_ov_core_dw <- intersect(tempdw1$`strsplit(temp_dw, ";")[[1]]`,coregenes$V1)
  #temp_ov_all_up <- intersect(tempup1$`strsplit(temp_up, ";")[[1]]`,allgenes$V1)
  #temp_ov_all_dw <- intersect(tempdw1$`strsplit(temp_dw, ";")[[1]]`,allgenes$V1)
  #M_up[i,1] = length(temp_ov_all_up)
  #M_up[i,2] = dim(tempup1)[1]
  #M_up[i,3] = Nall
  #M_up[i,4] = N_background
  
  #M_dw[i,1] = length(temp_ov_all_dw)
  #M_dw[i,2] = dim(tempdw1)[1]
  #M_dw[i,3] = Nall
  #M_dw[i,4] = N_background
  
  M_up_core[i,1] = length(temp_ov_core_up)
  M_up_core[i,2] = dim(tempup1)[1]
  M_up_core[i,3] = Ncore
  M_up_core[i,4] = N_background
  
  M_dw_core[i,1] = length(temp_ov_core_dw)
  M_dw_core[i,2] = dim(tempdw1)[1]
  M_dw_core[i,3] = Ncore
  M_dw_core[i,4] = N_background
}

#write.table(M_up, file="M_up.txt", 
            #row.names=FALSE, col.names=FALSE, sep = "\t")
#write.table(M_dw, file="M_dw.txt", 
            #row.names=FALSE, col.names=FALSE, sep = "\t")
write.table(M_up_core, file="output/M_up_core.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")
write.table(M_dw_core, file="output/M_dw_core.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")
#####then combined with "ForStep2_treatmentID.txt"
#####p-value was calculated in EXCEL using HYPGEOM.DIST
####p-value was then adjusted in Prism (version 9.0) to obtain q-value
##### ratio was also calculated in EXCEL
####the whole result tables were "Result_step2_up.txt" and "Result_step2_dw.txt"

###step3-1 code for overlap analysis for combination
options(stringsAsFactors = FALSE)
Sig_dw <- read.table("step2_dn.txt",sep = "\t", header = FALSE) #
coregenes <- read.table("step1_21coregenes.txt",sep = "\t", header = FALSE) #
Index <- read.table("Index_treatmentID.txt",sep = "\t", header = FALSE) #
Index_Treatment <- read.table("id_treatment.txt",sep = "\t", header = FALSE) #
length = dim(Index)[1]
N_all = length*(length-1)/2
M_dw_core_ov_genes_two <- matrix(0,N_all,9)
N_background = 8706
Ncore <- dim(coregenes)[1]
Narray <- combn(length, 2)
######for core genes######
for (i in 1:N_all){
  M_dw_core_ov_genes_two[i,1] = Index_Treatment$V2[Narray[1,i]] # treatment 1
  M_dw_core_ov_genes_two[i,2] = Index_Treatment$V2[Narray[2,i]] # treatment 2
  Index1 <- Index_Treatment$V1[Narray[1,i]]
  Index2 <- Index_Treatment$V1[Narray[2,i]]
  temp1 <- Sig_dw[Index1,1]
  tempindex1 <-  as.data.frame(strsplit(temp1, ";")[[1]])
  temp2 <- Sig_dw[Index2,1]
  tempindex2 <-  as.data.frame(strsplit(temp2, ";")[[1]])
  temp_union <- union(tempindex1$`strsplit(temp1, ";")[[1]]`, tempindex2$`strsplit(temp2, ";")[[1]]`)
  temp_ov <- intersect(temp_union,coregenes$V1)
  temp_ov1 <- intersect(tempindex1$`strsplit(temp1, ";")[[1]]`,coregenes$V1)
  temp_ov2 <- intersect(tempindex2$`strsplit(temp2, ";")[[1]]`,coregenes$V1)
  M_dw_core_ov_genes_two[i,3] = length(temp_ov)
  M_dw_core_ov_genes_two[i,4] = length(temp_union)
  M_dw_core_ov_genes_two[i,5] = Ncore
  M_dw_core_ov_genes_two[i,6] = N_background
  M_dw_core_ov_genes_two[i,7] = paste(temp_ov, collapse = "; ")
  M_dw_core_ov_genes_two[i,8] = paste(temp_ov1, collapse = "; ")
  M_dw_core_ov_genes_two[i,9] = paste(temp_ov2, collapse = "; ")
}

write.table(M_dw_core_ov_genes_two, file="output/CombinationsResults.txt", 
            row.names=FALSE, col.names=FALSE, sep = "\t")
#####p-value was calculated in EXCEL using HYPGEOM.DIST
####the whole result tables were "Result_step3_Combinations.txt"


