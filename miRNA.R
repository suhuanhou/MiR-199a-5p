if(T){
  rm(list=ls())
  library(openxlsx)
  path = dirname(rstudioapi::getActiveDocumentContext()$path)  
  setwd(path)
}

#-----------------------------------------------
choose = 1
#-----------------------------------------------
if(choose == 1){
  miRNAdata <- read.csv(paste0(path, "\\GSE41655_new_exprSet.csv"))
  exprSet = miRNAdata
} else { 
  mRNAdata <- read.csv(paste0(path, "\\GSE41657_new_exprSet.csv"))
  exprSet = mRNAdata
}

vs = "Tumor-Normal"
if(T){
  rownames(exprSet)=exprSet[,1] 
  exprSet = exprSet[2:ncol(exprSet)] 
  meta <- read.csv(paste0(path, "\\meta.csv"))
  group_List = meta[which((meta[,1]) %in% colnames(exprSet)),8]  
}

if(F){
  exprSet[exprSet == 0] <- 1 
  exprSet = log2(exprSet + 1) 
}


# name
sample = c(paste0('sample',1:ncol(exprSet)))

group_list = group_List

ids <- read.xlsx(paste0(path, "\\ids.xlsx"))
names(ids)[1:2] <- c("probe_id","symbol")
ids = ids[match(rownames(exprSet),ids$probe_id),]  

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))])
  probes = as.character(tmp)
  print(dim(exprSet))  
  exprSet = exprSet[rownames(exprSet) %in% probes,]
  
  print(dim(exprSet)) 
  rownames(exprSet) = ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}


new_exprSet <- jimmy(exprSet, ids)  
exprSet = new_exprSet  
# write.csv(exprSet, "miRNA_expression.csv")


if(T){  
  group_list = group_List
  library(reshape2)  
  a = exprSet  
  a$"gene"=rownames(a)
  exprSet_L=melt(a,id.vars=c("gene"))
  colnames(exprSet_L)=c('probe','sample','value')
  group_list = exprSet_L$group=rep(group_list,each=nrow(exprSet))
  group_list = group_List
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='')
  
  
  # Define nodePar
  dev.off() 
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  plot(hc) 
  par(mar=c(5,5,5,10))
  # png(paste0(path,'_hclust.png',sep = ""),res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  # PCA
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_List
  

  dev.off()
  # png(paste0(path,'_pca.png',sep = ""),res=100)
  autoplot(prcomp(df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  dev.off()
}

#
# Difference analysis
#

if(T){
  # Grouping matrix
  library(limma)
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=group_List  
  design 
}



# contrast.matrix
contrast.matrix<-makeContrasts(vs,levels = design)
contrast.matrix 

if (T){
  # lmFit+eBayes+topTable
  deg = function(exprSet,design,contrast.matrix){
    #step1 
    fit <- lmFit(exprSet,design)
    #step2 
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    #step3 
    fit2 <- eBayes(fit2) 
    #step4 
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 

    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
}

re = deg(exprSet,design,contrast.matrix)
write.table(re, file = paste0(path,'/DEmiRNA.csv'),row.names=TRUE,col.names=TRUE,sep=",")



#
# Volcano map
#
library(ggplot2)
nrDEG = re; colnames(nrDEG); DEG=nrDEG
logFC_cutoff <- with(DEG,1)

DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))

g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value)), family = 'Arial' )+
  geom_point(size=5, shape = 21, aes(fill = change), color = "white", na.rm = T, position = "jitter") +
  scale_fill_manual(values = c("steelblue", "grey60", "orangered")) +

  geom_hline(yintercept = -log10(0.05), linetype = "dashed",    size = 1,    color = "grey50"  ) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed",    size = 1,    color = "grey50"  ) +
  xlab("log2(Fold Change)") + ylab("-log10(P-value)") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 

    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position="none", 
  ) +
  guides(fill = guide_legend(
    override.aes =
      list(size = 5)
  )
  ) 

print(g)
gse = ifelse(choose==1,"miRNA","mRNA")
ggsave(g,filename = paste0(gse, '_volcano.pdf'))
#save(new_exprSet,group_list,nrDEG,DEG, file=paste0(gse, '_DEG.Rdata'))




