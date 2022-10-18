get_subcat_perc <- function(i,d_meta,g_sets,topG){
  d.id <- as.character(d_meta[i,'drugid'])
  d.target <- as.character(d_meta[i,'gene_symbol_of_protein_target'])
  gs_name <- g_sets[g_sets$gene_symbol%in%d.target,c('gs_subcat','gs_name')]
  
  df <- g_sets[g_sets$gs_name%in%gs_name[[2]] & g_sets$ensembl_gene%in%topG[,d.id],c('gs_subcat','gs_name')]
  df <- df[!duplicated(df),]
  
  tbl_target <- table(gs_name$gs_subcat)
  tbl_markers <- table(df$gs_subcat)
  
  missing <- setdiff(names(tbl_target),names(tbl_markers))
  
  if(length(missing)>0) tbl_markers[missing] <- 0
  
  return(round(tbl_markers[names(tbl_target)]/tbl_target*100,1))
}



all_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens")
all_gene_sets <- all_gene_sets[all_gene_sets$gs_cat=='C5',]

gs <- split(all_gene_sets$ensembl_gene,all_gene_sets$gs_name)

pcc <- readRDS('rds/PCC_ccle_CTRPv2_filtered.rds')

gs_solid <- sapply(gs,function(x)intersect(x,rownames(pcc$solid)))
gs_liquid <- sapply(gs,function(x)intersect(x,rownames(pcc$liquid)))

#####

d_meta <- readRDS('rds/CTRPv2_sens.rds')$drug_meta

l1 <- !grepl(':',d_meta$drugid)
l2 <- d_meta$drugid %in% colnames(pcc$solid)
l2.2 <- d_meta$drugid %in% colnames(pcc$liquid)
l3 <- !d_meta$gene_symbol_of_protein_target==''
l4 <- !grepl(';',d_meta$gene_symbol_of_protein_target)

d_meta_solid <- d_meta[l1&l2&l3&l4,]
d_meta_liquid <- d_meta[l1&l2.2&l3&l4,]

pcc_solid <- pcc$solid[,d_meta_solid$drugid]
pcc_liquid <- pcc$liquid[,d_meta_liquid$drugid]

res_solid <- apply(pcc_solid,2,function(x)fgseaSimple(pathways = gs_solid,stats = na.omit(x),nperm = 1000,minSize = 10,maxSize = 500,nproc = 20),simplify = F)
res_liquid <- apply(pcc_liquid,2,function(x)fgseaSimple(pathways = gs_liquid,stats = na.omit(x),nperm = 1000,minSize = 10,maxSize = 500,nproc = 20),simplify = F)

#####

res_solid <- readRDS('rds/res_solid.rds')
res_liquid <- readRDS('rds/res_liquid.rds')

ptw_solid <- sapply(res_solid,function(x)x$pathway)
NES_solid <- sapply(res_solid,function(x)x$NES)

ptw_liquid <- sapply(res_liquid,function(x)x$pathway)
NES_liquid <- sapply(res_liquid,function(x)x$NES)

rm(res_solid,res_liquid); gc()

trg_solid <- lapply(d_meta_solid$gene_symbol_of_protein_target,function(i)all_gene_sets$gs_name[all_gene_sets$gene_symbol==i])
trg_liquid <- lapply(d_meta_liquid$gene_symbol_of_protein_target,function(i)all_gene_sets$gs_name[all_gene_sets$gene_symbol==i])

logi_solid <- mapply(function(x,y)x%in%y,ptw_solid,trg_solid,SIMPLIFY = F)
logi_liquid <- mapply(function(x,y)x%in%y,ptw_liquid,trg_liquid,SIMPLIFY = F)

NES_splt_solid <- mapply(function(x,i)split(x,i),NES_solid,logi_solid)
NES_splt_liquid <- mapply(function(x,i)split(x,i),NES_liquid,logi_liquid)

df <- list(solid=NES_splt_solid,liquid=NES_splt_liquid)
df <- reshape2::melt(df)

pdf('boxplot.pdf',width = 2.5,height = 5,fonts = 'Arial')
ggplot(df,aes(x=L1,y=value,fill=L3)) + geom_boxplot(position=position_dodge(.9),outlier.shape = NA) + ylab('NES') +xlab('') + theme_minimal() +
  coord_flip() + geom_hline(yintercept = 0, size = 1,linetype = "dashed",color='black') + scale_fill_manual(values = c('#b2182b','azure4')) + theme(legend.position = c(0.95, 0.8))
dev.off()

#####

topG_solid <- apply(pcc_solid,2,function(x)names(sort(x))[1:250])
topG_liquid <- apply(pcc_liquid,2,function(x)names(sort(x))[1:250])

sapply(1:ncol(topG_solid),function(i){
  l1 <- all_gene_sets$gs_name%in%trg_solid[[i]]&
  all_gene_sets$gs_name[]
})


df <- data.frame(M=c(mean(),mean(df$value[df$L1=='solid'&df$L3=='FALSE'],na.rm=T),
mean(df$value[df$L1=='liquid'&df$L3=='TRUE']),mean(df$value[df$L1=='liquid'&df$L3=='FALSE'],na.rm=T)),

SD=c(sd(df$value[df$L1=='solid'&df$L3=='TRUE']),sd(df$value[df$L1=='solid'&df$L3=='FALSE'],na.rm=T),
  sd(df$value[df$L1=='liquid'&df$L3=='TRUE']),sd(df$value[df$L1=='liquid'&df$L3=='FALSE'],na.rm=T)),

L1 =rep(c('solid','liquid'),each=2),

L2 =rep(c('TRUE','FALSE'),2))

ggplot(df, aes(x=L1, y=M, fill=L2)) +
  geom_bar(position=position_dodge(.9), stat="identity",
           colour='black') +
  geom_errorbar(aes(ymin=M-SD, ymax=M+SD), width=.2)

wilcox.test(df$value[df$L1=='liquid'&df$L3=='TRUE'],df$value[df$L1=='liquid'&df$L3=='FALSE'],alternative='less')

all_gene_sets <- all_gene_sets[all_gene_sets$ensembl_gene %in% union(rownames(pcc_solid),rownames(pcc_liquid)),]
tb <- table(all_gene_sets$gs_name)
all_gene_sets <- all_gene_sets[all_gene_sets$gs_name%in%names(which(tb >=10 & tb<=500)),]
logi <- all_gene_sets$gs_subcat==''
all_gene_sets[logi,"gs_subcat"]<-all_gene_sets[logi,"gs_cat"]


topG_random <- replicate(1000,replicate(ncol(topG_solid),sample(rownames(pcc_solid),50)))
colnames(topG_random) <- colnames(topG_solid)

perc_solid <- lapply(1:nrow(d_meta_solid),function(i)get_subcat_perc(i,d_meta_solid,all_gene_sets,topG_solid))
perc_liquid <- lapply(1:nrow(d_meta_liquid),function(i)get_subcat_perc(i,d_meta_liquid,all_gene_sets,topG_liquid))

cl <- parallel::makeCluster(20)
parallel::clusterExport(cl,c('topG_random','d_meta_solid','all_gene_sets','get_subcat_perc'))

perc_random <- parallel::parLapply(cl,1:1000,function(i){
  perc <- lapply(1:nrow(d_meta_solid),function(j)get_subcat_perc(j,d_meta_solid,all_gene_sets,topG_random[,,i]))
})

parallel::stopCluster(cl)

