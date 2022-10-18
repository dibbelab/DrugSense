set.seed(1234)
eff_HB_broad <- readRDS('/mnt/d/data/predicted_efficacy/efficacy_HB_broad_solid.rds')
eff_HB_broad <- eff_HB_broad[-grep(':',eff_HB_broad$drug.id,fixed=T),]

drugID <- unique(eff_HB_broad$drug.id)

eff_HB_broad <- data.frame(split(eff_HB_broad$`G-value`,eff_HB_broad$sample.id))
rownames(eff_HB_broad) <- drugID

HB_sub <- openxlsx::read.xlsx('input/HB_subtypes.xlsx')
cl <- split(HB_sub$id,HB_sub$subtype)

x <- as.matrix(eff_HB_broad)
x_aggr <- apply(x,2,function(y)rownames(x)[order(y)])


res <- lapply(cl,function(i){mu=rownames(x)[order(rowMeans2(x[,i]))][1:20]
Me=rownames(x)[order(rowMedians(x[,i]))][1:20]
mn=rownames(x)[order(rowMins(x[,i]))][1:20]
sp_ce=RankAggreg::RankAggreg(t(x_aggr[,i]),k=20)$top.list
return(list(mu,Me,mn,sp_ce))})

res <- do.call(function(...)mapply(cbind,...,SIMPLIFY = F),args = res)

pdf('barplot.pdf',paper = 'a4')
par(mfrow=c(2,1))
drug_meta <- readRDS('rds/CTRPv2_sens.rds')$drug_meta
v <- drug_meta[match(res[[4]][,2],drug_meta$drugid),]
barplot(sort(table(v$target_or_activity_of_compound),decreasing = T),las=2,col='#e7e1ef',cex.names  =0.6)
v <- v[!v$drugid%in%res[[4]][,1],]
barplot(sort(table(v$target_or_activity_of_compound),decreasing = T),las=2,col='#e7e1ef',cex.names =0.6)
dev.off()

#####

library(clusterProfiler)
library(org.Hs.eg.db)

PCC <- readRDS('rds/PCC_ccle_CTRPv2_filtered.rds')$solid
HB_expr <- readRDS('rds/HB_expr.rds')

comm_gene <- intersect(rownames(PCC),rownames(HB_expr))
CDK_inh <- PCC[comm_gene,c('Alvocidib','Dinaciclib')]
mrk_CDK_inh <- apply(CDK_inh,2,function(x)names(sort(x))[1:250])

comm_mrk <- intersect(mrk_CDK_inh[,1],mrk_CDK_inh[,2])
comm_mrk <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=comm_mrk, column='SYMBOL', keytype='ENSEMBL')


#c2 <- read.gmt('input/c2.all.v7.4.symbols.txt')
c5 <- read.gmt('input/c5.go.v7.4.symbols.txt')
hll <- read.gmt('input/h.all.v7.4.symbols.txt')

#egmt_c2 <- enricher(comm_mrk, TERM2GENE=c2, pAdjustMethod = "BH")
egmt_c5 <- enricher(comm_mrk, TERM2GENE=c5, pAdjustMethod = "BH")
egmt_hll <- enricher(comm_mrk, TERM2GENE=hll, pAdjustMethod = "BH")

#p1 <- dotplot(egmt_c2,orderBy='GeneRatio',showCategory=10) + ggtitle("curated gene sets") +theme(legend.position = 'bottom')
p1 <- dotplot(egmt_c5,orderBy='GeneRatio',showCategory=10) + ggtitle("ontology gene sets") +theme(legend.position = 'bottom')
p2 <- dotplot(egmt_hll,orderBy='GeneRatio',showCategory=10) + ggtitle("hallmark gene sets") +theme(legend.position = 'bottom')

pdf('Supplementary_S3.pdf')
easyGgplot2::ggplot2.multiplot(p1,p2,cols = 1)
dev.off()

