library(magrittr)

get.entropy <- function(mtx){
  e = apply(mtx,1,function(x) entropy::entropy.empirical(entropy::discretize(x[!is.na(x)],numBins = 10)) )
  return(e)
}

filter_rows <- function(mtx,check_zeros=T){
  if(check_zeros) mtx<-mtx[!rowSums(mtx==0)/ncol(mtx)==1,]
  
  e<-get.entropy(mtx)
  Qe <- quantile(e,probs=0.05)

  mtx <- mtx[e>Qe,]

  mu <- rowMeans(mtx,na.rm = T)
  Qmu <- quantile(mu,probs=0.25)

  mtx <- mtx[mu>Qmu,]
  
  mtx[mtx==0] <- NA
  
  return(mtx)
}

hgu219probe2ensid <- function(mtx){
  
  mapped_probes <- AnnotationDbi::mappedkeys(hgu219.db::hgu219ENSEMBL)
  
  xx <- as.list(hgu219.db::hgu219ENSEMBL[mapped_probes])
  mtx <- mtx[rownames(mtx)%in%names(xx),]
  
  ensbl<-xx[match(rownames(mtx),names(xx))]
  mtx <- mtx[rep(rownames(mtx),sapply(ensbl,length)),]
  
  df <- data.frame(ENSEMBL=unlist(ensbl),mtx,stringsAsFactors = F,check.names = F)
  MAD<-apply(df[,2:ncol(df)],1,mad)
  
  df <- df[order(MAD,decreasing = T),]
  df <- df[!duplicated(df$ENSEMBL),]
  
  rownames(df) <- df$ENSEMBL
  
  return(as.matrix(df[,-1]))
}

readAffy_rma_convertId <- function(meta){
  expr <- affy::ReadAffy(filenames = meta$Array.Data.File,celfile.path = 'Expression Data/CEL_CCL/')
  expr <- Biobase::exprs(affy::rma(expr))
  expr <- hgu219probe2ensid(expr)
  colnames(expr) <- meta$COSMIC_ID[match(colnames(expr),meta$Array.Data.File)]
  return(expr)
}

import_CCLE_expression_data <- function(){
  
  file <- "Expression data/CCLE_RNAseq_reads.csv" #21Q4
  expr <- read.csv(file,check.names = F,stringsAsFactors = F)
  colnames(expr)[1] = "DepMap_ID"
  expr <- expr %>% tibble::column_to_rownames("DepMap_ID") %>% t()
  expr <- expr[-grep('ERCC-',rownames(expr)),]
  rownames(expr) <- stringr::str_extract(rownames(expr),'ENSG[0-9]+')
  
  file <- "Expression data/sample_info.csv"
  cell_meta <- read.csv(file,check.names = F,stringsAsFactors = F)
  cell_meta <- cell_meta[match(colnames(expr),cell_meta$DepMap_ID),]
  
  logi_eng_unkwn <- grepl('engineered|unknown',cell_meta$lineage,ignore.case = T)
  cell_meta <- cell_meta[!logi_eng_unkwn,]
  expr <- expr[,!logi_eng_unkwn]
  colnames(expr) <- cell_meta$stripped_cell_line_name
  
  logi_liquid <- cell_meta$sample_collection_site=='haematopoietic_and_lymphoid_tissue'
  expr_list <- list(solid=expr[,!logi_liquid],liquid=expr[,logi_liquid])
  expr_cpm_filtered <- lapply(expr_list,function(x)filter_rows(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=x),normalized.lib.sizes = T))))
  
  return(list(expr = expr_cpm_filtered, cell_meta = cell_meta))    
}

import_GDSC_expression_data <- function(){

  file <- 'Expression Data/E-MTAB-3610.sdrf.txt'
  sdrf <- read.csv(file,sep='\t',stringsAsFactors = F)
  sdrf$COSMIC_ID<- paste('COSM:',sapply(strsplit(sdrf$Source.Name,'_'),tail,n=1),sep='')
  sdrf <- sdrf[-grepl('unknown',sdrf$Extract.Name,ignore.case = T),]
  
  sdrf_splt <- split(sdrf,grepl('haematopoietic|blood',sdrf$Extract.Name,ignore.case = T))
  names(sdrf_splt) <- c('solid','liquid')
  expr_rma <- sapply(sdrf_splt,readAffy_rma_convertId)
  
  return(list(expr = expr_rma, cell_meta = sdrf))
}

import_HB_expression_data <- function(){
  
  HB <- affy::ReadAffy(celfile.path = 'Expression Data/CEL_HB/')
  HB <- Biobase::exprs(affy::rma(HB))
  colnames(HB) <- gsub('_U133A.CEL','',colnames(HB))
  
  CCL <- affy::ReadAffy(celfile.path = 'Expression Data/CEL_CCL/')
  CCL <- Biobase::exprs(affy::rma(CCL))
  
  U133plus_U219 <- read.csv('Expression Data/U133PlusVsU219_BestMatch.txt',sep='\t',check.names = F)
  U133plus_U219 <- U133plus_U219[U133plus_U219$`A Probe Set Name`%in%rownames(HB),]
  
  HB <- HB[U133plus_U219$`A Probe Set Name`,]
  HB <- aggregate(HB,by=list(U133plus_U219$`B Probe Set Name`),FUN=mean, na.rm=TRUE)
  
  rownames(HB) <- HB$Group.1
  HB <- HB[,-1]
  
  mtx <- as.matrix(cbind(HB,CCL[rownames(HB),]))
  mtx <- preprocessCore::normalize.quantiles(mtx,copy = F)
  
  btch<- rep(1,ncol(mtx))
  btch[1:29] <- 0
  
  mtx <- limma::removeBatchEffect(mtx, batch=factor(btch))
  mtx <- hgu219probe2ensid(mtx)
  
  return(mtx[,1:25])
}

import_TCGA_expression_data <- function(){
  
  expr.mtx <- NULL
  
  fls <- grep('assay',list.files('GDC_rsd'),value = T)
  
  for(filename in fls){
    path <- paste('GDC_rsd/',filename,sep='')
    i<- gsub('-assay.rds','',filename)
    expr.profile <- assays(readRDS(path))@listData$`HTSeq - Counts`
    colnames(expr.profile) <- paste(substr(colnames(expr.profile), 1, 12),i,sep='-')
    if(!is.null(expr.mtx))expr.mtx <- cbind(expr.mtx,expr.profile[rownames(expr.mtx),])
    else expr.mtx <- cbind(expr.mtx,expr.profile)
  }
  
  expr.mtx <- filter_rows(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=expr.mtx),normalized.lib.sizes = T)))
  
  return(expr.mtx)
}

import_PPTC_expression_data <- function(){
  file <- "Expression data/tumormap-exp-2019-07-30-pdx.tsv"
  expr <- read.table(file,sep='\t',check.names = F,header = T)
  rownames(expr) <- gsub('\\..*','',expr[[1]])
  expr <- expr[-1]
  
  file <- "Expression data/PCAT_meta.csv"
  meta <-  read.csv(file,stringsAsFactors = F)
  has_preclinical <- meta$sample_id[meta$HasDrugTestingData=='Yes']
  sel_ids <- sort(intersect(names(expr),has_preclinical))
  
  expr <- as.matrix(expr[,sel_ids])
  logi_liquid <- grepl('ALL',meta$Histology_Detailed[match(sel_ids,meta$sample_id)])
  expr_list <- list(solid=expr[,!logi_liquid],liquid=expr[,logi_liquid])
  
  return(list(expr = expr_list, meta = meta))
}

import_PPTC_expression_data2 <- function(){
  file <- "Expression data/tumormap-exp-2019-07-30-pdx.tsv"
  expr_PPTC <- read.table(file,sep='\t',check.names = F,header = T)
  rownames(expr_PPTC) <- gsub('\\..*','',expr_PPTC[[1]])
  expr_PPTC <- expr_PPTC[-1]
  
  file <- "Expression data/CCLE_RNAseq_reads.csv" #21Q4
  expr_CCLE <- read.csv(file,check.names = F,stringsAsFactors = F)
  colnames(expr_CCLE)[1] = "DepMap_ID"
  expr_CCLE <- expr_CCLE %>% tibble::column_to_rownames("DepMap_ID") %>% t()
  expr_CCLE <- expr_CCLE[-grep('ERCC-',rownames(expr_CCLE)),]
  rownames(expr_CCLE) <- stringr::str_extract(rownames(expr_CCLE),'ENSG[0-9]+')
  expr_CCLE <- filter_rows(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=expr_CCLE),normalized.lib.sizes = T)))
  
  comm_genes <- intersect(rownames(expr_PPTC),rownames(expr_CCLE))
  
  expr <- as.matrix(cbind(expr_PPTC[comm_genes,],expr_CCLE[comm_genes,]))
  expr <- preprocessCore::normalize.quantiles(expr,copy = F)
  
  btch<- rep(1,ncol(expr))
  btch[1:ncol(expr_PPTC)] <- 0
  
  expr <- limma::removeBatchEffect(expr, batch=factor(btch))
  expr <- expr[,1:ncol(expr_PPTC)]
  
  file <- "Expression data/PCAT_meta.csv"
  meta <-  read.csv(file,stringsAsFactors = F)
  has_preclinical <- meta$sample_id[meta$HasDrugTestingData=='Yes']
  sel_ids <- sort(intersect(colnames(expr),has_preclinical))
  
  expr <- expr[,sel_ids]
  logi_liquid <- grepl('ALL',meta$Histology_Detailed[match(sel_ids,meta$sample_id)])
  expr_list <- list(solid=expr[,!logi_liquid],liquid=expr[,logi_liquid])
  
  return(list(expr = expr_list, meta = meta))
}

import_CTRPv2_sensitivity <- function(){
  
  CTRPv2 <- PharmacoGx::downloadPSet("CTRPv2_2015",saveDir=file.path("Drug Sensitivity Data/CTRPv2"))
  sensitivityData <- PharmacoGx::summarizeSensitivityProfiles(CTRPv2,sensitivity.measure = "aac_published")
  
  cell_meta <- PharmacoGx::cellInfo(CTRPv2) %>% tibble::as_tibble()
  colnames(sensitivityData) <- cell_meta$ccl_name
  
  logi_liquid <- cell_meta$ccle_primary_site=='haematopoietic_and_lymphoid_tissue'
  sensitivityList <- list(solid=sensitivityData[,!logi_liquid],liquid=sensitivityData[,logi_liquid])
  
  drug_meta <- PharmacoGx::drugInfo(CTRPv2) %>% tibble::as_tibble()
  
  return(list(sens=sensitivityList,cell_meta=cell_meta,drug_meta=drug_meta))    
}

import_GDSCv1_sensitivity <- function(){
  
  GDSCv1 <- PharmacoGx::downloadPSet("GDSC_2020(v1-8.2)",saveDir=file.path("Drug Sensitivity Data/GDSC/"))
  
  file <- 'Drug Sensitivity Data/GDSC1_fitted_dose_response_25Feb20.xlsx'
  tbl <- openxlsx::read.xlsx(file)
  tbl$COSMIC_ID <- paste('COSM:',tbl$COSMIC_ID,sep='')
  cells <- sort(unique(tbl$COSMIC_ID))
  drugs <- sort(unique(tbl$DRUG_NAME))
  sensitivityData <- matrix(NA,nrow=length(drugs),ncol=length(cells),dimnames = list(drugs,cells))
  sensitivityData[as.matrix(tbl[,c('DRUG_NAME','COSMIC_ID')])] <- tbl$LN_IC50
  
  cell_meta <- PharmacoGx::cellInfo(GDSCv1) %>% tibble::as_tibble()
  cell_meta$COSMIC.identifier <-paste('COSM:',cell_meta$COSMIC.identifier,sep='')
  cell_meta <- cell_meta[match(colnames(sensitivityData),cell_meta$COSMIC.identifier),]
  
  logi_liquid <- cell_meta$unique.tissueid.fromstudies=='haematopoietic_and_lymphoid_tissue'
  sensitivityList <- list(solid=sensitivityData[,!logi_liquid],liquid=sensitivityData[,logi_liquid])
  
  drug_meta <- PharmacoGx::drugInfo(GDSCv1) %>% tibble::as_tibble()
  
  return(list(sens=sensitivityList,cell_meta=cell_meta,drug_meta=drug_meta))    
}

resampled_cor <- function(x,y){
  thr <- round(length(x)*0.4)
  pair_indx <- !(is.na(x) | is.na(y))
  pairs <- sum(pair_indx)
  if(pairs>thr){
    n_sets <- choose(n=pairs,k=thr)
    if(n_sets>1000) r <- mean(replicate(1000,{idx<-sample(which(pair_indx),thr);return(cor(x[idx],y[idx]))}))
    else{
      cbn <- combn(which(pair_indx),thr)
      r <- mean(apply(cbn,2,function(idx)cor(x[idx],y[idx])))
    }
  }
  else r <- NA
  return(r)
}

compute_PCC <- function(expr,sens,cl){
  comm_ccl <- intersect(colnames(expr),colnames(sens))

  PCC_mtx <- parSapply(cl, rownames(expr),function(i){
    PCC <- vector(mode = 'numeric',length = nrow(sens))
    x<-expr[i,comm_ccl]
    for(j in 1:nrow(sens)){
      y<-sens[j,comm_ccl]
      PCC[j] <- resampled_cor(x,y)
    }
    names(PCC) <- rownames(sens)
    return(PCC)
  })
  
  not_na <- !is.na(PCC_mtx)
  PCC_mtx <- t(PCC_mtx[!rowSums(not_na)<500,])
  
  return(PCC_mtx)
}

#####

setwd('/mnt/d/data/')

ccle_expr <- import_CCLE_expression_data()
gdsc_expr <- import_GDSC_expression_data()
HB_expr <- import_HB_expression_data()
PPTC_expr <- import_PPTC_expression_data()
PPTC_expr2 <- import_PPTC_expression_data2()
TCGA_expr <- import_TCGA_expression_data()

CTRPv2_sens <- import_CTRPv2_sensitivity()
gdsc_sens <- import_GDSCv1_sensitivity()


cl <- makeCluster(100)
clusterExport(cl, c("compute_PCC","resampled_cor"))

set.seed(123)
PCC_ccle_CTRPv2 <- mapply(function(x,y)compute_PCC(x,y,cl),ccle_expr$expr,CTRPv2_sens$sens,SIMPLIFY = F)
PCC_gdsc <- mapply(function(x,y)compute_PCC(x,y,cl),gdsc_expr$expr,gdsc_sens$sens,SIMPLIFY = F)

stopCluster(cl)

