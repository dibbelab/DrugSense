get_gene_sets<-function(prof_mtx,reference,limit){
  idx <- rownames(prof_mtx)%in%rownames(reference)
  prof_mtx <- prof_mtx[idx,]
  
  signature <- vector(mode = 'list', length = ncol(prof_mtx))
  names(signature) <- colnames(prof_mtx)
  
  for(i in colnames(prof_mtx)){
    signature[[i]] <- vector(mode = 'list', length = 2)
    names(signature[[i]]) <- c('up_genes','down_genes')
    
    expr_val <- prof_mtx[!is.na(prof_mtx[,i]),i]
    sorted_genes <- names(sort(expr_val,decreasing=T))
    
    signature[[i]]$up_genes <- sorted_genes[1:limit]
    signature[[i]]$down_genes <- rev(sorted_genes)[1:limit]
  }
  
  signature <- unlist(signature,recursive = F)
  return(signature)
}

calcGseaStatCol <- function(gene_sets,single_gene_list){
  res<-vector(mode = 'numeric',length = length(gene_sets))
  for(i in seq_along(gene_sets)){
    idx <- fastmatch::fmatch(gene_sets[[i]],names(single_gene_list))
    res[i] <- fgsea::calcGseaStat(stats = single_gene_list,selectedStats = idx,gseaParam = 0)
  }
  return(res)
}

compute_gsea_fast <- function(gene_sets,multiple_gene_lists){
  res <- matrix(0,nrow=length(gene_sets),ncol=ncol(multiple_gene_lists),dimnames=list(names(gene_sets),colnames(multiple_gene_lists)))
  for(i in colnames(multiple_gene_lists)){
    sorted_list <- sort(multiple_gene_lists[,i],decreasing=T)
    res[,i] <- calcGseaStatCol(gene_sets,sorted_list)
  }
  return(res)
}

get_single_strat_pval <- function(stat,stat_rand,direction='down',reverse=FALSE){
  stat_rand <- stat_rand[grepl(paste(direction,'_genes',sep=''),rownames(stat_rand)),]
  rownames(stat_rand) <- gsub(paste('.',direction,'_genes',sep=''),'',rownames(stat_rand))
  
  stat <- stat[grepl(paste(direction,'_genes',sep=''),rownames(stat)),]
  rownames(stat) <- gsub(paste('.',direction,'_genes',sep=''),'',rownames(stat))
  
  if(reverse){
    stat_rand <- t(stat_rand)
    stat <- t(stat)
  }
  
  p_emp <- stat
  for(i in seq_along(colnames(stat))){
    for(j in seq_along(stat[,i])){
      if (direction=='down') p_emp[j,i] <- (sum(stat_rand[,i]>stat[j,i])+1)/length(stat_rand[,i])
      else p_emp[j,i] <-(sum(stat_rand[,i]<stat[j,i])+1)/length(stat_rand[,i])
    }
  }
  
  p_emp <- reshape2::melt(p_emp)
  p_emp <- p_emp[order(p_emp$value),]
  p_emp[,1:2] <- lapply(p_emp[,1:2], as.character)
  return(p_emp)
}

generate_random <- function(exp_mtx){
  time <- ceiling(10000/ncol(exp_mtx))
  rand_mtx<-do.call(cbind,replicate(time,exp_mtx,simplify=F))
  
  dms <- prod(dim(rand_mtx))
  thr <- round(dms*0.25)
  
  set.seed(123)
  idx <- sample(1:dms,thr)
  rand_mtx[idx] <- sample(rand_mtx[idx])
  colnames(rand_mtx) <- paste('rand',1:ncol(rand_mtx),sep='') 
  return(rand_mtx)
}

get_all_strat_pval <- function(exp_mtx,rand_mtx,model,limit=250){
  
  ### Get signatures of up and down-regulated genes from random expression profiles
  
  gene_sets <- get_gene_sets(rand_mtx,model,limit)
  mrk_sets <- get_gene_sets(model,rand_mtx,limit)
  
  ### Get Enrichment Scores (ES)
  
  gsea_rand2drug <- compute_gsea_fast(gene_sets,model)
  gsea_drug2rand <- compute_gsea_fast(mrk_sets,rand_mtx)
  
  ### Get signatures of up and down-regulated genes from expression profiles of CCLs
  
  gene_sets <- get_gene_sets(exp_mtx,model,limit)
  
  ### Get Enrichment Scores (ES)
  
  gsea_data2drug <- compute_gsea_fast(gene_sets,model)
  gsea_drug2data <- compute_gsea_fast(mrk_sets,exp_mtx)
  
  ### 
  
  a <- get_single_strat_pval(gsea_data2drug,gsea_rand2drug,direction='up')
  b <- get_single_strat_pval(gsea_data2drug,gsea_rand2drug)
  c <- get_single_strat_pval(gsea_drug2data,gsea_drug2rand,reverse=T)
  d <- get_single_strat_pval(gsea_drug2data,gsea_drug2rand,direction='up',reverse=T)
  
  ###
  
  all_strat_pval <- list(a=a,b=b,c=c,d=d)
  for(i in letters[1:4]) names(all_strat_pval[[i]])[3] <- paste(names(all_strat_pval[[i]])[3],i,sep='_')
  
  return(all_strat_pval)
}

geo.mean<-function (x, na.rm = TRUE) {
  return(exp(mean(log(x), na.rm = TRUE)))
}

get_drug_efficacy<-function(all_strat_pval){
  df<-Reduce(function(x,y)merge(x,y,by=c('Var1','Var2')),all_strat_pval)
  df<-data.frame(df[,1:2],G=apply(df[,3:6],1,geo.mean))
  colnames(df)<-c('sample.id','drug.id','G-value')
  return(df)
}


#####

PCC_ccle_CTRPv2 <- readRDS('rds/PCC_ccle_CTRPv2_filtered.rds')
PCC_gdsc <- readRDS('rds/PCC_gdsc_filtered.rds')
PCC <- unlist(list(PCC_ccle_CTRPv2,PCC_gdsc),recursive = F)
names(PCC) <- paste(rep(c('broad','sanger'),each=2),names(PCC),sep = '_')
rm(PCC_ccle_CTRPv2,PCC_gdsc)

ccle_expr <- readRDS('rds/ccle_expr.rds')
gdsc_expr <- readRDS('rds/gdsc_expr.rds')
PPTC_expr <- readRDS('rds/PPTC_expr.rds')
HB_expr <- readRDS('rds/HB_expr.rds')
expr <- unlist(list(ccle_expr$expr,gdsc_expr$expr,PPTC_expr$expr,HB=list(HB_expr)),recursive = F)
names(expr)[1:6] <- paste(rep(c('ccle','gdsc','PPTC'),each=2),names(expr)[1:6],sep = '_')
rm(ccle_expr,gdsc_expr,PPTC_expr,HB_expr)

#####

grd <- as.matrix(expand.grid(expression=names(expr)[1:6], sensitivity=names(PCC)))
grd<- rbind(grd, c(names(expr)[7],'broad_solid'))


cl <- parallel::makeCluster(25)
parallel::clusterExport(cl,ls())

parallel::parSapply(cl,1:nrow(grd),function(i){
  m <- grd[i,'expression']
  n <-  grd[i,'sensitivity']
  rand_mtx <- generate_random(expr[[m]])
  all_strat_pval <- get_all_strat_pval(expr[[m]],rand_mtx,PCC[[n]])
  drug_efficacy <- get_drug_efficacy(all_strat_pval)
  saveRDS(drug_efficacy,file=paste('predicted_efficacy/efficacy_',m,'_',n,'.rds',sep=''))
})

parallel::stopCluster(cl)
