get_gold_stds <- function(sens_mtx){
  zscore <- apply(sens_mtx,1,mosaic::zscore,na.rm=T)
  perc5 <- zscore<=qnorm(0.05)
  idx <- which(perc5,arr.ind = T)
  df <- data.frame(cellid=rownames(zscore)[idx[,1]],drugid=colnames(zscore)[idx[,2]])
  return(list(zscore=zscore,gld_std=df))
}

ppvpercents <- function(gld_std,test,CCL_sel_ids=NULL,norm_ppv=T){
  CCL_ids <- intersect(rownames(gld_std$zscore),unique(test$sample.id))
  if(!is.null(CCL_sel_ids)) CCL_ids <- intersect(CCL_ids,CCL_sel_ids)
  test <- test[test$sample.id%in%CCL_ids,]
  test <- test[order(test$`G-value`),]
  ind <- prodlim::row.match(gld_std$gld_std,test[,1:2],nomatch = NA)
  ind <- sort(ind[!is.na(ind)])
  if(!is.null(ind)){
    nconn <- nrow(test)
    ppv <- vector()
    iters <- seq(1,100,1)
    for (i in 1:length(iters)){
      perc <- round(nconn/100 * iters[i])
      if(perc ==0) perc <- 1
      ppv[[i]] <- sum(ind<perc)/perc
      if(norm_ppv){
        ppv[[i]] <- ppv[[i]]/(length(ind)/nconn)
      }
    }
    return(ppv)
  }
  else return(NULL)
}

rename_cell <- function(pred,meta){
  pred$sample.id <- meta$Characteristics.cell.line.[match(pred$sample.id,meta$COSMIC_ID)]
  return(pred)
}

#####

setwd('/mnt/d/data')

CTRPv2_sens <- readRDS('rds/CTRPv2_sens.rds')
gdsc_cell_meta <- readRDS('rds/gdsc_expr.rds')$cell_meta
gdsc_cell_meta$Characteristics.cell.line. <- toupper(gsub("[^[:alnum:] ]", "", gdsc_cell_meta$Characteristics.cell.line.))
ccle_cell_meta <- readRDS('rds/ccle_expr.rds')$cell_meta

gld_stds <- lapply(CTRPv2_sens$sens,get_gold_stds)

fls <- list.files('predicted_efficacy',full.names = T)
sel_fls <- grep('broad',fls,value = T)[1:8]

eff <- lapply(sel_fls,readRDS)
names(eff) <- gsub("predicted_efficacy/efficacy_|.rds","",sel_fls)

idx_gdsc <- grep('gdsc',names(eff),value=T)
eff[idx_gdsc] <- lapply(eff[idx_gdsc],function(x)rename_cell(x,gdsc_cell_meta))

gld_type <- stringi::stri_extract_first(sel_fls,regex = 'solid|liquid')

res <- mapply(function(x,i)ppvpercents(gld_stds[[i]],x),eff,gld_type,SIMPLIFY = F)

####################################### VALIDATION CCLs

#General CCLE

pos_control <- list(ccle_model=list(liquid=res$ccle_liquid_broad_liquid,solid=res$ccle_solid_broad_solid))
df0 <- reshape2::melt(pos_control)

neg_control <- list(ccle_model=list(liquid = res$ccle_liquid_broad_solid,solid=res$ccle_solid_broad_liquid))
df1 <- reshape2::melt(neg_control)

#General GDSC

pos_control <- list(gdsc_model=list(liquid = res$gdsc_liquid_broad_liquid,solid=res$gdsc_solid_broad_solid))
df2 <- reshape2::melt(pos_control)

neg_control <- list(gdsc_model=list(liquid = res$gdsc_liquid_broad_solid,solid=res$gdsc_solid_broad_liquid))
df3 <- reshape2::melt(neg_control)

#####

l_map <- data.frame(extended=c("Ewing_sarcoma","hepatocellular_carcinoma","medulloblastoma","neuroblastoma"),abbreviation=c('EW','LIHC','MB','NB'))
for(i in 1:nrow(l_map))ccle_cell_meta$lineage_subtype <- gsub(l_map[i,1],l_map[i,2],ccle_cell_meta$lineage_subtype)

paediatric_tumor <- c('AML','ALL','EW','LIHC','MB','NB')

ccle_paediatric <- ccle_cell_meta[ccle_cell_meta$lineage_subtype%in%paediatric_tumor,]
ccle_paediatric <- split(ccle_paediatric$stripped_cell_line_name,ccle_paediatric$lineage_subtype)

gdsc_paediatric <- gdsc_cell_meta[gdsc_cell_meta$ct%in%paediatric_tumor,]
gdsc_paediatric <- split(gdsc_paediatric$Characteristics.cell.line.,gdsc_paediatric$ct)

ccle_gdsc_paediatric <- mapply(c,gdsc_paediatric,ccle_paediatric)

res <- mapply(function(x,i)lapply(ccle_gdsc_paediatric,function(y)ppvpercents(gld_stds[[i]],x,CCL_sel_ids = y)),eff,gld_type,SIMPLIFY = F)

#Pediatric Tumor CCLE

pt_ccle <- list(ccle_model=list(liquid = res$ccle_liquid_broad_liquid,solid=res$ccle_solid_broad_solid))
df4 <- reshape2::melt(pt_ccle)
df4 <- df4[!is.na(df4$value),]

#Pediatric Tumor GDSC

pt_gdsc <- list(gdsc_model=list(liquid = res$gdsc_liquid_broad_liquid,solid=res$gdsc_solid_broad_solid))
df5 <- reshape2::melt(pt_gdsc)
df5 <- df5[!is.na(df5$value),]

####################################### PLOT GENERATION

dfs <- list(df0,df1,df2,df3)

p <- lapply(dfs,function(df){
  ggplot(data=df,aes(x=rep(1:100,nrow(df)/100), y=value,colour = L2)) +
    geom_hline(yintercept = 1, colour='black') + xlab('% CCL-Drug pairs ranked by G-score') + 
    ylab('Normalized PPV') + geom_line(size=1) + theme_minimal() +
    scale_color_manual(values=c('#b2182b','azure4')) +
    theme(legend.title = element_blank(),panel.border = element_blank(), panel.grid.minor = element_blank(),
          legend.position = 'none', text = element_text(size=rel(3)),strip.text =  element_text(face='bold',size=rel(3.5)),
          strip.background = element_blank(),strip.text.x = element_blank(),
          legend.text = element_text(size=rel(2.5)))
})

p[[1]] <- p[[1]] + scale_y_continuous(breaks=c(4,8,12,16),limits=c(0.5,16))
p[[2]] <- p[[2]] + scale_y_continuous(breaks=c(4,8,12,16),limits=c(0.5,16))
p[[3]] <- p[[3]] + scale_y_continuous(breaks=c(3,6,9),limits=c(0.5,11))
p[[4]] <- p[[4]] + scale_y_continuous(breaks=c(3,6,9),limits=c(0.5,11))
 


dfs_pt <- list(df4,df5)

p_pt <- lapply(dfs_pt,function(df){
  ggplot(data=df,aes(x=rep(1:100,nrow(df)/100), y=value,colour = L2)) +
    geom_hline(yintercept = 1, colour='black') + xlab('% CCL-Drug pairs ranked by G-score') + 
    ylab('Normalized PPV') + geom_line(size=1) + facet_wrap(~L3,scales = 'free') + theme_minimal() +
    scale_color_manual(values=c('#b2182b','azure4')) +
    theme(legend.title = element_blank(),panel.border = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size=rel(3)),strip.text =  element_text(face='bold',size=rel(3.5)), legend.position = 'none',
          legend.text = element_text(size=rel(2.5)))
})
  


pdf('Figures/Figure2.pdf',width = 5,height = 2.5,fonts = 'Arial')
easyGgplot2::ggplot2.multiplot(p[[1]],p_pt[[1]])
dev.off()


pdf('Figures/Supplementary_S1_part1.pdf',width = 2.5,height = 2.5,fonts = 'Arial')
p[[2]]
dev.off()


pdf('Figures/Supplementary_S1_part2.pdf',width = 7.5,height = 2.5,fonts = 'Arial')
easyGgplot2::ggplot2.multiplot(p[[3]],p_pt[[2]],p[[4]],cols = 3)
dev.off()


####################################### VENN DIAGRAM CCLs

library('ggvenn')

ccle_expr <- readRDS('rds/ccle_expr.rds')
gdsc_expr <- readRDS('rds/gdsc_expr.rds')
CTRPv2_sens <- readRDS('rds/CTRPv2_sens.rds')

gdsc_expr$cell_meta$Characteristics.cell.line. <- toupper(gsub("[^[:alnum:] ]", "", gdsc_expr$cell_meta$Characteristics.cell.line.))

venn_ccle_ctrpv2 <- mapply(function(i,j)list(ccle=colnames(i),ctrpv2=colnames(j)),ccle_expr$expr,CTRPv2_sens$sens,SIMPLIFY=F)

venn_gdsc_ctrpv2 <- mapply(function(i,j)list(gdsc=gdsc_expr$cell_meta$Characteristics.cell.line.[match(colnames(i),gdsc_expr$cell_meta$COSMIC_ID)],
                                             ctrpv2=colnames(j)),gdsc_expr$expr,CTRPv2_sens$sens,SIMPLIFY=F)

venn_lists <- c(venn_ccle_ctrpv2,venn_gdsc_ctrpv2)
v <- lapply(venn_lists,function(x)ggvenn(x,fill_color = c("black","white"),show_percentage = F))

pdf('Figures/Supplementary_S1_part1_venn_solid.pdf',width = 2.5,height = 2.5)
easyGgplot2::ggplot2.multiplot(v[[1]],v[[3]],cols = 1)
dev.off()

pdf('Figures/Supplementary_S1_part1_venn_liquid.pdf',width = 2.5,height = 2.5)
easyGgplot2::ggplot2.multiplot(v[[2]],v[[4]],cols = 1)
dev.off()

