ppv_PDXs <- function(reference,test,norm_ppv=T,percentiles=T){
  reference$sample.id <- gsub('\\.[0-9]+','',rownames(reference))
  reference <- reference[,c(3,1,2)]
  ind <- prodlim::row.match(reference[,1:2],test[,1:2])
  reference <- reference[!is.na(ind),]
  test <- test[na.omit(ind),]
  logi <- reference[,3]%in%c('CR','MCR','PR')
  logi <- logi[order(test$`G-value`)]
  
  ind <- which(logi)
  if(percentiles){
    nconn <- nrow(test)
    ppv <- vector()
    iters <- seq(1,100,1)
    rand <- length(ind)/nconn
    for (i in 1:length(iters)){
      perc <- round(nconn/100 * iters[i])
      if(perc ==0) perc <- 1
      ppv[i] <- sum(ind<perc)/perc
    }
  }
  else {
    ppv <- cumsum(logi)/(1:length(logi))
    rand <- sum(logi)/length(logi)
  }
  if(norm_ppv){
    ppv <- ppv/rand
  }
  return(ppv)
}

#####

setwd('/mnt/d/data')

preclinical_fls <- list.files('/mnt/d/PCAT/Preclinical')
dir_path <- '/mnt/d/PCAT/Preclinical/'
preclinical_data <- lapply(preclinical_fls,function(x)read.csv(paste(dir_path,x,sep = ''),stringsAsFactors = F))
names(preclinical_data) <- gsub('.csv','',preclinical_fls)

PDX_drug <- openxlsx::read.xlsx('PDX_drug.xlsx')
PDX_drug <- PDX_drug[!is.na(PDX_drug$CTRPv2),]

preclinical_data <- lapply(preclinical_data,function(x)merge(PDX_drug,x,by.x='PCAT',by.y='Drug_Name'))
preclinical_data <- preclinical_data[sapply(preclinical_data,nrow)>0]

PCAT_ids <- read.csv('/mnt/d/data/Expression Data/PCAT_meta.csv')
PCAT_ids <- PCAT_ids[match(names(preclinical_data),PCAT_ids$sample_id),]

hist_det. <- split(PCAT_ids$sample_id,PCAT_ids$Histology_Detailed)
preclinical_data_TT <- lapply(hist_det.,function(i)do.call(rbind,preclinical_data[i]))

preclinical_data_TT <- preclinical_data_TT[sapply(preclinical_data_TT,function(x)sum(x$Response_Level_annotation%in%c('CR','PR','MCR')))>=6]

eff_PPTC <- list(solid=readRDS('predicted_efficacy/efficacy_PPTC_solid_broad_solid.rds'),
                           liquid=readRDS('predicted_efficacy/efficacy_PPTC_liquid_broad_liquid.rds'))

#Histology Detailed 

ppv <- vector(mode = 'list')

for(i in names(preclinical_data_TT)){
  reference <- preclinical_data_TT[[i]][,c('CTRPv2','Response_Level_annotation')]
  if(grepl('ALL',i)) test <- eff_PPTC$liquid
  else test <- eff_PPTC$solid
  ppv[[i]] <- ppv_PDXs(reference,test,percentiles = F)
}

names(ppv) <- c('BCP-ALL','EW','MLL-ALL','NB','T-ALL','WT')

df0 <- reshape2::melt(ppv)
df0$L1 <- factor(df0$L1,levels=c('BCP-ALL','MLL-ALL','T-ALL','EW','NB','WT'))
df0$L2 <- grepl('ALL',df0$L1)
df0$value.x <- unlist(sapply(unique(df0$L1),function(x)1:sum(df0$L1==x)))

#####

preclinical_data <- split(preclinical_data,grepl('ALL',names(preclinical_data)))
preclinical_data <- lapply(preclinical_data,function(x)(do.call(rbind,x)))

#####

#General

ppv <- vector(mode = 'list')

for(i in names(preclinical_data)){
  reference <- preclinical_data[[i]][,c('CTRPv2','Response_Level_annotation')]
  if(i=='TRUE') test <- eff_PPTC$liquid
  else test <- eff_PPTC$solid
  ppv[[i]] <- ppv_PDXs(reference,test)
}

names(ppv) <- c('solid','liquid')

df1 <- reshape2::melt(ppv)

#####

eff_PPTC <- list(solid=readRDS('predicted_efficacy/efficacy_PPTC_solid_broad_liquid.rds'),
                 liquid=readRDS('predicted_efficacy/efficacy_PPTC_liquid_broad_solid.rds'))

##### 

#Negative Control

ppv <- vector(mode = 'list')

for(i in names(preclinical_data)){
  reference <- preclinical_data[[i]][,c('CTRPv2','Response_Level_annotation')]
  if(i=='TRUE') test <- eff_PPTC$liquid
  else test <- eff_PPTC$solid
  ppv[[i]] <- ppv_PDXs(reference,test)
}

names(ppv) <- c('solid','liquid')

df2 <- reshape2::melt(ppv)

####################################### PLOT GENERATION

dfs <- list(df1,df2) 

p <- lapply(dfs,function(df){
  ggplot(df,aes(x=rep(1:100,nrow(df)/100),y=value,colour = L1)) + geom_line(size=1) + geom_hline(aes(yintercept = 1)) +
    labs(x='% PDX-Drug pairs ranked by G-score', y='Normalized PPV') + theme_minimal() +
    scale_color_manual(values = c('#b2182b','azure4')) +
    theme(legend.title = element_blank(),panel.border = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size=rel(3)),strip.text =  element_text(face='bold',size=rel(3.5)), legend.position = 'none',
          legend.text = element_text(size=rel(2.5)))
})


p[[1]] <- p[[1]] + scale_y_continuous(breaks=c(0:4),limits=c(0,4.25)) 
p[[2]] <- p[[2]] + scale_y_continuous(breaks=c(0:4),limits=c(0,4.25)) 



p_hist_det <- ggplot(df0,aes(x=value.x,y=value,colour = L2)) + geom_line(size=1) + geom_hline(aes(yintercept = 1)) +
  facet_wrap(~ L1, scales = 'free')+ labs(x='PDX-Drug pairs ranked by G-score', y='Normalized PPV') + theme_minimal() +
  scale_color_manual(values = c('azure4','#b2182b')) + 
  theme(legend.title = element_blank(),panel.border = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=rel(3)),strip.text =  element_text(face='bold',size=rel(3.5)), legend.position = 'none',
        legend.text = element_text(size=rel(2.5))) 



pdf('Figures/Figure3.pdf',width = 5,height = 2.5,fonts = 'Arial')
easyGgplot2::ggplot2.multiplot(p[[1]],p_hist_det)
dev.off()


pdf('Figures/Supplementary_S1_part3.pdf',width = 2.5,height = 2.5,fonts = 'Arial')
p[[2]]
dev.off()

####################################### VENN DIAGRAM DRUGs

venn_drug <- list(list(solid_model=1:445,pcat=c(1:23,446:(446+47))),list(liquid_model=1:414,pcat=c(1:21,415:(415+34))))

v <- lapply(venn_drug,function(x)ggvenn(x,fill_color = c("black","white"),show_percentage = F))

pdf('Figures/Supplementary_S1_part3_venn_drugs.pdf',width = 2.5,height = 2.5)
easyGgplot2::ggplot2.multiplot(v[[1]],v[[2]],cols = 1)
dev.off()

