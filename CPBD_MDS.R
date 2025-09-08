#Cellphone DB for MDS vs. NBM

write_dgCMatrix <- function(mat,
                            filename,
                            col1_name = "gene",
                            chunk_size = 1000) {
  
  
  # Transpose so retrieval of "rows" is much faster
  mat <- Matrix::t(mat)
  
  # Row names
  row_names <- colnames(mat)
  
  # gene names are now columns
  col_names <- rownames(mat)
  
  n_row <- length(row_names)
  n_col <- length(col_names)
  
  n_chunks <- floor(n_row/chunk_size)
  
  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  names(chunk_df)[1] <- col1_name
  data.table::fwrite(chunk_df, file = filename, append = F, sep = '\t')
  
  # chunkation over chunksm
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Writing rows ",chunk_start," to ", chunk_end))
    chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
    chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
    data.table::fwrite(chunk_df, file = filename, append = T, sep = '\t')
  }
  
  # Remaining samples
  chunk_start <- (n_chunks*chunk_size + 1)
  chunk_end <- n_row
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  data.table::fwrite(chunk_df, file = filename, append = T, sep = '\t')
  
}

MDS_NBM_all_mt<-readRDS("~/single_cell_sequencing/MDS single cell//MDS_NBM_all_clustify_mutant_WT.RDS")

MDS_NBM_all_mt$type_for_LS<-as.character(MDS_NBM_all_mt$type_for_LS)
Idents(MDS_NBM_all_mt)<-"type_for_LS"  
DimPlot(MDS_NBM_all_mt, split.by = "Group", repel=T, label=T)
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="CD14+ Mono"]<-"CD14Mono"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="CD16+ Mono"]<-"CD16Mono"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="HSC/MPPs"]<-"HSCMPP"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="GMP/Promonocytes"]<-"Mprog"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="Late erythroid progenitors"]<-"Lprog"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="Plasma cells"]<-"Plasma"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="CD4 T"]<-"CD4T"
MDS_NBM_all_mt$type_for_LS[MDS_NBM_all_mt$type_for_LS=="CD8 T"]<-"CD8T"


MDS_NBM_all_mt_CPDB<-subset(MDS_NBM_all_mt, cells=sample(names(MDS_NBM_all_mt$Group),40000, replace = F))
MDS_NBM_all_mt_CPDB <- NormalizeData(MDS_NBM_all_mt_CPDB, normalization.method = "RC", scale.factor = 10000)
MDS_NBM_all_mt_CPDB@meta.data<-MDS_NBM_all_mt_CPDB@meta.data[c(1,6,25)]
MDS_NBM_all_mt_CPDB$type_for_LS<-paste0(MDS_NBM_all_mt_CPDB$type_for_LS,"_",MDS_NBM_all_mt_CPDB$Group)
Idents(MDS_NBM_all_mt_CPDB)<-"type_for_LS"



count_MDS_NBM_all<- MDS_NBM_all_mt_CPDB@assays$RNA@data
write_dgCMatrix(count_MDS_NBM_all, filename="count_MDS_NBM_all_counts.txt")
meta_data_MDS_NBM_all<- cbind(rownames(MDS_NBM_all_mt_CPDB@meta.data), MDS_NBM_all_mt_CPDB@meta.data[,"type_for_LS", drop=F])   #####  cluster is the user’s specific cluster column
write.table(meta_data_MDS_NBM_all, file="MDS_NBM_all_metadata.txt", sep="\t", quote=F, row.names=F)


#CellphondBD processing (Using python in terminal)
conda activate cpdb 

cellphonedb method statistical_analysis CellphoneDB/MDS_NBM_all_metadata.txt  CellphoneDB/count_MDS_NBM_all_counts.txt --counts-data=gene_name --output-format=txt --output-path=CellphoneDB/CPBD_MDS_NBM_all --project-name=CPBD_MDS_NBM_all --subsampling --subsampling-log true  
cellphonedb plot heatmap_plot CellphoneDB/AML_NBM_all_metadata.txt --output-path=CellphoneDB/CPBD_MDS_NBM_all --pvalues-path=CellphoneDB/CPBD_MDS_NBM/CPBD_MDS_NBM_all/pvalues.txt



##### Analysis 

#CellphondBD data analysis
x<-CPBD_all_analysis_no_integrin[[1]]
L_R_Recombination<-function(x){
  r_l<-x%>%filter(interacting_pair%in%receptor_ligand)
  l_r<-x%>%filter(!interacting_pair%in%receptor_ligand)
  pairs_to_reverse<-as.vector(sapply(receptor_ligand,function(x){
    x=strsplit(x,"_")
    x=paste0(x[[1]][2],"_",x[[1]][1])
    return(x)
  }))
  r_l<-r_l%>%mutate(interacting_pair=pairs_to_reverse)
  names(r_l)[2: ncol(r_l)]<-as.vector(
    sapply(names(r_l)[2: ncol(r_l)],function(x){
      x=strsplit(x,"\\.")
      x=paste0(x[[1]][2],"-",x[[1]][1])
      return(x)
    }))
  
  names(l_r)[2: ncol(l_r)]<-as.vector(
    sapply(names(l_r)[2: ncol(l_r)],function(x){
      x=strsplit(x,"\\.")
      x=paste0(x[[1]][1],"-",x[[1]][2])
      return(x)
    }))
  reversed_data<-rbind(l_r,r_l)
  return(reversed_data)
}

############# AML vs Ctrl BMSCs to other components
###all
CPBD_fils_path_all<-list.files(path="CellphoneDB/CPBD_MDS_NBM_all/CPBD_MDS_NBM_all",full.names=T)
CPBD_fils_name_all<-list.files(path="CellphoneDB/CPBD_MDS_NBM_all/CPBD_MDS_NBM_all",full.names=F)
CPBD_fils_all<-lapply(CPBD_fils_path_all,function(x){
  data.frame(read.table(x,header=T,sep="\t"))
})
names(CPBD_fils_all)<-CPBD_fils_name_all

CPBD_all_analysis<-CPBD_fils_all[c("means.txt","pvalues.txt")]


CPBD_all_analysis_mean<-CPBD_all_analysis[[1]]
CPBD_all_analysis_mean<-CPBD_all_analysis_mean[-c(1,3:10)]

CPBD_all_analysis_pval<-CPBD_all_analysis[[2]]
CPBD_all_analysis_pval<-CPBD_all_analysis_pval[-c(1,3:10)]

sig_all_index<-apply(CPBD_all_analysis_pval[-c(1,2)],1,function(x){
  any(as.numeric(x)<0.05)
})

CPBD_all_analysis_pval<-CPBD_all_analysis_pval[sig_all_index,]
CPBD_all_analysis_pval_no_integrin<-CPBD_all_analysis_pval%>%filter(is_integrin=="False")%>%dplyr::select(!is_integrin)
CPBD_all_analysis_mean_no_integrin<-CPBD_all_analysis_mean%>%
  filter(interacting_pair%in%CPBD_all_analysis_pval_no_integrin$interacting_pair)%>%dplyr::select(!is_integrin)

CPBD_all_analysis_no_integrin<-list(CPBD_all_analysis_pval_no_integrin,CPBD_all_analysis_mean_no_integrin)


CPBD_all_analysis_pairs<-CPBD_all_analysis_pval_no_integrin$interacting_pair

receptor_ligand<-CPBD_all_analysis_pval_no_integrin$interacting_pair[c(1:3,8:15,17:30,32:51,54:72,89,99:110,
                                                                       145:148,158:159,169,176,180,197:207,208:213,214,215,217:222,225:233)]

L_R_Recombination<-function(x){
  r_l<-x%>%filter(interacting_pair%in%receptor_ligand)
  l_r<-x%>%filter(!interacting_pair%in%receptor_ligand)
  pairs_to_reverse<-as.vector(sapply(receptor_ligand,function(x){
    x=strsplit(x,"_")
    x=paste0(x[[1]][2],"_",x[[1]][1])
    return(x)
  }))
  r_l<-r_l%>%mutate(interacting_pair=pairs_to_reverse)
  names(r_l)[2: ncol(r_l)]<-as.vector(
    sapply(names(r_l)[2: ncol(r_l)],function(x){
      x=strsplit(x,"\\.")
      x=paste0(x[[1]][2],"-",x[[1]][1])
      return(x)
    }))
  
  names(l_r)[2: ncol(l_r)]<-as.vector(
    sapply(names(l_r)[2: ncol(l_r)],function(x){
      x=strsplit(x,"\\.")
      x=paste0(x[[1]][1],"-",x[[1]][2])
      return(x)
    }))
  reversed_data<-rbind(l_r,r_l)
  return(reversed_data)
}

CPBD_all_analysis_no_integrin<-lapply(CPBD_all_analysis_no_integrin,L_R_Recombination)


# BMSC to ALL

source<-"BMSC"
target<-"NK"

interaction_with<-function(source, target){
  CPBD_MSC_to_all<-lapply(CPBD_all_analysis_no_integrin,function(x){
    x<-x%>%
      select(interacting_pair,starts_with(source))%>%select(interacting_pair,contains(target))
    x_NBM<-x%>%select(!contains("MDS"))
    x_MDS<-x%>%select(!contains("NBM"))
    x<-x_NBM%>%left_join(x_MDS)
  })
  
  CPBD_MSC_to_all<-lapply(CPBD_MSC_to_all,function(x){
    index<-apply(CPBD_MSC_to_all[[1]][-1],1,function(y){
      any(as.numeric(y)<0.05)})
    x<-x[index,]
    return(x)
  })
  
  CPBD_MSC_to_all_pval_no_integrin<-CPBD_MSC_to_all[[1]]%>%melt(id.vars="interacting_pair",variable.name="Cell_type",value.name="pval")
  CPBD_MSC_to_all_mean_no_integrin<-CPBD_MSC_to_all[[2]]%>%melt(id.vars="interacting_pair",variable.name="Cell_type",value.name="mean")
  CPBD_all_analysis_final<-CPBD_MSC_to_all_pval_no_integrin%>%left_join(CPBD_MSC_to_all_mean_no_integrin)
  CPBD_all_analysis_final$Cell_type<-gsub(pattern = "\\.", replace="-",x= CPBD_all_analysis_final$Cell_type)
  CPBD_all_analysis_final["target"]<- str_split(CPBD_all_analysis_final$Cell_type, "-", simplify = T)[,2]
  CPBD_all_analysis_final["target"]<- str_split(CPBD_all_analysis_final$target, "_", simplify = T)[,1]
  CPBD_all_analysis_final<-CPBD_all_analysis_final[CPBD_all_analysis_final$target %in% target,]
  CPBD_all_analysis_final<-CPBD_all_analysis_final%>%mutate(Group=case_when(str_detect(Cell_type,"NBM")~"NBM",TRUE~"MDS"))
  CPBD_all_analysis_final$target<-factor(CPBD_all_analysis_final$target,levels=c(target))
  CPBD_all_analysis_final<-CPBD_all_analysis_final%>%group_by(interacting_pair,target=target)%>%
    summarise(Group=Group,Cell_type=Cell_type,source=source,mean=mean,pval=pval,relative_mean=mean/sum(mean))
  #CPBD_all_analysis_final$Type<-factor(CPBD_all_analysis_final$Type,levels=c("BMSC","EC","HSC","MPP","MEP","GMP","CLP","Mk","DC","CD14+ Mono","CD16+ Mono","CD4 T","CD8 T","NK","B","Plasma cells"))         
  
  CPBD_all_analysis_final%>%
    mutate(Group=factor(Group,levels=c("NBM","MDS")))%>%mutate(relative_mean=case_when(relative_mean>0.75~0.75, relative_mean<0.25~0.25, TRUE~relative_mean))%>%
    ggplot(aes(x=Group,y=interacting_pair,size= (-log10(pval+0.0001)),col=relative_mean))+geom_point()+
    scale_size_continuous(range = c(0.1, 3))+facet_wrap(~target,scales="free_x",ncol=16)+theme_classic()+
    theme(axis.text.x=element_text(angle=90,vjust=1,size=10),axis.text.y=element_text(size=7),axis.title=element_blank(),strip.text = element_text(size=10),strip.background = element_blank())+
    scale_colour_viridis_c(limits = c(0.25, 0.75))+labs(size="-log10(pval)")+labs(title=paste0("source: ",source))
}

interaction_with(source="BMSC",target=c("HSCMPP","LMPPs","Mprog","MEPs","CLPs"))

interaction_with(source="BMSC",target=c("NK"))

interaction_with(target = "CD14Mono",source=c("CD8T"))
interaction_with(target = "HSCMPP",source=c("NK"))
interaction_with(target = "HSCMPP",source=c("CD8T"))

interaction_with(target = "NK",source=c("Mprog"))
interaction_with(target = "NK",source=c("CD14Mono"))

