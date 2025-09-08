suppressMessages(library('parallel'))
options(future.rng.onMisuse="ignore")
options(future.globals.maxSize = 2e4 * 1024^2)
options(future.fork.enable = TRUE)

#plan(strategy = 'multicore', workers = 16)
plan(strategy = 'sequential')


densityplot<-function(dataset=dataset,assay="RNA",features=features,split.by=NULL){
  DefaultAssay(dataset)<-"RNA"
  density_plot<-plot_density(dataset, features)
  density_plot$data<-cbind(density_plot$data,dataset@meta.data[split.by])
  if(is.null(split.by)){density_plot}else{density_plot+facet_wrap(.~Group)}
}

library("Nebulosa")
densityplot<-function(dataset=dataset,assay="RNA",features=features,split.by=NULL){
  DefaultAssay(dataset)<-"RNA"
  density_plot<-plot_density(dataset, features)
  density_plot$data<-cbind(density_plot$data,dataset@meta.data[split.by])
  density_plot
}




data=Niche_NBM_MDS
contrasts='condition,MDS,NBM'
condition_col="Group"
extr_col=NULL
cluster_col="type"
design='~ condition'
sample_col="Sample"
do_prefiltering=T

run_pseudobulkDE_new <- function(data,
                                 contrasts,
                                 design,
                                 sample_col,
                                 condition_col,
                                 cluster_col,
                                 extr_col = NULL, # + Extra
                                 batch_correction=F,
                                 do_prefiltering = TRUE) {
  
  
  # Pull the necessary columns from the metadata table
  meta_data<-data@meta.data
  anno_raw <- meta_data[, c(sample_col, condition_col, extr_col)]
  
  # Pull names from the additional column(s)
  extr_name <- colnames(anno_raw)[3:ncol(anno_raw)]
  # Replace forbidden chrs and add a chr in columns if values start with int
  anno_raw <- apply(anno_raw, 2, function(x) {
    x <- ifelse(grepl('^[0-9]', x), paste0('x', x), x)
    gsub('-|_', '.', x)
  })
  
  
  # Create cell-level annoation dataframe with columns as factors
  anno_cell <- data.frame(
    samples          = anno_raw[, sample_col],
    condition        = anno_raw[, condition_col],
    Extra            = anno_raw[, extr_col],
    row.names        = rownames(meta_data),
    stringsAsFactors = TRUE
  )
  # Add column name(s) to the additional column(s)
  
  # Add factor column that is the conjunction of samples and conditions
  anno_cell$sc <- factor(paste0(anno_cell$samples, '.', anno_cell$condition))
  # split the contrast into a strings
  contrast <- strsplit(contrasts, ',')[[1]]
  
  # Turn the design parameter into a formula object
  design <- formula(design)
  
  # Create groups that will be used to aggregate the expr. in the sparse matrix
  
  agg_g <- S4Vectors::DataFrame(anno_cell[, 'sc'])
  
  # Aggregate the RNA counts by summation; making a sample-level count-matrix
  data[["Sample1"]]<-paste(data[[sample_col]][,1],data[[condition_col]][,1])
  
  agg_m_seurat <- AggregateExpression(data, assays = "RNA", return.seurat = T, group.by = "Sample1")
  
  agg_m<-agg_m_seurat@assays$RNA$counts
  colnames(agg_m) <- as.character(sapply(colnames(agg_m), function(x) {
    x <- ifelse(grepl('^[0-9]', x), paste0('x', x), x)
    x<-gsub('-|_', '.', x)
    gsub(' ', '.', x)
  }))
  
  # Simply transpose agg.m
  
  # Create the intial dataframe used by DESeq2 for determining the conditions
  anno_con <- unique(anno_cell[,-c(1)])
  
  
  # Subset the condition annotation based on the condition being present in the contrast of interest
  anno_con <- anno_con[anno_con$condition %in% contrast, ]
  
  
  # Create the rownames of the condition annotation by pasting samples and condition
  rownames(anno_con) <- anno_con$sc
  
  # change the name of the agg.m object to adhere to the rest of the script
  counts_obj <- agg_m
  
  
  # Subset the full gene-count matrix based on the samples in the condition annotation
  count_name <- colnames(counts_obj)
  anno_name <- rownames(anno_con)
  counts_obj <- counts_obj[, which(count_name %in% anno_name)]
  
  # Reorder the condition annotation based on the column names in counts_obj
  anno_con <- anno_con[match(colnames(counts_obj), anno_name), ]
  
  # Batch effect correction using RUSeq
  dds <- DESeq2::DESeqDataSetFromMatrix(counts_obj, colData = anno_con,
                                        design = design)
  
  
  # Remove outliers, i.e. genes only expressed in a few samples
  if (do_prefiltering) {
    
    message('Info: perform pre-filtering')
    
    con_id <- anno_con[, contrast[1]] == contrast[2]
    
    fragments <- DESeq2::fpm(dds, robust = TRUE)
    
    keep_g1 <- rowSums(fragments[, con_id] >= 2) >= round(sum(con_id) / 2)
    keep_g2 <- rowSums(fragments[, !con_id] >= 2) >= round(sum(!con_id) / 2)
    
    keep_tot <- keep_g1 | keep_g2
    
    dds <- DESeq2::DESeq(dds[keep_tot, ], parallel = FALSE)
    
  } else {
    dds <- DESeq2::DESeq(dds, parallel = FALSE, minReplicatesForReplace = Inf)
  }
  
  # Compute normalized counts
  norm_counts <- as.matrix(DESeq2::counts(dds, normalized = TRUE))
  
  # Generate result table and compute shrunken log fold change
  res <- DESeq2::results(dds, contrast = contrast)
  res_lfc <- DESeq2::lfcShrink(dds, contrast = contrast, type = 'normal',
                               quiet = TRUE)
  res_lfc_ashr <- DESeq2::lfcShrink(dds, contrast = contrast, type = 'ashr',
                                    quiet = TRUE)
  
  # Perpare all the output in a list
  output_l <- list(
    dds          = dds,
    norm_counts  = norm_counts,
    res          = res,
    res_lfc      = res_lfc,
    res_lfc_ashr = res_lfc_ashr
  )
  # Return the DE analyis output list
  Output<-data.frame(cbind(norm_counts,res_lfc_ashr))
  return(Output)
}


#Data process for each samples
###############

#NBM 
Niche_NBM<-readRDS("~NBM_niche_merge.RDS")

#MDS

MDS337_niche_immune <- CreateSeuratObject(counts = Read10X(data.dir = "~/share_epronk/MDS_337_Niche_Immune/outs/filtered_feature_bc_matrix/") , min.cells = 3 , min.features = 200 , project = 'MDS337')
MDS345_niche_immune <- CreateSeuratObject(counts = Read10X(data.dir = "~/share_epronk/MDS_345_Niche_Immune/outs/filtered_feature_bc_matrix/") , min.cells = 3 , min.features = 200 , project = 'MDS345')
MDS10209_niche_immune <- CreateSeuratObject(counts = Read10X(data.dir = "~/share_epronk/MDS_10209_Niche_Immune/outs/filtered_feature_bc_matrix/") , min.cells = 3 , min.features = 200 , project = 'MDS10209')
MDS_SF3B1_niche_immune <- CreateSeuratObject(counts = Read10X(data.dir = "~/share_epronk/MDS_SF3B1_Niche_Immune/outs/filtered_feature_bc_matrix/") , min.cells = 3 , min.features = 200 , project = 'MDS_SF3B1')

Niche_MDS_individual<-list(MDS337_niche_immune,MDS345_niche_immune,MDS10209_niche_immune,MDS_SF3B1_niche_immune)

Niche_MDS_individual<-lapply(Niche_MDS_individual,function(x){
  x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000);
  x<- NormalizeData(x, normalization.method  = "LogNormalize", scale.factor = 10000); 
  x<-ScaleData(x,verbose = FALSE)
  x <- RunPCA(x, npcs = 30,verbose = FALSE);
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  x <- FindClusters(x, resolution = 0.3);
  x <- RunUMAP(x,dims = 1:20)
})


FeaturePlot(Niche_MDS_individual[[4]], reduction = "umap", label=T ,pt.size = 1,features=c("CXCL12"),split.by="orig.ident",min.cutoff = 0.25)



barcode_MDS337_niche<-CellSelector(DimPlot(Niche_MDS_individual[[1]],reduction="umap"))
barcode_MDS345_niche<-CellSelector(DimPlot(Niche_MDS_individual[[2]],reduction="umap"))
barcode_MDS10209_niche<-CellSelector(DimPlot(Niche_MDS_individual[[3]],reduction="umap"))
barcode_SF3B1_niche<-CellSelector(DimPlot(Niche_MDS_individual[[4]],reduction="umap"))



MDS337_niche<-subset(MDS337_niche_immune, cells=barcode_MDS337_niche)
MDS345_niche<-subset(MDS345_niche_immune,cells=barcode_MDS345_niche)
MDS10209_niche<-subset(MDS10209_niche_immune,cells=barcode_MDS10209_niche)
MDS_SF3B1_niche<-subset(MDS_SF3B1_niche_immune,cells=barcode_SF3B1_niche)


MDS337_niche[["percent.mt"]] <- PercentageFeatureSet(MDS337_niche, pattern = "^MT-");
MDS337_niche[["nRatio"]]=MDS337_niche$nCount_RNA/MDS337_niche$nFeature_RNA;    
VlnPlot(MDS337_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS337_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS337_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
plot(density(MDS337_niche$nRatio,to = 8));
MDS337_niche<-subset(MDS337_niche,subset = nRatio <=2.8 & percent.mt < 5);
VlnPlot(MDS337_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS337_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS337_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
plot(density(MDS337_niche$nRatio,to = 8));
MDS337_niche<- FindVariableFeatures(MDS337_niche, selection.method = "vst", nfeatures = 2000);
MDS337_niche<- NormalizeData(MDS337_niche, normalization.method  = "LogNormalize", scale.factor = 10000)

MDS345_niche[["percent.mt"]] <- PercentageFeatureSet(MDS345_niche, pattern = "^MT-");
MDS345_niche[["nRatio"]]=MDS345_niche$nCount_RNA/MDS345_niche$nFeature_RNA;    
VlnPlot(MDS345_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS345_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS345_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
plot(density(MDS345_niche$nRatio,to = 8));
MDS345_niche<-subset(MDS345_niche,subset = nRatio <=5 & percent.mt < 5);
VlnPlot(MDS345_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS345_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS345_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
MDS345_niche<- FindVariableFeatures(MDS345_niche, selection.method = "vst", nfeatures = 2000);
MDS345_niche<- NormalizeData(MDS345_niche, normalization.method  = "LogNormalize", scale.factor = 10000)



MDS10209_niche[["percent.mt"]] <- PercentageFeatureSet(MDS10209_niche, pattern = "^MT-");
MDS10209_niche[["nRatio"]]=MDS10209_niche$nCount_RNA/MDS10209_niche$nFeature_RNA;    
VlnPlot(MDS10209_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS10209_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS10209_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
plot(density(MDS10209_niche$nRatio,to = 8));
MDS10209_niche<-subset(MDS10209_niche,subset = nRatio <=3.6 & percent.mt < 5);
VlnPlot(MDS10209_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS10209_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS10209_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
MDS10209_niche<- FindVariableFeatures(MDS10209_niche, selection.method = "vst", nfeatures = 2000);
MDS10209_niche<- NormalizeData(MDS10209_niche, normalization.method  = "LogNormalize", scale.factor = 10000)


MDS_SF3B1_niche[["percent.mt"]] <- PercentageFeatureSet(MDS_SF3B1_niche, pattern = "^MT-");
MDS_SF3B1_niche[["nRatio"]]=MDS_SF3B1_niche$nCount_RNA/MDS_SF3B1_niche$nFeature_RNA;    
VlnPlot(MDS_SF3B1_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS_SF3B1_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS_SF3B1_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
plot(density(MDS_SF3B1_niche$nRatio,to = 8));
MDS_SF3B1_niche<-subset(MDS_SF3B1_niche,subset = nRatio <=3.2 & percent.mt < 5);
MDS_SF3B1_niche<-subset(MDS_SF3B1_niche,subset = nFeature_RNA <=2500& nCount_RNA < 6000);

VlnPlot(MDS_SF3B1_niche, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
CombinePlots(plots = list(FeatureScatter(MDS_SF3B1_niche, feature1 = "nCount_RNA", feature2 = "percent.mt"), FeatureScatter(MDS_SF3B1_niche, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")));
MDS_SF3B1_niche<- FindVariableFeatures(MDS_SF3B1_niche, selection.method = "vst", nfeatures = 2000);
MDS_SF3B1_niche<- NormalizeData(MDS_SF3B1_niche, normalization.method  = "LogNormalize", scale.factor = 10000)

MDS337_niche[["Group"]]<-"MDS"
MDS345_niche[["Group"]]<-"MDS"
MDS10209_niche[["Group"]]<-"MDS"
MDS_SF3B1_niche[["Group"]]<-"MDS"

Niche_MDS_individual<-list(MDS337_niche,MDS345_niche,MDS10209_niche,MDS_SF3B1_niche)

features_MDS <- SelectIntegrationFeatures(object.list = Niche_MDS_individual) #Integratation MDS
MDS_anchors <- FindIntegrationAnchors(object.list = Niche_MDS_individual, anchor.features = features_MDS, dims=1:30)
Niche_MDS <- IntegrateData(anchorset = MDS_anchors,dims=1:30,k.weight=45)

Niche_MDS[["Group"]]<-"MDS"

saveRDS(Niche_MDS,"~Niche_MDS.RDS")

#####Integration NBM MDS and AML 
Niche_NBM<-readRDS("~NBM_niche_merge.RDS")
Niche_AML<-readRDS("~AML_niche_merge.RDS")

NBM_AML_MDS_Niche<-list(Niche_NBM,Niche_AML,Niche_MDS) #Integratation Ctrl and AML
features_NBM_MDS_AML<-SelectIntegrationFeatures(object.list = NBM_AML_MDS_Niche)
Niche_NBM_MDS_AML_anchors <- FindIntegrationAnchors(object.list = NBM_AML_MDS_Niche, anchor.features = features_NBM_MDS_AML, dims=1:30)
Niche_NBM_MDS_AML <- IntegrateData(anchorset = Niche_NBM_MDS_AML_anchors)

Niche_NBM_MDS_AML[["Sample"]]<-Niche_NBM_MDS_AML$orig.ident
Niche_NBM_MDS_AML$Group[grepl("X",Niche_NBM_MDS_AML$Sample)]<-"AML"
Niche_NBM_MDS_AML$Group[grepl("NBM",Niche_NBM_MDS_AML$Sample)]<-"NBM"
Niche_NBM_MDS_AML$Group[grepl("MDS",Niche_NBM_MDS_AML$Sample)]<-"MDS"
Niche_NBM_MDS_AML$Group<-factor(Niche_NBM_MDS_AML$Group,levels=c("NBM","MDS","AML"))

saveRDS(Niche_NBM_MDS_AML,"~Niche_NBM_MDS_AML.RDS")

Niche_NBM_MDS_AML<-readRDS("scRNAseq_dataset/RDSfile/chenMDS/MDS single cell/~Niche_NBM_MDS_AML.RDS")
#Preprocessing

Niche_NBM_MDS_AML <- ScaleData(Niche_NBM_MDS_AML,verbose = FALSE);
Niche_NBM_MDS_AML <- RunPCA(Niche_NBM_MDS_AML, npcs = 30,verbose = FALSE);


Niche_NBM_MDS_AML<- FindVariableFeatures(Niche_NBM_MDS_AML, selection.method = "vst", nfeatures = 2000);

ElbowPlot(Niche_NBM_MDS_AML)
nDims=10
Niche_NBM_MDS_AML <- FindNeighbors(Niche_NBM_MDS_AML, reduction = "pca", dims = 1:nDims)
Niche_NBM_MDS_AML <- FindClusters(Niche_NBM_MDS_AML, resolution = 0.4);
set.seed(101)
Niche_NBM_MDS_AML <- Seurat::RunUMAP(Niche_NBM_MDS_AML,dims = 1:nDims);

DimPlot(Niche_NBM_MDS_AML, reduction = "umap",label=T,split.by = "Group" ,pt.size = 0.3)


Idents(Niche_NBM_MDS_AML)<-"seurat_clusters"

DefaultAssay(Niche_NBM_MDS_AML)<-"integrated"

FeaturePlot(Niche_NBM_MDS_AML, reduction = "umap", split.by="Group",features=c("NFKBIA"),sort=T,min.cutoff = 0.25)
DotPlot(Niche_NBM_MDS_AML,split.by = "Group", features=c("WNT4"),cols=rep("blue",3))
DotPlot(Niche_NBM_MDS_AML,group.by = "Sample", features=c("CXCL14"))

Niche_NBM_MDS_AML[["Clusters"]]<-Niche_NBM_MDS_AML$seurat_clusters

Niche_NBM_MDS_AML$Clusters<-ifelse(Niche_NBM_MDS_AML$seurat_clusters %in% c("1","3"),"BMSC-2",Niche_NBM_MDS_AML$Clusters)
Niche_NBM_MDS_AML$Clusters<-ifelse(Niche_NBM_MDS_AML$seurat_clusters %in% c("4"),"BMSC-0",Niche_NBM_MDS_AML$Clusters)
Niche_NBM_MDS_AML$Clusters<-ifelse(Niche_NBM_MDS_AML$seurat_clusters %in% c("2","0","5"),"BMSC-1",Niche_NBM_MDS_AML$Clusters)
Niche_NBM_MDS_AML$Clusters<-ifelse(Niche_NBM_MDS_AML$seurat_clusters %in% c("6"),"BMSC-3",Niche_NBM_MDS_AML$Clusters)
Niche_NBM_MDS_AML$Clusters<-factor(Niche_NBM_MDS_AML$Clusters,levels=c("BMSC-0","BMSC-1","BMSC-2","BMSC-3"))

saveRDS(Niche_NBM_MDS_AML,"scRNAseq_dataset/RDSfile/chenMDS/MDS single cell/~Niche_NBM_MDS_AML.RDS")

Idents(Niche_NBM_MDS_AML)<-"Clusters"
##################

#MDS vs. NBM
####################
Niche_NBM_MDS<-readRDS("scRNAseq_dataset/RDSfile/chenMDS/MDS single cell/~Niche_NBM_MDS.RDS")

Niche_NBM_MDS$Group<-ifelse(Niche_NBM_MDS$Group=="NBM","NBM","LRMDS")
Niche_NBM_MDS$Group<-factor(Niche_NBM_MDS$Group,levels=c("NBM","LRMDS"))

Niche_NBM_MDS$Clusters<-factor(Niche_NBM_MDS$Clusters,levels=c("BMSC-0","BMSC-1","BMSC-2","BMSC-3"))

DimPlot(Niche_NBM_MDS, reduction = "umap",label=F,split.by = "Group" ,pt.size = 0.3)

DimPlot(Niche_NBM_MDS, reduction = "umap",label=F,split.by = "Group" ,pt.size = 0.3,cols=c("grey60","darkgreen","tomato3","firebrick4"))

FeaturePlot(Niche_NBM_MDS, reduction = "umap",label=F,features=c("ICAM1","JAG1"),split.by = "Group",min.cutoff = 0.25,order=T)
FeaturePlot(Niche_NBM_MDS, reduction = "umap",label=F,features=c("NFKBIA"),split.by = "Group",min.cutoff = 0.25,order=T)

densityplot(Niche_NBM_MDS,assay="RNA",features="JAG1",split.by="Group")+facet_wrap(.~Group)


Niche_freq<-Niche_NBM_MDS@meta.data%>%group_by(Group,Clusters,Sample)%>%
  summarize(n=n())%>%group_by(Group,Sample)%>%summarize(Clusters=Clusters,per=n/sum(n)*100)
Sample1<-unique(Niche_freq$Sample)
Niche_freq<-Niche_freq%>%mutate(Sample=factor(Sample,levels=c(Sample1[grepl("NBM",Sample1)],Sample1[grepl("MDS",Sample1)])))

Niche_freq%>%ggplot(aes(Sample,per,fill=Clusters))+geom_bar(stat="identity",Position="stack")+theme_classic()+
  scale_fill_manual(values=c("steelblue4","darkgreen","tomato3","firebrick4","tan4","gray64","azure4","royalblue2"))+theme(axis.text.x=element_text(size=12,angle=45,hjust=0.8,vjust=0.8))

DotPlot(subset(Niche_NBM_MDS,Group=="MDS"),features=c("CXCL5"),assay="RNA",split.by = "Group",cols=c("blue","blue"))&
  theme(axis.title = element_blank())&coord_flip()&theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))&theme(legend.position = "right")


Niche_freq%>%filter(Clusters %in% c("BMSC-1","BMSC-2"))%>%ggplot(aes(Group,per,color=Group))+geom_boxplot(outlier.size = -1)+geom_jitter(width=0.3)+
  facet_wrap(.~Clusters,ncol=2)+theme_bw()+scale_color_manual(values=c("steelblue4","tomato4"))+
  theme(strip.background = element_blank())+NoLegend()+stat_compare_means(method = "t.test",label = "p.format")

Niche_freq_statistics<-Niche_freq%>%group_by(Group,Clusters)%>%summarise(Mean=mean(per),SD=sd(per))

Niche_freq_statistics_type<-Niche_freq%>%filter(Clusters=="BMSC-1")
t.test(Niche_freq_statistics_type$per[which(Niche_freq_statistics_type$Group=="NBM")],Niche_freq_statistics_type$per[which(Niche_freq_statistics_type$Group=="MDS")],var.equal = F)

VlnPlot(Niche_NBM_MDS,features=c("NFKBIA","ANXA1","CD44","CXCL8","CXCL12","KITLG","IL7","ANGPT1"),
        ,cols=c("steelblue4","tomato4"),ncol=2,pt.size=0,split.by = "Group", assay="RNA")&
  theme(axis.title = element_blank())

DotPlot(Niche_NBM_MDS,features=c("PDGFRA","LEPR","LPL","CXCL12","KITLG","ANGPT1","IL7","BGLAP","SPP1","SOX9","ACAN","COL2A1","S100A4","SEMA3C",
                                  "NES","CSPG4","ACTA2"),split.by = "Group",cols=c("blue","blue"))&
  theme(axis.title = element_blank())&coord_flip()&theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))&theme(legend.position = "right")


VlnPlot(subset(Niche_NBM_MDS,subset=Group=="NBM"),features=c("PDGFRA","LEPR","CXCL12","KITLG","ANGPT1","IL7"),
        ,cols=c("steelblue4"),ncol=2,pt.size=0,split.by = "Group", assay="RNA")&
  theme(axis.title = element_blank())

VlnPlot(Niche_NBM_MDS,features=c("PDGFRA"),
        ,cols=c("steelblue4","tomato4"),ncol=2,pt.size=0,group.by = "Group", assay="RNA")&
  theme(axis.title = element_blank())

MDS_NBM_BMSC_marker<-run_pseudobulkDE(meta_data=Niche_NBM_MDS@meta.data,expr_mat=Niche_NBM_MDS@assays$RNA@counts,
                                    contrasts='condition,MDS,NBM',design='~ condition',
                                    sample_col="Sample", condition_col="Group",cluster_col="Clusters" ,extr_col=NULL,do_prefiltering = TRUE,use_celltype = FALSE)
MDS_NBM_BMSC_marker["Gene"]<-row.names(MDS_NBM_BMSC_marker)
MDS_NBM_BMSC_marker["rank"]<- -sign(MDS_NBM_BMSC_marker$log2FoldChange) * log10(MDS_NBM_BMSC_marker$pvalue)


MDS_NBM_BMSC_marker_new<-run_pseudobulkDE_new(data=Niche_NBM_MDS,
                                      contrasts='condition,LRMDS,NBM',design='~ condition',
                                      sample_col="Sample", condition_col="Group",extr_col=NULL,do_prefiltering = TRUE)


MDS_NBM_BMSC_marker_new["Gene"]<-row.names(MDS_NBM_BMSC_marker_new)
MDS_NBM_BMSC_marker_new["rank"]<- -sign(MDS_NBM_BMSC_marker_new$log2FoldChange) * log10(MDS_NBM_BMSC_marker_new$pvalue)


MDS_NBM_BMSC_marker_new$pvalue<-ifelse(MDS_NBM_BMSC_marker_new$pvalue<1e-30,1e-30,MDS_NBM_BMSC_marker_new$pvalue)
MDS_NBM_BMSC_marker_new["Enrichment"]<-ifelse(MDS_NBM_BMSC_marker_new$rank>0,"MDS","NBM")
MDS_NBM_BMSC_marker_new["Enrichment"]<-factor(MDS_NBM_BMSC_marker_new$Enrichment,levels=c("NBM","MDS"))

MDS_NBM_BMSC_marker_sig<-MDS_NBM_BMSC_marker_new[which(MDS_NBM_BMSC_marker_new$padj<0.05),]

MDS_NBM_BMSC_marker_top<-MDS_NBM_BMSC_marker_sig%>%filter(abs(log2FoldChange)>0.585)

MDS_NBM_BMSC_seq_sig_lables<-MDS_NBM_BMSC_marker_sig%>%
  filter(Gene %in% c("BMP4","GADD45A","GADD45B","GDF15","IRF1","CDKN1A","BMP5","RELA","ANGPT1","KITLG","IL7","LEPR","LPL","ANXA1","NFKBIA","TNFAIP8","CD44","CXCL8","SOX9","FN1","COL5A3"))

MDS_NBM_BMSC_marker_new%>%ggplot(aes(log2FoldChange,-log10(padj)))+geom_point(size=0.1,col="grey")+
  geom_point(data=MDS_NBM_BMSC_marker_sig,mapping=aes(log2FoldChange,-log10(padj)),size=0.5,alpha=0.3,col="grey40")+
  geom_vline(xintercept= c(0.5,-0.5),col="grey",type=2)+geom_hline(yintercept = -log10(0.05),col="grey",type=2)+
  geom_text_repel(MDS_NBM_BMSC_seq_sig_lables,mapping=aes(log2FoldChange,-log10(padj),label=Gene,col=Enrichment),max.overlaps = 20)+
  theme_bw()+scale_color_manual(values=c("steelblue4","tomato4"))



write_xlsx(MDS_NBM_BMSC_marker,"MDS_NBM_BMSC_deseq.xlsx")


MDS_NBM_BMSC_marker_sig<-MDS_NBM_BMSC_marker[
  which(MDS_NBM_BMSC_marker$padj<0.05),]


rank_NBM_MDS_BMSC<-MDS_NBM_BMSC_marker_new$rank
names(rank_NBM_MDS_BMSC)<-row.names(MDS_NBM_BMSC_marker_new)
rank_NBM_MDS_BMSC<-rank_NBM_MDS_BMSC[is.finite(rank_NBM_MDS_BMSC)]

NBM_MDS_BMSC_Hall <- fgsea(hall_set, rank_NBM_MDS_BMSC, minSize=15, maxSize = 500, nperm=1000)
NBM_MDS_BMSC_Hall<-data.frame(NBM_MDS_BMSC_Hall)
NBM_MDS_BMSC_Hall["Enrichment"]<-ifelse(NBM_MDS_BMSC_Hall$NES>0,"MDS","NBM")
NBM_MDS_BMSC_Hall_sig<-NBM_MDS_BMSC_Hall%>%filter(padj<0.05)
NBM_MDS_BMSC_Hall_sig<-NBM_MDS_BMSC_Hall_sig[order(NBM_MDS_BMSC_Hall_sig$NES),]

NBM_MDS_BMSC_Hall_sig%>%group_by(Enrichment)%>%top_n(n=10,wt=abs(NES))%>%mutate(pathway=factor(pathway, levels=NBM_MDS_BMSC_Hall_sig$pathway))%>%
  ggplot(aes(NES,pathway,fill= padj))+geom_col()+theme_classic()+scale_fill_gradient(low="red",high="blue")

#Subpopulation

MDS_NBM_BMSC_0<-subset(Niche_NBM_MDS,subset=Clusters=="BMSC-0")


MDS_NBM_BMSC_marker_0<-run_pseudobulkDE_new(data=MDS_NBM_BMSC_0,
                                      contrasts='condition,LRMDS,NBM',design='~ condition',
                                      sample_col="Sample", condition_col="Group",cluster_col="Clusters" ,extr_col=NULL,do_prefiltering = TRUE)


MDS_NBM_BMSC_marker_0["Gene"]<-row.names(MDS_NBM_BMSC_marker_0)
MDS_NBM_BMSC_marker_0[MDS_NBM_BMSC_marker_0$Gene=="IL7",]


#Inflammatory cluster analysis in MDS

Niche_MDS<-subset(Niche_NBM_MDS,subset=Group=="MDS")
Niche_MDS[["ClusterGroup"]]<-ifelse(Niche_MDS$Clusters=="BMSC-2", "BMSC2","Others")

MDS_BMSC_marker<-run_pseudobulkDE(meta_data=Niche_MDS@meta.data,expr_mat=Niche_MDS@assays$RNA@counts,
                                      contrasts='condition,BMSC2,Others',design='~ condition + Extra',
                                      sample_col="Sample", condition_col="ClusterGroup",cluster_col="Clusters" ,extr_col="Sample",do_prefiltering = TRUE,use_celltype = FALSE)


MDS_BMSC_marker["Gene"]<-row.names(MDS_BMSC_marker)
MDS_BMSC_marker["rank"]<- -sign(MDS_BMSC_marker$log2FoldChange) * log10(MDS_BMSC_marker$pvalue)

MDS_BMSC_marker_new<-run_pseudobulkDE_new(data=Niche_MDS,
                                  contrasts='condition,BMSC2,Others',design='~ condition + Extra',
                                  sample_col="Sample", condition_col="ClusterGroup",extr_col="Sample",do_prefiltering = TRUE)

MDS_BMSC_marker_new["Gene"]<-row.names(MDS_BMSC_marker_new)

write_xlsx(MDS_BMSC_marker,"BMSC2_markergene.xlsx")

MDS_BMSC_marker_sig<-MDS_BMSC_marker[
  which(MDS_BMSC_marker$padj<0.05),]


MDS_BMSC_marker%>%ggplot(aes(-log10(padj),log2FoldChange))+geom_hline(yintercept= -0.5,linetype=3,color="grey")+
  geom_hline(yintercept=0.5,linetype=3,color="grey")+geom_vline(xintercept= -log10(0.05),linetype=3,color="grey")+
  geom_point(alpha=0.2)+
  geom_text_repel(data=MDS_BMSC_marker_sig,mapping=aes(-log10(padj),log2FoldChange,label=Gene),size=3,max.overlaps=30)+
  theme_classic()

rank_MDS_BMSC<-MDS_BMSC_marker$rank
names(rank_MDS_BMSC)<-row.names(MDS_BMSC_marker)
rank_MDS_BMSC<-rank_MDS_BMSC[is.finite(rank_MDS_BMSC)]

MDS_BMSC_Hall <- fgsea(hallmakers_set, rank_MDS_BMSC, minSize=15, maxSize = 500, nperm=1000)
MDS_BMSC_Hall_sig<-MDS_BMSC_Hall%>%filter(padj<0.05)%>%top_n(n=20,wt=NES)
MDS_BMSC_Hall_sig<-MDS_BMSC_Hall_sig[order(NES),]

MDS_BMSC_Hall_sig%>%mutate(pathway=factor(pathway, levels=MDS_BMSC_Hall_sig$pathway))%>%
  ggplot(aes(NES,pathway,fill= padj))+geom_col()+theme_classic()+scale_fill_gradient(low="red",high="blue")

#Comparision between BMSC subcluster in NBM
###############

Niche_NBM<-subset(Niche_NBM_MDS_AML,subset=Group=="NBM")
Niche_NBM[["ClusterGroup"]]<-ifelse(Niche_NBM$Clusters %in% c("BMSC-0","BMSC-1"), "BMSC01","BMSC23")

Niche_freq<-Niche_NBM@meta.data%>%group_by(Group,ClusterGroup,Sample)%>%
  summarize(n=n())%>%group_by(Group,Sample)%>%summarize(ClusterGroup=ClusterGroup,per=n/sum(n)*100)
Sample1<-unique(Niche_freq$Sample)
Niche_freq<-Niche_freq%>%mutate(Sample=factor(Sample,levels=c(Sample1[grepl("NBM",Sample1)],Sample1[grepl("MDS",Sample1)])))
Niche_freq%>%group_by(ClusterGroup)%>%summarize(Mean=mean(per),Sd=sd(per))

NBM_BMSC_marker<-run_pseudobulkDE_new(data=Niche_NBM,
                                  contrasts='condition,BMSC01,BMSC23',design='~ condition + Extra',
                                  sample_col="Sample", condition_col="ClusterGroup",cluster_col="Clusters" ,extr_col="Sample",do_prefiltering = TRUE)
NBM_BMSC_marker<-NBM_BMSC_marker[order(NBM_BMSC_marker$padj),]
NBM_BMSC_marker["Gene"]<-row.names(NBM_BMSC_marker)
NBM_BMSC_marker["rank"]<- -sign(NBM_BMSC_marker$log2FoldChange) * log10(NBM_BMSC_marker$pvalue)
write_xlsx(NBM_BMSC_marker,"NBM_BMSC-01_vs_34_markergene.xlsx")



NBM_BMSC_marker_sig<-NBM_BMSC_marker[
  which(NBM_BMSC_marker$padj<0.05),]


NBM_BMSC_marker%>%ggplot(aes(-log10(padj),log2FoldChange))+geom_hline(yintercept= -0.5,linetype=3,color="grey")+
  geom_hline(yintercept=0.5,linetype=3,color="grey")+geom_vline(xintercept= -log10(0.05),linetype=3,color="grey")+
  geom_point(alpha=0.2)+
  geom_text_repel(data=NBM_BMSC_marker_sig,mapping=aes(-log10(padj),log2FoldChange,label=Gene),size=3,max.overlaps=30)+
  theme_classic()

rank_NBM_BMSC<-NBM_BMSC_marker$rank
names(rank_NBM_BMSC)<-row.names(NBM_BMSC_marker)
rank_NBM_BMSC<-rank_NBM_BMSC[is.finite(rank_NBM_BMSC)]

NBM_BMSC_Hall <- fgsea(hallmakers_set, rank_NBM_BMSC, minSize=15, maxSize = 500, nperm=1000)
NBM_BMSC_Hall_sig<-NBM_BMSC_Hall%>%filter(padj<0.05)%>%top_n(n=20,wt=NES)
NBM_BMSC_Hall_sig<-NBM_BMSC_Hall_sig[order(NES),]

NBM_BMSC_Hall_sig%>%mutate(pathway=factor(pathway, levels=NBM_BMSC_Hall_sig$pathway))%>%
  ggplot(aes(NES,pathway,fill= padj))+geom_col()+theme_classic()+scale_fill_gradient(low="red",high="blue")

###############
#addmodelscore

hallmakers<-GSA.read.gmt("~/single_cell_sequencing/GSEA/hallmakers_genesymble.gmt")
hallmakers_set<-hallmakers$genesets
names(hallmakers_set)<-hallmakers$geneset.names

DefaultAssay(Niche_NBM_MDS)<-"RNA"
GSEA_score <- AddModuleScore(object = Niche_NBM_MDS, features = hallmakers_set[1], ctrl=50,name = 'TNFa_signaling')
GSEA_score <- AddModuleScore(object = GSEA_score, features = hallmakers_set[31], ctrl=50,name = 'Inflammatory_response')
GSEA_score <- AddModuleScore(object = GSEA_score, features = hallmakers_set[10], ctrl=50,name = 'Apoptosis')



Idents(GSEA_score)<-"Clusters"

FeaturePlot(GSEA_score, reduction = "umap", label=F,pt.size=0.5 ,min.cutoff = 0.05,cols=c("grey","red4"),ncol=4,order=T,
            features=c("Inflammatory_response1","TNFa_signaling1"),split.by="Group")&theme(axis.title.x=element_blank(),axis.title.y=element_blank())

VlnPlot(GSEA_score,features=c("Inflammatory_response1","TNFa_signaling1"),split.by = "Group", assay="RNA",pt.size=0,cols=c("steelblue4","tomato4"))

