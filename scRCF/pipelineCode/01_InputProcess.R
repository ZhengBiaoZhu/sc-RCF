# 01_inputProcess.R
# A step is to process the scRNA count matrix: 
# 1. Process the input data, convert each dataset into Seurat objects, 
# 2. perform differential gene expression analysis, 
# 3. and then conduct dimensionality reduction and clustering to visualize the results.
# 4. Next, use a list of TFs to filter the differentially expressed genes, retaining only the differentially expressed transcription factors.

library(Seurat)
### Load the scRNA-seq data from a specific path
## @scRNADataLoc: The location of the scRNA-seq data
## @scRNACelltype: Cell types in the input single-cell data
scRNALoad <- function(scRNADataLoc, scRNACelltype){
  scRNACount = readCountData(scRNADataLoc)
  scRNAData = CreateSeuratObject(scRNACount)
  scRNAData$celltype = scRNACelltype
  return(scRNAData)
}

### Read a raw count matrix (genes x cells) from a Seurat object, sparse matrix, data frame (.rds) or table (.txt/.csv)
## @scRNADataLoc: The location of the scRNA-seq data
readCountData <- function(scRNADataLoc){
  if(grepl('\\.rds$', scRNADataLoc, ignore.case = TRUE)){
    scRNACount = readRDS(scRNADataLoc)
    if(inherits(scRNACount, 'Seurat')){
      if(inherits(scRNACount[['RNA']], 'Assay5')){
        scRNACount = GetAssayData(JoinLayers(scRNACount[['RNA']]), layer = 'counts')
      }else{
        scRNACount = GetAssayData(scRNACount, assay = 'RNA', slot = 'counts')
      }
    }
  }else if(grepl('\\.txt$', scRNADataLoc, ignore.case = TRUE)){
    scRNACount = read.table(scRNADataLoc, sep = '\t', stringsAsFactors = F, header = T, row.names = 1, check.names = F)
  }else if(grepl('\\.csv$', scRNADataLoc, ignore.case = TRUE)){
    scRNACount = read.csv(scRNADataLoc, sep = ',', stringsAsFactors = F, header = T, row.names = 1, check.names = F)
  }else{
    stop("The input scRNA-seq data can only be in Seurat (.rds) or count matrix (.txt/.csv) format.")
  }
  return(scRNACount)
}

### Identify differentially expressed genes between initial and target cell types
## @iniCellData: The scRNA-seq data of initial cells
## @tarCellData: The scRNA-seq data of target cells
gainDEGs <- function(seu_merge, iniCellType, tarCellType){
  seu_merge[["percent.mt"]] <- PercentageFeatureSet(seu_merge, pattern = "^(MT|mt|Mt)-") # human genes start with "MT-", mouse with "mt-"
  seu_merge <- subset(seu_merge, subset = nFeature_RNA > 200 & percent.mt < 20)
  seu_merge = NormalizeData(seu_merge)
  if(packageVersion('Seurat') >= '5.0.0'){
    seu_merge = JoinLayers(seu_merge)
  }
  tarInitMarker =  FindMarkers(seu_merge, ident.1 = tarCellType, ident.2 = iniCellType, min.pct = 0.1, logfc.threshold = 0.25, only.pos = FALSE, group.by = 'celltype')
  return(tarInitMarker)
}

### Identify differentially expressed transcription factors from the differentially expressed genes
### The TFs information coming from TFDB_V3: https://guolab.wchscu.cn/AnimalTFDB#!/download
## @tarInitMarker: The differentially expressed genes bewteen initial and target cell types
## @spices: Species used in the experiment
gainDETFs <- function(tarInitMarker, spices, log2FC_t){
  if(spices == 'human'){
    TFs_list1 = read.table("../data/publicInfo/Homo_sapiens_TF.txt", header = T, sep = '\t')
    TFs_list2 = read.table("../data/publicInfo/Homo_sapiens_Cof.txt", header = T, sep = '\t')
  }else{
    TFs_list1 = read.table("../data/publicInfo/Mus_musculus_TF.txt", header = T, sep = '\t')
    TFs_list2 = read.table("../data/publicInfo/Mus_musculus_TF_cofactors.txt", header = T, sep = '\t')
  }
  TFs_all <- unique(c(TFs_list1$Symbol, TFs_list2$Symbol))
  DEGs_df_f = tarInitMarker[abs(tarInitMarker$avg_log2FC) > log2FC_t, ]
  DEGs_df_f$X = rownames(DEGs_df_f)
  DEGs_df_f = DEGs_df_f[DEGs_df_f$p_val_adj < 0.05, ]
  DEGs_df_f <- DEGs_df_f[DEGs_df_f$X %in% TFs_all, ]
  DEGs_df_f$Sign <- ifelse(DEGs_df_f$avg_log2FC > 0, 1, -1)
  DETFs = DEGs_df_f$Sign
  names(DETFs) <- DEGs_df_f$X
  return(DETFs)
}

### Obtain the count matrix of initial cell type
## @scRNADataLoc: The location of the scRNA-seq data
getCountMatrix <- function(scRNADataLoc){
  scRNACount = data.frame(as.matrix(readCountData(scRNADataLoc)), check.names = FALSE)
  return(scRNACount)
}

### Perform single-cell processing and dimensionality reduction clustering
## @seu_merge: The Seurat object contains initial cell type and target cell type
SeuratProcess <- function(seu_merge){
  library(harmony)
  seu_merge <- NormalizeData(seu_merge)
  seu_merge <- FindVariableFeatures(seu_merge, nfeatures = 3000)
  seu_merge <- ScaleData(seu_merge, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
  seu_merge <- RunPCA(seu_merge)
  seu_merge <- RunHarmony(seu_merge, group.by.vars = "celltype")
  seu_merge <- FindNeighbors(seu_merge, dims = 1:30, reduction = "harmony")
  seu_merge <- FindClusters(seu_merge, resolution = 0.3)
  seu_merge <- RunUMAP(seu_merge, dims = 1:30, reduction = "harmony")
  return(seu_merge)
}
