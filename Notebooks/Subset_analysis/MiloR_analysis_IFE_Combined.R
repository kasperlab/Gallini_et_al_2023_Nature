#MiloR on Sara's combined IFE dataset

#https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#1_Load_data

#Load libraries
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(zellkonverter)

#Path to saved and annotated h5ad file from Scanpy analysis
h5_path <- "G:/My Drive/MKA/Papers/Greco/Sara and Nur/Loom files/Processed/"
filename = 'IFE_Combined_datasets_Annotated.h5ad'
file = paste(h5_path, filename, sep = '')

#Read in file with readH5AD from zellkonverter package
adata = readH5AD(file, X_name = 'raw', layers = FALSE)


#Plot umaps to check if all went fine
plotReducedDim(adata, colour_by="leiden_comb", dimred = "X_umap")

adata@int_colData@listData$reducedDims$umap = adata@int_colData@listData$reducedDims$X_umap
milo@int_colData@listData$reducedDims$UMAP = milo@int_colData@listData$reducedDims$X_umap

#Begin  for differential abundance testing
##Add milo to SingleCellExperiment object
milo <- Milo(adata)
print(milo)

##Construct KNN graph
###Use the precomputed PCA from scanpy for the graph. 
k = 20
d = 30
milo <- buildGraph(milo, k = k, d = d, reduced.dim = "X_pca")

##Defining representative neighbourhoods on the KNN graph
milo <- makeNhoods(milo, prop = 0.2, k = k, d=d, refined = TRUE, reduced_dims = "X_pca")

##As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples. In this case 5*2 conditions = 10
plotNhoodSizeHist(milo)

##Counting cells in neighbourhoods. Count how many cells from biological replicates are in each neighborhood
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="orig.ident")
head(nhoodCounts(milo))

##Defining experimental design
###We want to detect differential abundance between the wound states stored in 'injury' and sample name in "orig.ident"
design <- data.frame(colData(milo))[,c("orig.ident", "injury")]

#Convert columns from integer to factor
design <- distinct(design)
rownames(design) <- design$orig.ident

##Computing neighbourhood connectivity
###Note, most time-consuming step, could take some minutes for large datasets
milo <- calcNhoodDistance(milo, d=d, reduced.dim = "X_pca")

##DA-testing
da_results <- testNhoods(milo, design = ~ injury, design.df = design, reduced.dim = 'X_pca')
head(da_results)

da_results %>%
  arrange(- SpatialFDR) %>%
  head()


#Inspecting DA testing results

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "X_umap", colour_by="injury", text_by = "batch", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="X_umap",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


da_results <- annotateNhoods(milo, da_results, coldata_col = "orig.ident")
da_results <- annotateNhoods(milo, da_results, coldata_col = "batch")
#head(da_results)

#Automatic grouping of neighbourhoods
milo <- buildNhoodGraph(milo)

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 2 )
head(da_results)

plotNhoodGroups(milo, da_results, layout="umap") 

plotDAbeeswarm(da_results, "NhoodGroup")

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 1)
plotDAbeeswarm(da_results, "NhoodGroup")

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 0.5)
plotDAbeeswarm(da_results, "NhoodGroup")

plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 3) , group.by = "NhoodGroup") 

plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 3, overlap = 1) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 3, overlap = 5) , group.by = "NhoodGroup") 

plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 1, overlap = 5) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 1, overlap = 1) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 0.5, overlap = 1) , group.by = "NhoodGroup") 

plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 1) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 5) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 3) , group.by = "NhoodGroup") 
plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 2) , group.by = "NhoodGroup") 

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 3, overlap=1)
plotNhoodGroups(milo, da_results, layout="umap") + scale_color_brewer(palette="tab20")
plotDAbeeswarm(da_results, group.by = "NhoodGroup")

