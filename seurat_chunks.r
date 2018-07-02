.libPaths(rev(.libPaths()))
library(Seurat)
library(dplyri)
library(Matrix)

load("~/dropseq/summary/R_Seurat_objects.rdata")
load("~/R_Seurat.rdata")


myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[sribo.gene.names, ])/col.total.umi, "pct.sribo")
gg <- VlnPlot(myumi,
              c("nUMI", "nGene", "top50", "umi.per.gene", "pct.Ribo", "pct.mito"),
              x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_UMI.pdf"),width=18,height=18)
ggsave(gg, file  = snakemake@output$pdf_violine, width = 18, height = 18)



# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


# standard log-normalization
myumi <- NormalizeData(myumi)

# choose ~1k variable genes
myumi <- FindVariableGenes(myumi, do.plot = FALSE, y.cutoff = 0.5)

# standard scaling (no regression)
myumi <- ScaleData(myumi, display.progress = FALSE)

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
myumi <- RunPCA(myumi, pcs.print = 0)
PCElbowPlot(myumi)

myumi <- FindClusters(myumi, dims.use = 1:2, k.param = 4, print.output = FALSE)
myumi <- RunTSNE(myumi, perplexity=7)

# in this plot, protein (ADT) levels are on top, and RNA levels are on the
# bottom
FeaturePlot(myumi, features.plot = c("pct.mito", "nUMI", "pct.Ribo", "top50"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 5.5, no.legend = F)












# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
gg <- GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
gg2 <- GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 15, genes.print = 5)
# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)

:

ggplot(t2, aes(x = nUMI, y = nGene)) +
  geom_point() +
  geom_smooth() +
coord_trans(x = "log10")

p1 <- TSNEPlot(pbmc, group.by = "louvain", do.return = TRUE)

##
