
seurat_preprocess <- function(obj) {
    NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData() |>
    RunPCA(features = VariableFeatures(object = x)) |>
    identity() 
}

seurat_cluster <- function(obj, res, dims) {
  FindNeighbors(obj, dims = 1:dims) |>
    FindClusters(resolution = res) |>
    RunUMAP(dims = 1:dims) |>
    identity()
}

seurat_jackstraw_plot <- function(obj) {
  JackStraw(obj, num.replicate = 100) |>
  ScoreJackStraw(dims = 1:20) |>
  JackStrawPlot(dims = 1:15)
}

seurat_add_metadata <- function(obj, metafile) {
  obj <- AddMetaData(obj, metafile)
}