CellChatDB.use <- CellChatDB

cellchat_create_from_seurat <- function(obj, sample_id) {
  sub_obj <- subset(obj, subset = sample_id == sample_id)
  cellChatobj <- createCellChat(object = sub_obj, meta = sub_obj@meta.data, group.by = "yulia_annot")
  cellChatobj@DB <- CellChatDB.use
  cellChatobj <- subsetData(cellChatobj)
  #future::plan("multiprocess", workers = 4)
  cellChatobj <- identifyOverExpressedGenes(cellChatobj)
  cellChatobj <- identifyOverExpressedInteractions(cellChatobj)
  # project gene expression data onto PPI network (optional)
  cellChatobj <- projectData(cellChatobj, PPI.human)
  cellChatobj <- computeCommunProb(cellChatobj)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellChatobj <- filterCommunication(cellChatobj, min.cells = 10)
  cellChatobj <- computeCommunProbPathway(cellChatobj)
  cellChatobj <- aggregateNet(cellChatobj)
  cellChatobj
}

SPP1_CD44_bubble <- function(cc_obj) {
    if("SPP1" %in% cc_obj@netP$pathways) {
    netVisual_bubble(cc_obj, remove.isolate = FALSE, pairLR.use = "SPP1_CD44",
                     + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
  }
}

