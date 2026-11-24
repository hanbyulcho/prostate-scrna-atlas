sc.proc = function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE, assay = "RNA") {
  DefaultAssay(sc) = assay
  if (normalize) {
    sc <- NormalizeData(sc, verbose = FALSE)
  }
  if (hvg) {
    sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = n, verbose = FALSE)
  }
  if (pca) {
    sc <- ScaleData(sc, verbose = FALSE)
    sc <- RunPCA(sc, npcs = pc, verbose = FALSE)
  }
  if (umap) {
    sc <- RunUMAP(sc, reduction = "pca", dims = 1:pc)
  }
  #sc <- RunTSNE(sc, reduction = "pca", dims = 1:pc)
  if (cluster) {
    sc <- FindNeighbors(sc, dims=1:pc)
    sc <- FindClusters(sc, resolution = res)
  }
  return(sc)
}
get.psk = function(obj, anno, sample, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  require(edgeR)
  obj$anno_sample = paste(as.character(obj@meta.data[,anno]), obj@meta.data[,sample], sep = ".")
  a = table(obj$anno_sample)
  a = a[a>min.cell]
  obj = obj[, obj$anno_sample %in% names(a)]
  whole.psk = Seurat::AggregateExpression(object = obj, assays = assay, return.seurat = F, group.by = "anno_sample")
  # whole.psk = lapply(names(a), function(x) {
  #   if (sum(obj$anno_sample==x)>1) {
  #     s = rowSums(obj[[assay]]$counts[,obj$anno_sample==x])
  #   } else {
  #     s = obj[[assay]]$counts[,obj$anno_sample==x]
  #   }
  #   return(s)
  # })
  # whole.psk = do.call(cbind, whole.psk)
  whole.psk = whole.psk[[1]]
  colnames(whole.psk) = names(a)
  if (keep.all.genes) {
    whole.psk = DGEList(counts = whole.psk)
  } else {
    whole.psk = DGEList(counts = whole.psk[rowSums(whole.psk)>0, ])
  }
  
  whole.psk = calcNormFactors(whole.psk)
  logtpm = log2(edgeR::cpm(whole.psk)+1)
  s = sapply(colnames(whole.psk), function(x) strsplit(x, "[.]")[[1]])
  s = as.data.frame(t(s))
  colnames(s) = c("anno", "sample")
  whole.psk$samples = cbind.data.frame(whole.psk$samples, s)
  return(list(dgelist = whole.psk, logtpm = logtpm))
}
process_singleR <- function(x, ref, labels) {
  print(names(x))
  p = SingleR(test = x[[2]], ref = ref[[1]], labels = labels[[2]], fine.tune = TRUE)
  s = lapply(1:nrow(p$scores), function(i) sort(p$scores[i,], decreasing = TRUE)[1:2])
  s1 = do.call(rbind, lapply(s, function(y) names(y)))
  s2 = do.call(rbind, s)
  s = as.data.frame(cbind(s1, s2))
  colnames(s) = c("hit1", "hit2", "hit1.score", "hit2.score")
  s$label = p$labels
  s$label.prune = p$pruned.labels
  s$label.final = s$hit1
  s$label.final[s$hit1 != s$label] = NA
  s = cbind(x[[1]]$samples, s)
  return(list(result = p, summary = s))
}
