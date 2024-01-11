source("../scripts/setup.R")

run_clustering = function(celltype){
  set.seed(0)

  cat("Setting up data.\n")
  single_cell_experiments = set_up_sce(celltype)

  # Cluster etc
  gene_variance_model = scran::modelGeneVar(single_cell_experiments$rna_sce)
  plot(gene_variance_model$mean, gene_variance_model$total,
       xlab="Mean log-expression", ylab="Variance")
  curve(metadata(gene_variance_model)$trend(x), col="blue", add=TRUE)
  top.hvgs <- scran::getTopHVGs(gene_variance_model, n=2000)
  single_cell_experiments$rna_sce <- scater::runPCA(single_cell_experiments$rna_sce, subset_row=top.hvgs)
  cat("Running minibatch kmeans.\n")
  clusters <- mbkmeans::mbkmeans(single_cell_experiments$rna_sce,
                                 clusters = 100,
                                 reduceMethod = "PCA")
  colData(single_cell_experiments$rna_sce)[["cluster"]] = clusters$Clusters
  if(celltype=="pbmc"){
    pseudo_bulk = scater::sumCountsAcrossCells(single_cell_experiments$rna_sce, ids = single_cell_experiments$rna_sce[["cluster"]]) %>%
      assay("sum")
    celltype_by_cluster = -pseudo_bulk[pbmc_markers$marker,] %>%
      apply(1, scale) %>% set_rownames(colnames(pseudo_bulk)) %>%
      apply(1, rank) %>%  set_colnames(colnames(pseudo_bulk)) %>%
      reshape2::melt() %>%
      set_colnames(c("marker", "cluster", "rank")) %>%
      subset(rank<=2) %>%
      merge(pbmc_markers) %>% 
      dplyr::group_by(cluster) %>%
      dplyr::summarise(celltypes = paste(names(table(celltype)[table(celltype)==max(table(celltype))]), collapse = ","), markers = paste(marker, collapse = ",")) 
    write.csv(celltype_by_cluster, "pbmc_celltype_by_cluster.csv")
    markers_by_cluster = pseudo_bulk[rowSums(pseudo_bulk)>100, ] %>%
      add(1) %>%
      log2 %>%
      apply(1, scale, center = TRUE, scale = FALSE) %>% set_rownames(colnames(pseudo_bulk)) %>%
      t %>%
      reshape2::melt() %>%
      set_colnames(c("marker", "cluster", "log2fc")) %>%
      dplyr::group_by(cluster) %>% 
      dplyr::top_n(20) %>%
      dplyr::arrange(cluster, -log2fc)
    write.csv(markers_by_cluster, "pbmc_markers_by_cluster.csv")
    celltype_by_cluster = setNames(celltype_by_cluster$celltypes, celltype_by_cluster$cluster)
    colData(single_cell_experiments$rna_sce)[["celltype"]] = celltype_by_cluster[colData(single_cell_experiments$rna_sce)[["cluster"]]]
  }
  
  cat("Running tsne. \n")
  single_cell_experiments$rna_sce <- scater::runTSNE(single_cell_experiments$rna_sce,
                                                          num_dim = 2,
                                                          dimred = "PCA",
                                                          n_dimred = 50)

  # Add same clusters to the ATAC data
  colData(single_cell_experiments$atac_sce) =
    merge(colData(single_cell_experiments$atac_sce),
          colData(single_cell_experiments$rna_sce)[c("atac.bc", "cluster")],
          by = "atac.bc",
          all.x = T,
          all.y = F)

  # Examine results
  cat("Making metacells and saving various outputs. \n")
  dir.create(create_pseudobulk_path(celltype), recursive = T, showWarnings = F)

  withr::with_dir(create_pseudobulk_path(celltype), {
    dir.create("description", recursive = F, showWarnings = F)
    for( colour_by in rev( c( "cluster", "celltype", "sizeFactor", pbmc_markers$marker) ) ){
      try({scater::plotTSNE(single_cell_experiments$rna_sce) + ggplot2::coord_fixed()})
      try({scater::plotTSNE(single_cell_experiments$rna_sce, colour_by = colour_by) + ggplot2::coord_fixed()})
      ggplot2::ggsave(paste0("description/", colour_by, ".png"), width = 10, height = 6)
    }

    write.csv(
      Reduce(f=cbind,
             list(
               colData(single_cell_experiments$rna_sce),
               reducedDim(single_cell_experiments$rna_sce, "PCA"),
               reducedDim(single_cell_experiments$rna_sce, "TSNE")
             )
      ),
      "rna_metadata.csv"
    )
    pseudo_bulk = scater::sumCountsAcrossCells(single_cell_experiments$rna_sce, ids = single_cell_experiments$rna_sce[["cluster"]])
    write.csv(assay(pseudo_bulk, "sum"), "rna_pseudo_bulk.csv")

    # per-cluster variances
    clusters = read.csv("rna_metadata.csv")[["cluster"]]
    do_one_cluster = function(cluster_idx){
      totals = assay(single_cell_experiments$rna_sce, "logcounts")[, 2] %>%
        raise_to_power(2, .) %>%
        magrittr::subtract(1) %>%
        sum
      factor_to_get_back_to_cpm = 1e6 / totals
      if(sum(clusters==cluster_idx)>1){
        cluster_variances = assay(single_cell_experiments$rna_sce, "logcounts")[, clusters==cluster_idx] %>%
          raise_to_power(2, .) %>%
          magrittr::subtract(1) %>%
          multiply_by(factor_to_get_back_to_cpm) %>%
          DelayedMatrixStats::rowVars()
      } else {
        rep(NA, dim(single_cell_experiments$rna_sce)[[1]])
      }
    }
    pseudobulk_vars = sapply(1:100, do_one_cluster)
    write.csv(pseudobulk_vars, "rna_pseudo_bulk_vars.csv")


    # Log some basic metacell composition info
    cluster_composition =
      table(colData(single_cell_experiments$rna_sce)[["cluster"]],
            colData(single_cell_experiments$rna_sce)[["celltype"]])
    image(cluster_composition)
    metacell_metadata = data.frame(
      cluster = rownames(cluster_composition),
      total_cell_count            = apply(cluster_composition, 1, sum) ,
      largest_subpopulation_count = apply(cluster_composition, 1, max),
      largest_subpopulation       = apply(cluster_composition, 1, function(x) names(which.max(x))) ,
      total_cell_count_atac = c(table(single_cell_experiments$atac_sce$cluster)[rownames(cluster_composition)])
    )
    {
      pdf("pseudo_bulk_homogeneity.pdf")
      plot(largest_subpopulation_count ~ total_cell_count,
           data = metacell_metadata,
           main = "Pseudo-bulk internal homogeneity")
      abline(a=0, b=1)
      plot(total_cell_count_atac ~ total_cell_count,
           data = metacell_metadata,
           main = "Pseudo-bulk RNA vs ATAC counts")
      abline(a=0, b=1)
      dev.off()
    }
    write.csv(metacell_metadata, "metacell_metadata.csv")

    metacell_metadata = read.csv("metacell_metadata.csv", row.names = 1)

    # Per cluster sums for ATAC
    pseudo_bulk_atac = scater::sumCountsAcrossCells(single_cell_experiments$atac_sce,
                                                    ids = single_cell_experiments$atac_sce[["cluster"]])
    write.csv(assay(pseudo_bulk_atac, "sum"), "atac_pseudo_bulk.csv")
  })
}

run_clustering(celltype = "pbmc")
run_clustering(celltype = "skin")
run_clustering(celltype = "keratinocyte")
cat("Done.")