#' wl-02-07-2020, Thu: Put package files together
#' wl-06-07-2020, Mon:
#'  - Find '<<-'. Should be '<-'
#'  - Global data sets: data_GOslim and data_ORF2KEGG are used in
#'    gene_clustering.
#' wl-12-07-2020, Sun: Find where NAs come from (caused by reshape2:dcast)
#' wl-13-07-2020, Mon: Handle NAs in pre_processing
#' wl-14-07-2020, Tue:
#'   - Fix a bug in gene_network.
#'   - Option for PCA computation: use R core package stats
#' wl-17-07-2020, Fri: Find more bugs in network analysis.
#' wl-22-07-2020, Wed: Re-write pre_processing function. Add option to omit NAs
#'   in wide format.
#' wl-23-07-2020, Thu: Debug function for converting long to widedata format.
#'   Also check NAs in wide format. If use 'length', the zero is NAs in
#'   wide format.
#' wl-25-07-2020, Sat: Simplify gene_clustering
#' wl-27-07-2020, Mon: Simplify GeneNetWork
#' wl-30-07-2020, Thu:
#'   - gene_network: Keep the largest cluster
#'   - pre_processing: Change user provided sd format. Two columns: Ion and sd
#'   - Test on data subset: fewer batches and fewer Ion items
#' wl-04-08-2020, Tue: re-write gene_clustering. Correct one mistake
#' wl-14-08-2020, Fri: remove two R packages "pheatmap", "qgraph"
#' wl-01-09-2020, Tue: Change data prparation in pre_processing
#' wl-02-09-2020, Wed: change variable names in pre_processing

#' =======================================================================
#'
pre_processing <- function(data = NULL, stdev = NULL, var_id = 1,
                           batch_id = 2, data_id = 3) {

  ## -------------------> Import data
  ## get raw data stats summary
  ## wl-01-09-2020, Tue: update data set
  data <- data[, c(var_id, batch_id, data_id:ncol(data))]
  names(data)[1:2] <- c("Knockout", "Batch_ID")
  mat <- data[, -c(1:2)]
  #' mat <- data[, -c(1:(data_id - 1))]

  res <- as.data.frame(t(sapply(mat, function(x) {
    c(round(summary(x), 3), round(var(x), 3))
  })))
  names(res)[ncol(res)] <- "Variance"
  res <- cbind(Ion = names(mat), res)
  rownames(res) <- NULL
  df_raw <- res

  ## -------------------> Outlier detection
  ## wl-22-07-2020, Wed: more general for Ion contents
  data_long <- reshape2::melt(data,
    id = c("Knockout", "Batch_ID"),
    variable.name = "Ion",
    value.name = "Concentration"
  )

  ## wl-06-07-2020, Mon: convert to factors before using levels function.
  ## wl-30-07-2020, Thu: replace 'as.factor' with 'factor' in case level
  ##  updating
  data_long$Knockout <- factor(data_long$Knockout)
  data_long$Ion <- factor(data_long$Ion)
  ion_name <- levels(data_long$Ion)

  #' wl-23-07-2020, Thu: get Knockout outliers based on Ion
  data_long <- plyr::ddply(data_long, "Ion", function(x) {
    lowerq <- quantile(x$Concentration, na.rm = T)[2]
    upperq <- quantile(x$Concentration, na.rm = T)[4]
    iqr <- upperq - lowerq
    extreme_upper <- (iqr * 3) + upperq
    extreme_lower <- lowerq - (iqr * 3)
    x$Outlier <- ifelse((x$Concentration > extreme_upper) |
      (x$Concentration < extreme_lower), 1, 0)
    return(x)
  })

  df_outlier <-
    data.frame(cbind(
      levels(data_long$Ion),
      table(data_long$Ion, data_long$Outlier),
      round(table(data_long$Ion, data_long$Outlier)[, 2] /
        dim(data_long)[1] * 100, 2)
    ))
  rownames(df_outlier) <- c()
  colnames(df_outlier) <- c("Ion", "no_outlier", "outlier", "outlier(%)")

  #' wl-23-07-2020, Thu: remove Knockout outliers in terms of Ion.
  data_long <- data_long[data_long$Outlier < 1, ]
  data_long <- subset(data_long, select = -Outlier)
  ## wl-23-07-2020, Thu: NAs in wide format due to outlier removal.
  ## con.tab(data_long)  #' NAs: 28 in 1454 Knockout

  ## -------------------> Median batch correction
  data_long$log <- log(data_long$Concentration)

  #' wl-23-07-2020, Thu: remove median of each batch in each Ion
  data_long <- plyr::ddply(data_long, "Ion", function(x) {
    res <- plyr::ddply(x, "Batch_ID", function(y) {
      med <- median(y$log)
      y$log_corr <- y$log - med
      y
    })
  })

  #' get stats of log_corr
  res <- plyr::ddply(data_long, "Ion", function(x) {
    c(round(summary(x$log_corr), 3), round(var(x$log_corr), 3))
  })
  names(res)[ncol(res)] <- "Variance"
  df_bat <- res

  ## -------------------> Standardisation

  #' wl-08-07-2020, Wed: Use plyr::ddplyr. sds is for Ion
  if (is.null(stdev)) {
    sds <- plyr::ddply(data_long, "Ion", function(x) sd(x$log_corr))
    nam <- sds[, 1]
    sds <- as.numeric(as.vector(sds[, 2]))
    names(sds) <- nam
  } else if (ncol(stdev) == 1) {
    ## wl-30-07-2020, Thu: do NOT use one column. Alway with two columns:
    ## Ion and std.
    sds <- as.numeric(as.vector(stdev[, 1]))
    names(sds) <- ion_name ## problem if the size is not consistent.
  } else {
    sds <- stdev
    nam <- sds[, 1]
    sds <- as.numeric(as.vector(sds[, 2]))
    names(sds) <- nam
  }

  #' wl-21-07-2020, Tue: Normalise corr based ion std. Factor always gives
  #' trouble
  dat <- data_long[, c("Ion", "log_corr")]
  dat$Ion <- as.character(dat$Ion)
  tmp <- apply(dat, 1, function(x) {
    idx <- as.character(x[1]) == names(sds)
    as.numeric(x[2]) / sds[idx]
  })
  data_long <- cbind(data_long, log_corr_norm = tmp)

  #' really need to sort? (keep consistent with original code)
  data_long <- data_long[order(data_long$Knockout), ]
  rownames(data_long) <- NULL

  #' wl-22-07-2020, Wed: get summary of log_corr_norm
  res <- plyr::ddply(data_long, "Ion", function(x) {
    c(round(summary(x$log_corr_norm), 3), round(var(x$log_corr_norm), 3))
  })
  names(res)[ncol(res)] <- "Variance"
  df_std <- res

  ## -------------------> symbolization
  data_long$symb <-
    ifelse((data_long$log_corr_norm > -3) & (data_long$log_corr_norm < 3),
      0, ifelse(data_long$log_corr_norm >= 3, 1, -1)
    )

  ## -------------------> Aggregation of the batch replicas
  ## wl-13-07-2020, Mon: add prefix and change * as +
  dat <- data_long[, c("Knockout", "Ion", "log_corr_norm", "symb")]
  data_long_unique <-
    data.frame(stats::aggregate(. ~ Knockout + Ion, dat, median))
  ## update symb
  data_long_unique$symb <-
    ifelse((data_long_unique$symb < 0.5) & (data_long_unique$symb > -0.5),
      0, ifelse(data_long_unique$symb >= 0.5, 1, -1)
    )

  #' wl-23-07-2020, Thu: The missing values are from outlier detection
  #' wl-13-07-2020, Mon: Fill in structural(aggregation) missing values
  data_wide_unique <-
    reshape2::dcast(data_long_unique, Knockout ~ Ion,
      # fill = 0, #' wl: keep it or not?
      fun.aggregate = mean,
      value.var = "log_corr_norm"
    )
  data_wide_unique_symb <-
    reshape2::dcast(data_long_unique, Knockout ~ Ion,
      # fill = 0, #' wl: keep it or not?
      fun.aggregate = mean,
      value.var = "symb"
    )

  #' sum(is.na(data_wide_unique))
  #' dim(data_wide_unique)

  #' remove NAs
  data_wide_unique <- na.omit(data_wide_unique)
  data_wide_unique_symb <- na.omit(data_wide_unique_symb)

  #' wl-02-09-2020, Wed: change 'log_corr' to 'log_corr_norm'
  p1 <-
    ggplot(
      data = data_long,
      aes(
        x = factor(Batch_ID), y = log_corr_norm,
        col = factor(Batch_ID)
      )
    ) +
    geom_point(shape = 1) +
    facet_wrap(~Ion) +
    xlab("Batch.ID") +
    ylab("log(Concentration) (ppm)") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())

  p2 <-
    ggplot(data = data_long_unique, aes(x = log_corr_norm)) +
    geom_histogram(binwidth = .1) +
    facet_wrap(~Ion) +
    xlab("log(Concentration) (z-score)") +
    ylab("Frequency") +
    geom_vline(xintercept = c(-3, 3), col = "red")

  ## -------------------> Output
  res <- list()
  res$stats_raw_data <- df_raw # raw data
  res$stats_outliers <- df_outlier # outliers
  res$stats_batch_data <- df_bat # batch corrected data
  res$stats_stand_data <- df_std # standardised data
  res$data_long_bat <- data_long # with Batch_ID
  res$data_long <- data_long_unique # without Batch_ID
  res$data_wide <- data_wide_unique
  res$data_wide_symb <- data_wide_unique_symb
  res$plot_dot <- p1
  res$plot_hist <- p2
  return(res)
}

#' =======================================================================
#'
exploratory_analysis <- function(data = NULL) {

  ## -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4", "white", "firebrick"))

  corrplot.mixed(cor(data[, -1], use = "complete.obs"),
    number.cex = .7,
    lower.col = "black", upper.col = col3(100)
  )
  p_corr <- recordPlot()

  ## -------------------> PCA
  #' wl-14-07-2020, Tue: Original (trust) pca computation if there is no NAs.
  dat <- t(data[, -1])
  pca <- prcomp(dat, center = T, scale. = F)

  #' variance explained
  vars <- pca$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  dfn <- paste(names(vars), " (", vars[names(vars)], "%)", sep = "")

  #' PCA scores
  pca_scores <- data.frame(pca$x)
  #' names(pca_scores) <- dfn

  #' PCA loading
  pca_loadings <- data.frame(pca$rotation)
  rownames(pca_loadings) <- data$Knockout
  pca_loadings <- pca_loadings[, 1:2]

  #' PCA plot using ggplot2
  pca_p <-
    ggplot(data = pca_scores[, 1:2], aes(x = PC1, y = PC2)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.4) +
    geom_text_repel(aes(label = row.names(pca_scores)), size = 4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(dfn[1]) +
    ylab(dfn[2]) +
    labs(title = "PCA")

  ## pca_p

  ## -------------------> HEATMAP
  ## wl-14-08-2020, Fri:  use ggplots
  heatmap.2(as.matrix(data[, -1]),
    scale = "row", col = bluered(100),
    trace = "none", density.info = "none",
    hclustfun = function(x) hclust(x, method = "ward.D")
  )
  ## library(pheatmap)
  ## pheatmap(data[, -1], show_rownames = F, cluster_cols = T, cluster_rows = T,
  ##          legend = T, fontsize = 15, clustering_method = "ward.D",
  ##          scale = "row")
  pheat <- recordPlot()

  ## -------------------> PAIRWISE CORRELATION MAP
  col <- colorRampPalette(c("skyblue4", "white", "plum4"))(20)
  corr <- cor(na.omit(data[, -1]))
  heatmap(
    x = corr, col = col, symm = TRUE, cexRow = 1.4, cexCol = 1.4,
    main = ""
  )
  pcm <- recordPlot()

  #' -------------------> Regularized partial correlation network MAP
  #' wl-13-08-2020, Thu: there is no 'qgraph' in conda forge and bio conda.
  #' Have to plot the correlation network instead.
  if (T) {
    ## wl-14-08-2020, Fri: debug code only
    ## library(glasso)
    ## corr_reg <- glasso(corr, rho = 0.01)
    ## net <- network::network(corr_reg$w, directed = FALSE)
    net <- network::network(corr, directed = FALSE)

    #' set edge attributes
    net %e% "weight" <- corr
    net %e% "weight_abs" <- abs(corr) * 6
    net %e% "color" <- ifelse(net %e% "weight" > 0, "lightgreen", "coral")
    ## set.edge.value(net, "weight", corr)
    ## list.network.attributes(net)
    ## list.edge.attributes(net)
    ## list.vertex.attributes(net)
    net_p <-
      ggnet2(net,
        label = TRUE, mode = "spring",
        node.size = 10, edge.size = "weight_abs", edge.color = "color"
      )
    ## net_p
  } else {
    ## wl-06-07-2020, Mon: 'cor_auto' is from package qgraph(lavaan)
    ## wl-28-07-2020, Tue: cad and corr are the same
    cad <- cor_auto(data[, -1])
    suppressWarnings(qgraph(cad,
      graph = "glasso", layout = "spring",
      tuning = 0.25, sampleSize = nrow(data[, -1])
    ))
    Graph_lasso <- recordPlot()
  }

  ## -------------------> Output
  res <- list()
  res$plot_Pearson_correlation <- p_corr
  res$plot_pca_Individual <- pca_p
  res$data_pca_loadings <- pca_loadings
  res$plot_heatmap <- pheat
  res$plot_pairwise_correlation_map <- pcm
  ## res$plot_regularized_partial_correlation_network <- Graph_lasso
  res$plot_correlation_network <- net_p
  return(res)
}

#' =======================================================================
#'
gene_clustering <- function(data = NULL, data_symb = NULL, thres_clus = 10,
                           thres_anno = 5) {

  ## -------------------> Define clusters
  res_dist <- dist(data_symb[, -1], method = "manhattan") #' "euclidean"
  res_hc <- hclust(d = res_dist, method = "single")
  clus <- cutree(res_hc, h = 0) # distance 0

  data_symb$cluster <- clus

  ## -------------------> Subset cluster with more than 10 genes
  df <- as.data.frame(table(clus), stringsAsFactors = F)
  names(df) <- c("cluster", "nGenes")
  df_sub <- df[df$nGenes > thres_clus, ]
  rownames(df_sub) <- c()

  #' wl-24-07-2020, Fri: cluster index satisfing threshold of cluster number
  idx <- clus %in% df_sub$cluster
  #' sum(idx)

  mat <- data[idx, ]
  mat$cluster <- clus[idx]

  mat_long <-
    reshape2::melt(mat,
      id = c("Knockout", "cluster"), variable.name = "Ion",
      value.name = "log_corr_norm"
    )

  res <- sapply(mat_long$cluster, function(x) {
    tmp <- df_sub[df_sub$cluster == x, ]
    tmp <- paste0("Cluster ", tmp[1], " (", tmp[2], " genes)")
  })
  mat_long$cluster <- res #' update cluster with gene numbers

  clus_p <-
    ggplot(
      data = mat_long,
      aes(x = Ion, y = log_corr_norm)
    ) +
    facet_wrap(~cluster) +
    geom_line(aes(group = Knockout)) +
    stat_summary(fun.data = "mean_se", color = "red") +
    labs(x = "", y = "") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 10)
    )

  ## -------------------> KEGG AND GO SLIM ANNOTATION
  mat <- data_symb[idx, ]
  data_GOslim$Ontology <- as.character(data_GOslim$Ontology)

  kego <- plyr::dlply(mat, "cluster", function(x) {
    ## x <- subset(mat, cluster == "15")
    inputGeneSet <- as.character(x$Knockout)
    N <- length(inputGeneSet)

    res <- data_GOslim %>%
      dplyr::mutate(Ontology = setNames(
        c(
          "Biological process",
          "Cellular component",
          "Molecular function"
        ),
        c("P", "C", "F")
      )[Ontology]) %>%
      dplyr::filter(ORFs %in% inputGeneSet) %>%
      dplyr::group_by(GOslim, Ontology) %>%
      dplyr::filter(GOslim != "other") %>%
      dplyr::rename(Term = GOslim) %>%
      dplyr::summarise(Count = n()) %>%
      dplyr::mutate(Percent = Count / N * 100) %>%
      dplyr::bind_rows(data_ORF2KEGG %>%
        dplyr::filter(ORF %in% inputGeneSet) %>%
        dplyr::group_by(KEGGID, Pathway) %>%
        dplyr::summarise(Count = n()) %>%
        dplyr::mutate(Ontology = "KEGG") %>%
        dplyr::rename(Term = Pathway) %>%
        dplyr::ungroup() %>%
        dplyr::filter(KEGGID != "01100") %>%
        dplyr::select(-KEGGID) %>%
        dplyr::mutate(Percent = Count / N * 100)) %>%
      dplyr::filter(!Term %in% c(
        "molecular_function", "biological_process",
        "cellular_component"
      ))
  })
  names(kego) <- paste0("Cluster ", df_sub[[1]], " (", df_sub[[2]], " genes)")

  #' wl-25-07-2020, Sat: filter annotation results. Should set threshold for
  #' Percent?
  kego <- lapply(kego, function(x) {
    ## x[(x$Percent > 5) &
    ##   (x$Ontology %in% c("Biological process", "Cellular component",
    ##                      "Molecular function")), ]
    x[x$Percent > thres_anno, ]
  })

  #' wl-04-08-2020, Tue: bind together
  kego <- dplyr::bind_rows(kego, .id = "Cluster")

  ## -------------------> GO TERMS ENRICHMENT

  #' wl-04-08-2020, Tue: re-write
  universeGenes <- as.character(data_symb$Knockout)
  mat <- data_symb[idx, ]

  goen <- plyr::dlply(mat, "cluster", function(x) {
    ## x <- subset(mat, cluster == "10")
    inputGeneSet <- as.character(x$Knockout)
    ont <- c("BP", "MF", "CC") ## wl-04-08-2020, Tue: why three?
    res <- lapply(ont, function(y) {
      params <- new("GOHyperGParams",
        geneIds = inputGeneSet,
        universeGeneIds = universeGenes,
        annotation = "org.Sc.sgd.db",
        categoryName = "GO",
        ontology = y,
        pvalueCutoff = 0.05,
        conditional = T,
        testDirection = "over"
      )
      hyperGTest(params)
    })
    names(res) <- ont

    #' extract some results and move out filtering
    res_1 <- lapply(ont, function(y) {
      hgOver <- res[[y]]
      tmp <- cbind(setNames(
        tibble(
          ID = names(pvalues(hgOver)),
          Term = Term(ID),
          pvalues = pvalues(hgOver),
          oddsRatios = oddsRatios(hgOver),
          expectedCounts = expectedCounts(hgOver),
          geneCounts = geneCounts(hgOver),
          universeCounts = universeCounts(hgOver)
        ),
        c(
          "GO_ID", "Description", "Pvalue", "OddsRatio",
          "ExpCount", "Count", "CountUniverse"
        )
      ),
      Ontology = y
      ) ## %>% dplyr::filter(Pvalue <= 0.05 & Count > 1)
    })
    do.call("rbind", res_1)
  })

  names(goen) <- paste0("Cluster ", df_sub[[1]], " (", df_sub[[2]], " genes)")

  #' binding and filtering
  goen <- lapply(goen, "[", -c(4, 5, 8)) %>%
    dplyr::bind_rows(.id = "Cluster") %>%
    dplyr::filter(Pvalue <= 0.05 & Count > 1)

  ## -------------------> Output
  res <- list()
  res$stats_clusters <- df_sub # selected clusters
  res$plot_profiles <- clus_p # plot cluster profiles
  res$stats_Kegg_Goslim_annotation <- kego # KEGG AND GO SLIM ANNOTATION
  res$stats_Goterms_enrichment <- goen # GO TERMS ENRICHMENT
  return(res)
}

#' =======================================================================
#'
gene_network <- function(data = NULL, data_symb = NULL,
                        thres_clus = 10, thres_corr = 0.6) {

  ## Cluster of gene with same profile
  res_dist <- dist(data_symb[, -1], method = "manhattan")
  res_hc <- hclust(d = res_dist, method = "single")
  clus <- cutree(res_hc, h = 0) # distance 0

  data_symb$cluster <- clus

  df <- as.data.frame(table(clus), stringsAsFactors = F)
  names(df) <- c("cluster", "nGenes")
  ## filter clusters with threshold
  df_sub <- df[df$nGenes > thres_clus, ]

  #' wl-26-07-2020, Sun: remove the largest clusters?
  ## Cluster 2 (largest cluster) contains genes with no phenotype hence not
  ## considered (input to 0)
  if (T) df_sub <- df_sub[-which.max(df_sub[, 2]), ]
  #' wl-09-09-2020, Wed: if it true, the small data set may fail in network
  #  analysis

  rownames(df_sub) <- c()

  #' wl-24-07-2020, Fri: cluster index satisfing threshold of cluster number
  index <- clus %in% df_sub$cluster

  ## cluster labels with info of accumulation/decumulation of Ions
  ## (high/lower abundance)
  ## Assign label
  df_symb <- data_symb[index, ]
  lab <- plyr::ddply(df_symb, "cluster", function(x) {
    mat <- x[, !names(df_symb) %in% c("Knockout", "cluster")]
    res <- NULL
    if (length(names(which(colSums(mat == 1) > 0))) > 0) {
      res <- paste0(names(which(colSums(mat == 1) > 0)), "(+)")
    }
    if (length(names(which(colSums(mat == -1) > 0))) > 0) {
      res <- c(res, paste0(names(which(colSums(mat == -1) > 0)), "(-)"))
    }
    res <- paste(res, collapse = ", ")
  })
  names(lab)[2] <- "label"

  #' update df_sub with labels
  df_sub <- cbind(df_sub, label = lab[, 2])

  #' get cluster+gene+label names
  label <- sapply(df_symb$cluster, function(x) {
    tmp <- df_sub[df_sub$cluster == x, ]
    res <- paste0("Cluster ", tmp[1], " (", tmp[2], " genes)")
    res <- paste(res, tmp[3], sep = ": ")
  })
  df_symb$Label <- label

  ## Compute empirical correlation matrix
  corrGenes <- cor(t(as.matrix(data[, -1])),
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  ## Subset correlation matrix based on the cluster filtering
  A <- corrGenes[index, index]
  ## Diagonal value (1's) put to 0 to avoid showing edges from/to the same gene
  diag(A) <- 0

  ## Subset correlation matrix based on threshold=0.6
  ## wl-27-07-2020, Mon: need another threshold?
  A <- (A > thres_corr)
  A <- ifelse(A == TRUE, 1, 0)

  ## Generate network
  net <- network::network(A, directed = FALSE)

  #' wl-28-07-2020, Tue: add an vertex attribute and use 'Set2' in
  #'  RColorBrewer but the max. number of colors is 8 in 'Set2'
  #' wl-29-07-2020, Wed: Some layouts: circle, fruchtermanreingold,
  #'  kamadakawai, spring
  net %v% "Label" <- df_symb$Label
  tmp <- unique(df_symb$Label)
  ## fix a bug (may from ggnet2)
  if (length(tmp) != 1) {
    cpy <- rainbow(length(tmp))
    names(cpy) <- tmp
  } else {
    cpy <- "Set2"
  }
  net_p <- GGally::ggnet2(net,
    mode = "fruchtermanreingold",
    color = "Label",
    palette = cpy,
    edge.alpha = 0.5, size = 2, color.legend = "Label",
    legend.size = 10, legend.position = "right"
  )
  #' net_p

  ## Impact and betweenness
  btw <- sna::betweenness(A) # or use 'net' instead of 'A'
  impact <- apply(data[index, -1], 1, norm, type = "2") # L2 norm

  df_res <- data.frame(
    Knockout = data$Knockout[index],
    impact = round(impact, 3),
    betweenness = round(btw, 3),
    log.betweenness = round(log(btw + 1), 3),
    pos = factor(ifelse((impact < quantile(impact, .75)) &
      (log(btw + 1) < quantile(log(btw + 1), .75)), 1,
    ifelse((impact < quantile(impact, .75)) &
      (log(btw + 1) > quantile(log(btw + 1), .75)), 2,
    ifelse((impact > quantile(impact, .75)) &
      (log(btw + 1) < quantile(log(btw + 1), .75)), 3, 4)
    )
    )),
    pos.label = factor(ifelse((impact < quantile(impact, .75)) &
      (log(btw + 1) < quantile(log(btw + 1), .75)),
    "Low impact, low betweenness",
    ifelse((impact < quantile(impact, .75)) &
      (log(btw + 1) > quantile(log(btw + 1), .75)),
    "Low impact, high betweenness",
    ifelse((impact > quantile(impact, .75)) &
      (log(btw + 1) < quantile(log(btw + 1), .75)),
    "High impact, low betweenness",
    "High impact, high betweenness"
    )
    )
    ))
  )
  rownames(df_res) <- data$Knockout[index]

  q1 <- row.names(subset(df_res, (impact < quantile(impact, .75)) &
    (log.betweenness < quantile(log.betweenness, .75))))
  q2 <- row.names(subset(df_res, (impact < quantile(impact, .75)) &
    (log.betweenness > quantile(log.betweenness, .75))))
  q3 <- row.names(subset(df_res, (impact > quantile(impact, .75)) &
    (log.betweenness < quantile(log.betweenness, .75))))
  q4 <- row.names(subset(df_res, (impact > quantile(impact, .75)) &
    (log.betweenness > quantile(log.betweenness, .75))))

  #' idx <- unique(c(sample(q1,6),sample(q2,6),sample(q3,6),sample(q4,6)))
  #' wl-27-07-2020, Mon: random choose at least 24 genes to show in plot
  #' wl-17-07-2020, Fri: any or all of q1~q4 may be character(0)
  #' wl-14-07-2020, Tue: potential bug in sample replacement
  N <- 6
  lst <- list(q1, q2, q3, q4) #' sapply(lst, length)
  idx <- lapply(lst, function(x) {
    if (length(x) > N) sample(x, N) else x
  })
  idx <- unique(unlist(idx))

  df_idx <- df_res[idx, ]

  im_be_p <-
    ggplot(data = df_res, aes(x = impact, y = log.betweenness)) +
    geom_point(aes(col = pos.label), alpha = .3, size = 3) +
    scale_color_manual(values = c(
      "plum4", "palegreen4", "indianred",
      "cornflowerblue"
    )) +
    theme_linedraw() +
    theme_light() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title = element_blank()) +
    geom_text_repel(data = df_idx, aes(label = Knockout), size = 3.5) +
    geom_vline(xintercept = quantile(df_res$impact, .75), linetype = "dashed") +
    geom_hline(
      yintercept = quantile(df_res$log.betweenness, .75),
      linetype = "dashed"
    ) +
    xlab("Impact") +
    ylab("Log(betweenness+1)")

  rownames(df_res) <- c()
  df_res2 <- df_res[, -c(4, 5)]
  names(df_res2) <- c("Knockout", "Impact", "Betweenness", "Position")

  gene.cluster <- df_symb[, c("Knockout", "Label")]
  names(gene.cluster) <- c("Knockout", "Cluster")
  df_res3 <- merge(df_res2, gene.cluster, by = "Knockout", all.x = TRUE)

  #' wl-28-07-2020, Tue: better to return df_tab instead of df_tab2
  df_tab <- data.frame(table(df_res3$Cluster, df_res3$Position))
  names(df_tab) <- c("Cluster", "Position", "nGenes")
  df_tab <- dplyr::arrange(df_tab, desc(nGenes))
  df_tab2 <- df_tab %>%
    dplyr::group_by(Cluster) %>%
    top_n(1, nGenes)

  ## -------------------> Output
  res <- list()
  res$plot_pnet <- net_p # plot gene network
  res$plot_impact_betweenness <- im_be_p # plot impact betweenees
  res$stats_impact_betweenness <- df_res3 # impact betweenees data
  res$stats_impact_betweenness_by_cluster <- df_tab2 # plot position by cluster
  ## wl-28-07-2020, Tue: return this one as well
  res$stats_impact_betweenness_tab <- df_tab # contingency table
  return(res)
}
