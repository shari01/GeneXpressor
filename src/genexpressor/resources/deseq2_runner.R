suppressPackageStartupMessages({
  library(DESeq2); library(BiocParallel); library(dplyr); library(ggplot2)
  library(ggrepel); library(pheatmap); library(readr); library(tidyr)
  library(tibble); library(rmarkdown); library(RColorBrewer)
})

options(encoding = "UTF-8", useFancyQuotes = FALSE)

sanitize_ascii <- function(x) iconv(paste0(x, collapse=""), from="", to="ASCII//TRANSLIT", sub="?")
write_log <- function(path, ...) {
  raw <- paste0(...)
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", sanitize_ascii(raw), "\n")
  cat(msg); flush.console(); cat(msg, file=path, append=TRUE)
}

pkg_avail <- function(p) { suppressWarnings(requireNamespace(p, quietly = TRUE)) }
read_any <- function(path, n_max = Inf) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) return(readr::read_csv(path, show_col_types=FALSE, n_max=n_max))
  if (ext %in% c("tsv","txt")) return(readr::read_tsv(path, show_col_types=FALSE, n_max=n_max))
  if (ext %in% c("xlsx","xls")) {
    if (!pkg_avail("readxl")) stop("readxl not installed for Excel: ", path)
    return(readxl::read_excel(path, sheet=1))
  }
  if (ext %in% c("rds")) return(readRDS(path))
  if (ext %in% c("feather","parquet")) {
    if (!pkg_avail("arrow")) stop("arrow not installed for feather/parquet: ", path)
    if (ext == "feather") return(arrow::read_feather(path)) else return(arrow::read_parquet(path))
  }
  stop("Unsupported file type: ", path)
}

classify_file <- function(path) {
  tryCatch({
    x <- read_any(path, n_max = 100)
    if (!is.data.frame(x)) return(NA_character_)
    cn <- tolower(colnames(x))
    if ("condition" %in% cn) return("metadata")
    if (all(c("row","sample","count") %in% cn)) return("counts")
    if (ncol(x) >= 3) {
      first_char <- x[[1]]
      numeric_rate <- mean(vapply(x[-1], function(col) is.numeric(col) || is.integer(col), logical(1)))
      if ((is.character(first_char) || is.factor(first_char)) && numeric_rate > 0.6) return("counts")
    }
    NA_character_
  }, error=function(e) NA_character_)
}

find_dataset_files <- function(folder) {
  all_files <- list.files(folder, full.names=TRUE, recursive=FALSE)
  if (length(all_files) == 0) {
    subdirs <- list.dirs(folder, recursive=FALSE)
    if (length(subdirs) == 1) { folder <- subdirs[1]; all_files <- list.files(folder, full.names=TRUE, recursive=FALSE) }
  }
  cand <- all_files[grepl("\\.(csv|tsv|txt|xlsx|xls|rds|feather|parquet)$", tolower(all_files))]
  if (length(cand) == 0) return(list(ok=FALSE, reason="No supported files", folder=folder))

  cls <- vapply(cand, classify_file, FUN.VALUE=character(1))
  if (all(is.na(cls))) {
    meta_like  <- grepl("meta|pheno|annot|sample", tolower(basename(cand)))
    count_like <- grepl("count|expr|matrix|feature|abundance", tolower(basename(cand)))
    cls[meta_like]  <- ifelse(is.na(cls[meta_like]), "metadata", cls[meta_like])
    cls[count_like] <- ifelse(is.na(cls[count_like]), "counts", cls[count_like])
  }
  meta_files  <- cand[cls == "metadata"]
  count_files <- cand[cls == "counts"]

  if (length(meta_files) == 0) {
    for (p in cand) {
      x <- try(read_any(p, n_max=100), silent=TRUE)
      if (!inherits(x,"try-error") && is.data.frame(x)) {
        if ("condition" %in% tolower(colnames(x))) meta_files <- c(meta_files, p)
      }
    }
  }
  if (length(count_files) == 0) count_files <- setdiff(cand, meta_files)
  if (length(meta_files) == 0) return(list(ok=FALSE, reason="No metadata (needs 'condition' col)", folder=folder))
  if (length(count_files) == 0) return(list(ok=FALSE, reason="No counts file", folder=folder))

  pick_pref <- function(files, keys) {
    if (length(files) == 1) return(files[1])
    score <- vapply(basename(files), function(f) sum(grepl(keys, tolower(f))), numeric(1))
    files[order(-score)][1]
  }
  meta_path  <- pick_pref(meta_files,  "meta|pheno|annot|sample")
  count_path <- pick_pref(count_files, "count|expr|matrix|feature|abundance")
  list(ok=TRUE, metadata=meta_path, counts=count_path, folder=folder, all=cand)
}

assert_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"condition" %in% colnames(metadata)) stop("metadata must contain 'condition'")
  if (!"sample" %in% tolower(colnames(metadata))) {
    if (!is.null(rownames(metadata))) metadata$sample <- rownames(metadata)
    else stop("metadata needs 'sample' column")
  }
  names(metadata) <- make.names(names(metadata))
  if (!"sample" %in% names(metadata)) {
    sc <- names(metadata)[tolower(names(metadata)) == "sample"]
    if (length(sc)==1) names(metadata)[names(metadata)==sc] <- "sample"
  }
  if (!"condition" %in% names(metadata)) {
    cc <- names(metadata)[tolower(names(metadata)) == "condition"]
    if (length(cc)==1) names(metadata)[names(metadata)==cc] <- "condition"
  }
  metadata$sample    <- as.character(metadata$sample)
  metadata$condition <- factor(metadata$condition)
  metadata
}

clean_names <- function(v) {
  v <- as.character(v); v <- trimws(v); v <- gsub("\\s+", "_", v)
  v <- gsub("\\.+", "_", v); v <- gsub("-", "_", v); v
}

to_deseq <- function(counts, metadata) {
  counts   <- as.data.frame(counts, stringsAsFactors = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  metadata$sample    <- as.character(metadata$sample)
  metadata$condition <- factor(metadata$condition)

  lc <- tolower(colnames(counts))
  tidy_cols <- c("row","sample","count")
  if (all(tidy_cols %in% lc)) {
    colnames(counts) <- lc
    counts$sample <- clean_names(counts$sample)
    metadata$sample <- clean_names(metadata$sample)
    DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition, tidy=TRUE)
  } else {
    gene_col <- colnames(counts)[1]
    mat <- as.data.frame(counts, stringsAsFactors = FALSE)
    rownames(mat) <- as.character(mat[[gene_col]]); mat[[gene_col]] <- NULL
    num_ok <- vapply(mat, function(x) is.numeric(x) || is.integer(x), logical(1))
    if (!any(num_ok)) stop("No numeric columns found in counts")
    mat <- cbind.data.frame(mat[, num_ok, drop=FALSE]); mat[] <- lapply(mat, function(x) as.integer(round(x,0)))
    colnames(mat) <- clean_names(colnames(mat)); metadata$sample <- clean_names(metadata$sample)

    common <- intersect(colnames(mat), metadata[["sample"]])
    if (length(common) < 2) {
      common <- intersect(tolower(colnames(mat)), tolower(metadata[["sample"]]))
      if (length(common) >= 2) {
        idx <- match(tolower(colnames(mat)), common, nomatch=0) > 0
        colnames(mat)[idx] <- tolower(colnames(mat)[idx]); metadata[["sample"]] <- tolower(metadata[["sample"]])
      }
    }
    common <- intersect(colnames(mat), metadata[["sample"]])
    if (length(common) < 2) stop("No sample overlap between counts columns and metadata$sample")

    mat <- mat[, common, drop=FALSE]
    rownames(metadata) <- metadata[["sample"]]
    metadata <- metadata[common, , drop=FALSE]
    metadata$condition <- factor(metadata$condition)

    DESeqDataSetFromMatrix(countData=as.matrix(mat), colData=metadata, design=~condition)
  }
}

save_plot <- function(p, path, w=8, h=6) ggsave(path, plot=p, dpi=300, width=w, height=h)

run_dataset <- function(subfolder, idx_prefix,
                        CASE_LEVEL, CONTROL_LVL, ALPHA, LFC_THR,
                        TOP_LABELS, TOP_HEATMAP, MAKE_REPORT, DEBUG) {

  pretty <- basename(subfolder)
  found <- find_dataset_files(subfolder)
  if (!isTRUE(found$ok)) { cat(sprintf("[X] [%02d] %s  ->  SKIPPED (%s)\n", idx_prefix, pretty, found$reason)); return(invisible(FALSE)) }

  metadata   <- try(read_any(found$metadata), silent=TRUE)
  counts_obj <- try(read_any(found$counts),   silent=TRUE)
  if (inherits(metadata,"try-error") || !is.data.frame(metadata)) { cat(sprintf("[X] [%02d] %s  ->  SKIPPED (cannot read metadata)\n", idx_prefix, pretty)); return(invisible(FALSE)) }
  metadata <- try(assert_metadata(metadata), silent=TRUE)
  if (inherits(metadata,"try-error")) { cat(sprintf("[X] [%02d] %s  ->  SKIPPED (bad metadata)\n", idx_prefix, pretty)); return(invisible(FALSE)) }

  outdir     <- file.path(subfolder, sprintf("%02d_%s_DESeq2_results", idx_prefix, pretty))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  dir_tables <- file.path(outdir, "tables"); dir_plots <- file.path(outdir, "plots")
  dir.create(dir_tables, showWarnings=FALSE, recursive=TRUE)
  dir.create(dir_plots,  showWarnings=FALSE, recursive=TRUE)
  log_path   <- file.path(outdir, "run_log.txt")

  if (DEBUG) {
    dbg <- file.path(subfolder, sprintf("%02d_%s_debug_summary.txt", idx_prefix, pretty))
    cat(file=dbg, "FOLDER: ", normalizePath(subfolder), "\n",
        "METADATA: ", found$metadata, "\n",
        "COUNTS: ", found$counts, "\n",
        "ALL_CANDIDATES:\n  - ", paste(found$all, collapse="\n  - "), "\n", sep="")
  }
  write_log(log_path, "[", sprintf("%02d", idx_prefix), "] Start dataset: ", pretty)

  dds <- NULL
  if (inherits(counts_obj, "DESeqDataSet")) {
    dds <- counts_obj
    if (!"condition" %in% colnames(colData(dds))) { write_log(log_path, "DESeqDataSet lacks 'condition' -> skip"); return(invisible(FALSE)) }
    design(dds) <- ~ condition
  } else if (inherits(counts_obj, "SummarizedExperiment")) {
    if (!"condition" %in% colnames(colData(counts_obj))) { write_log(log_path, "SummarizedExperiment lacks 'condition' -> skip"); return(invisible(FALSE)) }
    colData(counts_obj)$condition <- factor(as.character(colData(counts_obj)$condition))
    dds <- DESeqDataSet(counts_obj, design = ~ condition)
  } else if (is.data.frame(counts_obj)) {
    dds <- try(to_deseq(counts_obj, metadata), silent=TRUE)
    if (inherits(dds,"try-error")) { write_log(log_path, "Failed to construct DESeqDataSet"); return(invisible(FALSE)) }
  } else { write_log(log_path, "Unknown counts object -> skip"); return(invisible(FALSE)) }

  if (!all(c(CASE_LEVEL, CONTROL_LVL) %in% unique(colData(dds)$condition))) { write_log(log_path, "Condition factor missing case/control -> skip"); return(invisible(FALSE)) }
  colData(dds)$condition <- factor(colData(dds)$condition)
  if (levels(colData(dds)$condition)[1] != CONTROL_LVL) colData(dds)$condition <- stats::relevel(colData(dds)$condition, ref=CONTROL_LVL)

  write_log(log_path, "Running DESeq2 (estimateSizeFactors -> DESeq)...")
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds, parallel = (.Platform$OS.type == "unix"))

  write_log(log_path, "Extracting results & LFC shrinkage (apeglm)...")
  res <- results(dds, contrast = c("condition", CASE_LEVEL, CONTROL_LVL), alpha = ALPHA)
  resLFC <- tryCatch(lfcShrink(dds, contrast=c("condition", CASE_LEVEL, CONTROL_LVL), type="apeglm"),
                     error=function(e){ message("lfcShrink failed; using unshrunken results"); res })

  res_df <- as.data.frame(resLFC) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      is_significant = !is.na(padj) & padj < ALPHA & !is.na(log2FoldChange) & abs(log2FoldChange) >= LFC_THR,
      regulation     = dplyr::case_when(
        is_significant & log2FoldChange >  0 ~ "upregulated",
        is_significant & log2FoldChange <  0 ~ "downregulated",
        TRUE ~ "not_significant"
      )
    )

  total_genes <- nrow(res_df)
  sig_only    <- dplyr::filter(res_df, is_significant)
  up_tbl      <- dplyr::filter(sig_only, regulation == "upregulated")   %>% dplyr::arrange(dplyr::desc(log2FoldChange))
  down_tbl    <- dplyr::filter(sig_only, regulation == "downregulated") %>% dplyr::arrange(log2FoldChange)

  readr::write_csv(res_df,  file.path(dir_tables, "DESeq2_full_results.csv"))
  readr::write_csv(sig_only, file.path(dir_tables, "DESeq2_significant.csv"))
  if (nrow(up_tbl)   > 0) readr::write_csv(up_tbl,   file.path(dir_tables, "DESeq2_upregulated.csv"))
  if (nrow(down_tbl) > 0) readr::write_csv(down_tbl, file.path(dir_tables, "DESeq2_downregulated.csv"))

  write_log(log_path, "Total genes: ", total_genes,
            " | Significant: ", nrow(sig_only),
            " | Up: ", nrow(up_tbl),
            " | Down: ", nrow(down_tbl))
  if (nrow(up_tbl) > 0)  write_log(log_path, "Top upregulated: ", paste(head(up_tbl$gene, 5), collapse=", "))
  if (nrow(down_tbl) > 0)write_log(log_path, "Top downregulated: ", paste(head(down_tbl$gene, 5), collapse=", "))

  write_log(log_path, "Making PCA...")
  vsd  <- vst(dds, blind = FALSE)
  pdat <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  varp <- attr(pdat, "percentVar") * 100
  p_pca <- ggplot(pdat, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(title="DESeq2 PCA",
         x=paste0("PC1 (", round(varp[1],2), "%)"),
         y=paste0("PC2 (", round(varp[2],2), "%)")) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(dir_plots, "PCA.png"), plot = p_pca, dpi = 300, width = 8, height = 6)

  write_log(log_path, "Making Volcano...")
  thr_line <- -log10(ALPHA)
  res_df$regulation <- factor(res_df$regulation, levels=c("upregulated","downregulated","not_significant"))
  p_vol <- ggplot(dplyr::filter(res_df, !is.na(padj)),
                  aes(log2FoldChange, -log10(padj), color = regulation)) +
    geom_point(alpha = 0.8, size = 1.6) +
    scale_color_manual(values = c(upregulated="#D73027", downregulated="#1E90FF", not_significant="grey75")) +
    geom_hline(yintercept = thr_line, linetype = "dashed", color="grey40") +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color="grey40") +
    labs(title="DESeq2 Volcano Plot", x="log2(Fold Change)", y="-log10(FDR)", color="Regulation") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(dir_plots, "Volcano.png"), plot = p_vol, dpi = 300, width = 8, height = 6)

  top_lbl <- sig_only %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = TOP_LABELS)
  if (nrow(top_lbl) > 0) {
    p_vol_lbl <- p_vol + ggrepel::geom_text_repel(data = top_lbl, aes(label = gene), size = 3, max.overlaps = 50)
    ggsave(file.path(dir_plots, "Volcano_labeled.png"), plot = p_vol_lbl, dpi = 300, width = 10, height = 7)
  }

  write_log(log_path, "Making sample-distance heatmap...")
  d <- stats::dist(t(assay(vsd))); m <- as.matrix(d)
  png(file.path(dir_plots,"SampleDistanceHeatmap.png"), width=1400, height=1100, res=150)
  pheatmap::pheatmap(m, clustering_distance_rows=d, clustering_distance_cols=d,
                     main="Sample distance", color = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255))
  dev.off()

  write_log(log_path, "Making top DEG heatmap...")
  sig_genes <- sig_only %>% dplyr::arrange(padj) %>% dplyr::pull(gene)
  tg <- head(sig_genes, TOP_HEATMAP)
  if (length(tg) >= 2) {
    hm <- assay(vsd)[tg, , drop = FALSE]
    hm <- t(scale(t(hm)))
    ann <- as.data.frame(SummarizedExperiment::colData(vsd)[, "condition", drop = FALSE])
    png(file.path(dir_plots,"TopDEG_Heatmap.png"), width=1400, height=1400, res=150)
    pheatmap::pheatmap(hm, annotation_col=ann, show_rownames=TRUE, main="Top DE genes (z-scored VST)",
                       color = colorRampPalette(brewer.pal(9,"PuOr"))(255))
    dev.off()
  }

  write_log(log_path, "Making Top30 bar charts...")
  topN <- 30L
  sig_only <- sig_only
  top_up   <- sig_only %>% dplyr::filter(regulation=="upregulated")   %>% dplyr::slice_head(n=topN)
  top_down <- sig_only %>% dplyr::filter(regulation=="downregulated") %>% dplyr::slice_head(n=topN)
  top_abs  <- sig_only %>% dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>% dplyr::slice_head(n=topN)

  if (nrow(top_up) > 0) {
    p_up <- ggplot(top_up, aes(x=reorder(gene, log2FoldChange), y=log2FoldChange, fill=log2FoldChange)) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low="#FFD92F", high="#D73027") +
      labs(title=paste0("Top ", topN, " Upregulated Genes (by log2FC)"), x="Gene", y="log2FC", fill="log2FC") +
      theme_minimal() + theme(panel.grid.minor = element_blank())
    ggsave(file.path(dir_plots, "Top30_Upregulated.png"), plot = p_up, dpi = 300, width = 9, height = 7)
  }
  if (nrow(top_down) > 0) {
    p_down <- ggplot(top_down, aes(x=reorder(gene, log2FoldChange), y=log2FoldChange, fill=log2FoldChange)) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low="#1E90FF", high="#002147") +
      labs(title=paste0("Top ", topN, " Downregulated Genes (by log2FC)"), x="Gene", y="log2FC", fill="log2FC") +
      theme_minimal() + theme(panel.grid.minor = element_blank())
    ggsave(file.path(dir_plots, "Top30_Downregulated.png"), plot = p_down, dpi = 300, width = 9, height = 7)
  }
  if (nrow(top_abs) > 0) {
    p_abs <- ggplot(top_abs, aes(x=reorder(gene, abs(log2FoldChange)), y=abs(log2FoldChange),
                                 fill=abs(log2FoldChange))) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low="#B2DF8A", high="#33A02C") +
      labs(title=paste0("Top ", topN, " Genes by |log2FC|"), x="Gene", y="|log2FC|", fill="|log2FC|") +
      theme_minimal() + theme(panel.grid.minor = element_blank())
    ggsave(file.path(dir_plots, "Top30_AbsLFC.png"), plot = p_abs, dpi = 300, width = 9, height = 7)
  }

  if (MAKE_REPORT) {
    write_log(log_path, "Rendering HTML report...")
    report_rmd <- file.path(outdir, "report.Rmd")
    writeLines(c(
      "---", 'title: "DESeq2 Differential Expression Report"', "output: html_document", "---", "",
      paste0("**Dataset:** ", sprintf("%02d", idx_prefix), " - ", pretty),
      paste0("**Contrast:** ", CASE_LEVEL, " vs ", CONTROL_LVL),
      paste0("**FDR (alpha):** ", ALPHA, " | **|log2FC| threshold:** ", LFC_THR), "",
      "## Summary",
      paste0("- Total genes: ", total_genes),
      paste0("- Significant: ", nrow(sig_only), " (Up: ", nrow(top_up), ", Down: ", nrow(top_down), ")"), "",
      "## Key Plots",
      "![](plots/PCA.png)", "", "![](plots/Volcano.png)", "", "![](plots/Volcano_labeled.png)", "",
      "![](plots/SampleDistanceHeatmap.png)", "", "![](plots/TopDEG_Heatmap.png)", "",
      "![](plots/Top30_Upregulated.png)", "", "![](plots/Top30_Downregulated.png)", "", "![](plots/Top30_AbsLFC.png)", "",
      "## Result Tables",
      "- `tables/DESeq2_full_results.csv`",
      "- `tables/DESeq2_significant.csv`",
      "- `tables/DESeq2_upregulated.csv`",
      "- `tables/DESeq2_downregulated.csv`", "",
      "## Session Info", "```{r}", "sessionInfo()", "```"
    ), con = report_rmd)
    rmarkdown::render(report_rmd, output_file=file.path(outdir, "DESeq2_report.html"), quiet=TRUE)
  }

  write_log(log_path, "Done. Results at: ", outdir)
  cat(sprintf("[OK] [%02d] %s  ->  %s\n", idx_prefix, pretty, outdir))
  invisible(TRUE)
}

run_all <- function(PARENT_DIR,
                    PICK="AUTO",
                    CASE_LEVEL="Disease",
                    CONTROL_LVL="Control",
                    ALPHA=0.05,
                    LFC_THR=1.0,
                    TOP_LABELS=20,
                    TOP_HEATMAP=50,
                    MAKE_REPORT=TRUE,
                    DEBUG=TRUE,
                    threads=1,
                    COUNTS_OVERRIDE=NULL,
                    METADATA_OVERRIDE=NULL,
                    DATASET_INDEX_OVERRIDE=NA_integer_) {

  stopifnot(dir.exists(PARENT_DIR))
  if (!is.na(DATASET_INDEX_OVERRIDE)) PICK <- as.character(DATASET_INDEX_OVERRIDE)

  if (!is.null(COUNTS_OVERRIDE) || !is.null(METADATA_OVERRIDE)) {
    if (is.null(COUNTS_OVERRIDE) || is.null(METADATA_OVERRIDE)) stop("Both --counts and --metadata must be provided when overriding.")
    tmp <- file.path(PARENT_DIR, sprintf("override_%s", as.integer(Sys.time())))
    dir.create(tmp, showWarnings=FALSE, recursive=TRUE)
    file.copy(COUNTS_OVERRIDE,   file.path(tmp, basename(COUNTS_OVERRIDE)),   overwrite=TRUE)
    file.copy(METADATA_OVERRIDE, file.path(tmp, basename(METADATA_OVERRIDE)), overwrite=TRUE)
    if (.Platform$OS.type == "unix") { try(BiocParallel::register(BiocParallel::MulticoreParam(workers = threads)), silent=TRUE)
    } else { BiocParallel::register(BiocParallel::SerialParam()) }
    ok <- try(run_dataset(tmp, 1L, CASE_LEVEL, CONTROL_LVL, ALPHA, LFC_THR, TOP_LABELS, TOP_HEATMAP, MAKE_REPORT, DEBUG), silent=TRUE)
    BiocParallel::register(BiocParallel::SerialParam())
    return(isTRUE(ok))
  }

  root_candidates <- list.files(PARENT_DIR, full.names=TRUE, recursive=FALSE)
  root_has_files <- any(grepl("\\.(csv|tsv|txt|xlsx|xls|rds|feather|parquet)$", tolower(root_candidates)))
  if (root_has_files) subfolders <- c(PARENT_DIR) else subfolders <- list.dirs(PARENT_DIR, recursive=FALSE)
  if (length(subfolders) == 0) stop("No dataset folders/files found in: ", PARENT_DIR)

  sf_tbl <- tibble::tibble(idx = seq_along(subfolders), path = subfolders, name = basename(subfolders))

  if (toupper(PICK) == "AUTO") {
    valid_idx <- integer()
    for (i in sf_tbl$idx) {
      found <- find_dataset_files(sf_tbl$path[sf_tbl$idx==i])
      if (isTRUE(found$ok)) valid_idx <- c(valid_idx, i)
    }
    chosen_idx <- valid_idx
    if (length(chosen_idx) == 1) {
      message("AUTO selected dataset: [", chosen_idx, "] ", basename(sf_tbl$path[sf_tbl$idx==chosen_idx]))
    } else if (length(chosen_idx) == 0) {
      message("AUTO found no valid datasets (check debug files in subfolders).")
    }
  } else if (toupper(PICK) == "ALL") {
    chosen_idx <- sf_tbl$idx
  } else {
    chosen_idx <- suppressWarnings(as.integer(unlist(strsplit(PICK, "\\s*,\\s*"))))
    chosen_idx <- chosen_idx[!is.na(chosen_idx)]
    if (length(chosen_idx) == 0 || any(!chosen_idx %in% sf_tbl$idx)) stop("Invalid PICK")
  }

  if (.Platform$OS.type == "unix") {
    try(BiocParallel::register(BiocParallel::MulticoreParam(workers = threads)), silent=TRUE)
  } else {
    BiocParallel::register(BiocParallel::SerialParam())
  }

  ran_ok <- logical()
  for (i in chosen_idx) {
    sf <- sf_tbl$path[sf_tbl$idx == i]
    ok <- try(run_dataset(sf, idx_prefix = i,
                          CASE_LEVEL=CASE_LEVEL, CONTROL_LVL=CONTROL_LVL,
                          ALPHA=ALPHA, LFC_THR=LFC_THR,
                          TOP_LABELS=TOP_LABELS, TOP_HEATMAP=TOP_HEATMAP,
                          MAKE_REPORT=MAKE_REPORT, DEBUG=DEBUG), silent=TRUE)
    ran_ok <- c(ran_ok, isTRUE(ok))
  }
  BiocParallel::register(BiocParallel::SerialParam())
  any(ran_ok)
}
