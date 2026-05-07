# MODULE: run multiple igraph community detection methods and write xlsx for each
# 修复后的 MODULE 函数（针对你报告的问题进行修正）
MODULE <- function(edges,
                   heatmap = TRUE,
                   tarGet = FALSE,
                   write_output = TRUE,
                   outdir = ".",
                   pval_cutoff = 0.05,
                   fdr_method = "BH") {
  
  # ----- deps -----
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Please install package 'igraph' before running MODULE().")
  # xlsx output support
  use_openxlsx <- FALSE
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    use_openxlsx <- TRUE
  } else if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Please install 'openxlsx' or 'writexl' to enable xlsx output.")
  }
  
  library(igraph)
  
  # ----- input checks -----
  if (!is.data.frame(edges)) stop("Param 'edges' must be a data.frame with at least two columns.")
  if (ncol(edges) < 2) stop("Param 'edges' must have at least two columns (interactorA, interactorB).")
  edges_df <- as.data.frame(edges[, 1:2], stringsAsFactors = FALSE)
  colnames(edges_df)[1:2] <- c("V1", "V2")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  time_start_all <- Sys.time()
  
  # ----- build graph and basic clean -----
  g <- graph_from_data_frame(edges_df, directed = FALSE)
  # ensure simple graph (remove loops/multiple) and ensure vertex names exist
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
  
  node_names <- V(g)$name
  
  methods <- c("louvain", "fast_greedy", "leiden", "infomap", "walktrap", "leading_eigen")
  results <- vector("list", length(methods))
  names(results) <- methods
  timings <- vector("list", length(methods))
  names(timings) <- methods
  
  for (m in methods) {
    message(sprintf("Starting community detection: %s ...", m))
    t0 <- Sys.time()
    df_out <- NULL
    
    tryCatch({
      if (m == "louvain") {
        com <- igraph::cluster_louvain(g)
      } else if (m == "fast_greedy") {
        com <- igraph::cluster_fast_greedy(g)
      } else if (m == "leiden") {
        if (exists("cluster_leiden", where = asNamespace("igraph"), mode = "function")) {
          com <- igraph::cluster_leiden(g)
        } else {
          stop("cluster_leiden is not available in your igraph namespace.")
        }
      } else if (m == "infomap") {
        com <- igraph::cluster_infomap(g)
      } else if (m == "walktrap") {
        com <- igraph::cluster_walktrap(g)
      } else if (m == "leading_eigen") {
        # leading_eigen may fail (ARPACK). We'll allow it to error and be caught.
        com <- igraph::cluster_leading_eigen(g)
      } else {
        stop("Unknown method")
      }
      
      # membership vector (ordered by vertex sequence)
      mem <- membership(com)
      # mem should be length vcount(g); if for some reason lengths mismatch, handle gracefully
      if (length(mem) != vcount(g)) {
        # try to coerce: create full-length vector using vertex sequence order
        # but safest is to error out
        stop(sprintf("membership length (%d) does not match number of vertices (%d)", length(mem), vcount(g)))
      }
      
      # build data.frame using V(g)$name to ensure consistent names
      df_out <- data.frame(name = as.character(V(g)$name),
                           module = as.integer(mem),
                           stringsAsFactors = FALSE)
      
      # sanity: verify lengths match
      if (nrow(df_out) != vcount(g)) stop("Internal error: output rows != number of vertices")
      
      # write excel if requested
      if (write_output && !is.null(df_out)) {
        fname <- file.path(outdir, paste0("modules_", m, ".xlsx"))
        if (use_openxlsx) {
          openxlsx::write.xlsx(df_out, file = fname, rowNames = FALSE)
        } else {
          writexl::write_xlsx(df_out, path = fname)
        }
        message(sprintf("Wrote %s (n = %d) to %s", m, nrow(df_out), fname))
      }
      
      results[[m]] <- df_out
      t1 <- Sys.time()
      timings[[m]] <- as.numeric(difftime(t1, t0, units = "secs"))
      message(sprintf("%s finished in %.2f sec.", m, timings[[m]]))
    }, error = function(e) {
      # capture error, mark result null and timing NA
      message(sprintf("Method %s failed or not available: %s", m, conditionMessage(e)))
      results[[m]] <- NULL
      timings[[m]] <- NA_real_
    })
  }
  
  time_end_all <- Sys.time()
  total_secs <- as.numeric(difftime(time_end_all, time_start_all, units = "secs"))
  message(sprintf("All methods done. Total elapsed: %.2f sec.", total_secs))
  
  return(list(results = results, timings = timings, total_time_sec = total_secs))
}



library(readxl)
edges_core <- read_xlsx("top10%_CoreNet.xlsx")
edges_core <- as.data.frame(edges_core)

res <- MODULE(edges_core,
              heatmap = TRUE,
              tarGet = FALSE,
              write_output = TRUE,
              outdir = "modules_output",
              pval_cutoff = 0.05,
              fdr_method = "BH")

# 查看结果
names(res$results)       # 各方法名
res$timings              # 每个方法的耗时（秒）