# execution

source(file = "") # point towards 1_functions.R
working_dir <- "" # point to the directory where you want your output to be there (write absolute path)

is_bipartite=FALSE # TRUE/FALSE
dir.create(path = working_dir,recursive = T)

setwd(working_dir)

library(dplyr)
library(purrr)
library(readr)

configs<-read_csv("") # point to the directory where your config.txt is there

run_pipeline_one <- function(cell_line,
                             condition,
                             loops_file,
                             prom_file,
                             se_file,
                             te_file,
                             out_dir,
                             relax_by = 10000) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("====================================\n")
  cat("Running pipeline for: ", cell_line, " | ", condition,"\n")
  cat("====================================\n")
  
  relax_by=relax_by
  cat("relaxing loops by ", relax_by)
  
  # 1. process input files
  a <- process_loops_se_te_pr(
    loops_file = loops_file,
    prom_file  = prom_file,
    se_file    = se_file,
    te_file    = te_file,
    relax_by   = relax_by
  )
  
  se    <- a$se
  te    <- a$te
  pr    <- a$pr
  loops <- a$loops
  
  # 2. find overlaps / truth table
  b <- give_truth_table(loops, te, se, pr)
  
  pr_gr    <- b$promoter_gr
  se_gr    <- b$super_enhancers_gr
  te_gr    <- b$typical_enhancers_gr
  loops_gi <- b$loops_gi
  
  hits1_p  <- b$hits1_p
  hits2_p  <- b$hits2_p
  hits1_se <- b$hits1_se
  hits2_se <- b$hits2_se
  hits1_te <- b$hits1_te
  hits2_te <- b$hits2_te
  
  # 3. attach loop IDs
  loops_gi <- giveloop_ids(
    loops_gi, pr_gr, te_gr, se_gr,
    hits1_p, hits2_p,
    hits1_se, hits2_se,
    hits1_te, hits2_te
  )
  
  # 4. classify anchor interaction categories
  loops_gi <- make_loop_interaction_categories(loops_gi)
  
  # 5. generate edges and nodes
  e <- give_total_edges_nodes(loops_g1 = loops_gi, loops = loops)
  nodes <- e$total_nodes
  edges <- e$total_edges
  
  # 6. construct network
  prefix <- file.path(out_dir, paste0(condition, "_"))
  f <- construct_nw(total_nodes = nodes, total_edges = edges, path=prefix, bipartite = is_bipartite)
  
  nodes_composition <- f$df
  cluster_type_count <- f$df1
  cluster_sizes <- f$df2
  
  cluster_classification <- f$cluster_classification
  g <- f$graph
  components_g <- f$components_g
  
  # 7. calculate network properties
  prefix <- file.path(out_dir, paste0(condition, "_"))
  enhancer_ppts <- calculate_nw_ppts(
    g,
    components_g,
    cluster_classification,
    path = prefix
  )
  
  # 8. attach cluster types
  enhancer_ppts_2 <- enhancer_ppts %>%
    left_join(
      cluster_classification %>% select(Cluster, cluster_type),
      by = "Cluster"
    ) %>%
    mutate(
      cell_line = cell_line,
      condition = condition
    )
  
  # 9. determine metrics
  metrics <- enhancer_ppts_2 %>%
    select(-any_of(c("name", "LA1_type", "Cluster", "cluster_type", "cell_line", "condition"))) %>%
    select(where(is.numeric)) %>%
    names()
  
  # 10. save tabular outputs
  write.csv(enhancer_ppts_2,
            file = file.path(out_dir, paste0("enhancer_properties_", condition, ".csv")),
            row.names = FALSE)
  
  write.csv(nodes,
            file = file.path(out_dir, paste0("nodes_", condition, ".csv")),
            row.names = FALSE, quote = F)
  
  write.csv(edges,
            file = file.path(out_dir, paste0("edges_", condition, ".csv")),
            row.names = FALSE, quote = F)
  
  list(
    cell_line = cell_line,
    condition = condition,
    loops = loops,
    loops_gi = loops_gi,
    nodes = nodes,
    edges = edges,
    graph = g,
    components_g = components_g,
    cluster_classification = cluster_classification,
    enhancer_ppts = enhancer_ppts,
    enhancer_ppts_2 = enhancer_ppts_2, 
    nodes_composition = nodes_composition,
    cluster_type_count = cluster_type_count,
    cluster_sizes = cluster_sizes
  )
}


log_file <- "results/pipeline.log"

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

all_results <- purrr::pmap(configs, run_pipeline_one)

sink(type = "message")
sink()
close(log_con)

names(all_results) <- paste(configs$cell_line, configs$condition, sep = "_")

combined_enhancer_ppts <- bind_rows(
  lapply(all_results, function(x) x$enhancer_ppts_2)
)

write.csv(
  combined_enhancer_ppts,
  "results/all_enhancer_properties.csv",
  row.names = FALSE
)

saveRDS(object = all_results,file = "results/all_results.rds")

