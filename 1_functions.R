# functions

# This contains all the functions required to build network Hi-C/HiChIP loops data

# Requirements: 
# Loops file: minimum with the coordinates, optionally value, (columns: chrom1, start1, end1, chrom2, start2, end2, and optionally, value)
# Promoters, Super Enhancer and Typical Enhancer coordinates pertaining to the cell line/ patient of interest. (Format: tsv, column names: chromosome, start, end, id)

# Abbreviations used: 
# Super Enhancer (SE)
# Typical Enhancer (TE)
# Promoter (P)

# download packages 
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(tibble)
library(purrr)  
library(igraph)
library(tidyverse)



# ---- 1. Process each of the files and create dataframes ----

# Function name: process_loops_se_te_pr()

# Input:
# loops_file -> path to the loops file (as csv)
# prom_file -> path to the promoter coordinates file (as csv)
# se_file -> path to the super enhancer coordinates file (as csv)
# te_file -> path to the typical enhancer coordinates file (as csv)
# relax_by -> number of bases to relax the loop coordinates by, for better instersection (in bp)

# What it does: 

# Output: 
# se, te and pr <- dataframe of super enhancers, typical enhancers and promoters respectively, having consistent column names and formatted in order
# loops <- dataframe of loops

process_loops_se_te_pr <- function(loops_file, prom_file, se_file, te_file,
                                   relax_by){
  
  # ---- process the loops file -----
  loops=as.data.frame(read.csv(loops_file))
  
  if (ncol(loops) < 6) {
    stop("Loops file has fewer than 6 columns.")
  }
  else {
    loops <- loops[, c(1:6, 7)]
    colnames(loops) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "value")
  }
  
  loops <- loops %>% mutate(loop_length=abs(start2-start1))
  
  # relax the loops if possible
  if(loops %>% filter(loop_length<relax_by) %>% nrow()==0){
    
    loops <- loops %>% 
      mutate(start1 = start1 - relax_by,
             end1   = end1 + relax_by,
             start2 = start2 - relax_by,
             end2   = end2 + relax_by)
    
    cat("\n the min looplength is : ", loops %>% select(loop_length) %>% min())
    
    cat("\n the loops have been relaxed by : ", relax_by)
    
  } else{
    
    cat("\n loops cant be relaxed")
    
  }
  
  # ---- make the SE, TE and promoter files in format wanted ----
  
  se <- as.data.frame(read.csv(se_file))
  
  colnames(se) <- c("chr", "start", "end", "id")
  
  te <- as.data.frame(read.csv(te_file))
  
  colnames(te) <- c("chr", "start", "end", "id")
  
  pr <- as.data.frame(read.csv(prom_file))
  
  colnames(pr) <- c("chr", "start", "end", "id")
  
  return(list(se=se,
              te=te,
              pr=pr,
              loops=loops))
}


# ---- 2. Find hits of the promoters, SE and TE with the loop anchors  ----

# Function name: give_truth_table()

# Input:
# loops, te, se, pr <- dataframe of loops, typical enhancers, super enhancers and promoters, respectively, created using process_loops_se_te_pr() function

# What it does: 

# Output:
# loops_gi <- GInteraction object of loop coordinates
# promoter_gr, super_enhancers_gr, typical_enhancers_gr <- GRanges object of P, SE and TE coordinates respectively
# hits1_p, hits2_p <- S4Vector object containing the hits between P and anchor1 and anchor2 of loops respectively
# hits1_se, hits2_se <- S4Vector object containing the hits between SE and anchor1 and anchor2 of loops respectively
# hits1_te, hits2_te <- S4Vector object containing the hits between TE and anchor1 and anchor2 of loops respectively

give_truth_table<-function(loops, te, se, pr){
  
  # library(GenomicRanges)
  
  # ---- Make a minimal copy of the loops depending on how many columns are there ----
  
  loops_min <- loops[, c("chr1","start1","end1","chr2","start2","end2","value")]
  
  loops_min$chr1   <- as.character(loops_min$chr1)
  loops_min$chr2   <- as.character(loops_min$chr2)
  loops_min$start1 <- as.integer(loops_min$start1)
  loops_min$end1   <- as.integer(loops_min$end1)
  loops_min$start2 <- as.integer(loops_min$start2)
  loops_min$end2   <- as.integer(loops_min$end2)
  
  n <- nrow(loops_min)
  
  cat("\n total loops =", n)
  
  # ---- Make GInteraction object from the loops ---- 
  
  anchor1 <- GenomicRanges::GRanges(
    seqnames = loops_min$chr1,
    ranges   = IRanges(loops_min$start1, loops_min$end1)
  )
  
  anchor2 <- GenomicRanges::GRanges(
    seqnames = loops_min$chr2,
    ranges   = IRanges(loops_min$start2, loops_min$end2)
  )
  
  loops_gi <- InteractionSet::GInteractions(anchor1, anchor2)
  
  mcols(loops_gi)$value <- loops_min$value
  
  
  
  # ---- Make the GRanges object of the P, SE and TE dfs ----
  
  promoter_gr<-GRanges(seqnames = pr$chr,ranges = IRanges(start=pr$start, end=pr$end))
  mcols(promoter_gr)$id<-pr$id
  
  typical_enhancers_gr<-GRanges(seqnames = te$chr,ranges = IRanges(start=te$start, end=te$end))
  mcols(typical_enhancers_gr)$id<-te$id
  
  super_enhancers_gr<-GRanges(seqnames = se$chr,ranges = IRanges(start=se$start, end=se$end))
  mcols(super_enhancers_gr)$id<-se$id
  
  # ---- Find overlaps between loop anchors and genomic elements (SE, TE and P) ----
  
  a1<-InteractionSet::anchors(loops_gi, type="first")
  a2<-InteractionSet::anchors(loops_gi, type="second")
  
  hits1_se<-findOverlaps(a1,super_enhancers_gr)
  hits2_se<-findOverlaps(a2,super_enhancers_gr)
  
  hits1_te<-findOverlaps(a1,typical_enhancers_gr)
  hits2_te<-findOverlaps(a2,typical_enhancers_gr)
  
  hits1_p<-findOverlaps(a1,promoter_gr)
  hits2_p<-findOverlaps(a2,promoter_gr)
  
  N <- length(loops_gi)
  
  is_P_a1  <- rep(FALSE, N)
  is_P_a2  <- rep(FALSE, N)
  is_TE_a1 <- rep(FALSE, N)
  is_TE_a2 <- rep(FALSE, N)
  is_SE_a1 <- rep(FALSE, N)
  is_SE_a2 <- rep(FALSE, N)
  
  is_P_a1[ unique(queryHits(hits1_p)) ]  <- TRUE
  is_P_a2[ unique(queryHits(hits2_p)) ]  <- TRUE
  is_TE_a1[ unique(queryHits(hits1_te)) ] <- TRUE
  is_TE_a2[ unique(queryHits(hits2_te)) ] <- TRUE
  is_SE_a1[ unique(queryHits(hits1_se)) ] <- TRUE
  is_SE_a2[ unique(queryHits(hits2_se)) ] <- TRUE
  
  mcols(loops_gi)$is_P_a1  <- is_P_a1
  mcols(loops_gi)$is_P_a2  <- is_P_a2
  mcols(loops_gi)$is_TE_a1 <- is_TE_a1
  mcols(loops_gi)$is_TE_a2 <- is_TE_a2
  mcols(loops_gi)$is_SE_a1 <- is_SE_a1
  mcols(loops_gi)$is_SE_a2 <- is_SE_a2
  
  return(list(loops_gi=loops_gi,
              promoter_gr=promoter_gr,
              super_enhancers_gr=super_enhancers_gr, 
              typical_enhancers_gr=typical_enhancers_gr,
              hits1_p=hits1_p,
              hits2_p=hits2_p,
              hits1_se=hits1_se,
              hits2_se=hits2_se,
              hits1_te=hits1_te,
              hits2_te=hits2_te))
}



# ---- 3. Extracting the ids of the hits ----

# Function name: giveloop_ids()

# Input:
# loops_gi <- GInteraction object of loop coordinates
# pr, te, se <- granges object of super enhancers, typical enhancers and promoters respectively, having consistent column names and formatted in order (from process_loops_se_te_pr() function)
# hits1_p, hits2_p, hits1_se, hits2_se, hits1_te, hits2_te <- S4Vector object containing the hits between the genomic elements (P, SE, TE) and the loop anchors

# What it does:

# Output
# loops_gi <- GInteraction object containing the ids of the genomic elements that have hits against the loop anchors


giveloop_ids <- function(loops_gi, pr_gr, te_gr, se_gr, 
                         hits1_p, hits2_p,
                         hits1_se, hits2_se,
                         hits1_te, hits2_te){
  
  anchor_wise_ids<-function(df_gr, hits){
    
    df_idx<-mcols(df_gr)$id[subjectHits(hits)]
    
    list<-split(df_idx, queryHits(hits))
    
    return(list)
  }
  
  # ---- extracting the ids of promoters, se and te for which loop anchors have hits ----
  # ---- for promoter ---- 
  
  anchor1_prom_list <- anchor_wise_ids(pr_gr, hits1_p)
  anchor2_prom_list <- anchor_wise_ids(pr_gr, hits2_p)
  
  N <- length(loops_gi)
  
  anchor1_promIDs <- vector("list", N)
  anchor2_promIDs <- vector("list", N)
  
  n <- length(anchor1_prom_list)
  
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  
  for (idx in seq_along(anchor1_prom_list)) {
    i <- as.integer(names(anchor1_prom_list)[idx])
    
    anchor1_promIDs[[i]] <- anchor1_prom_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor1_prom: %d%% (%d/%d)", round(idx/n*100), idx, n))
    
  }
  
  mcols(loops_gi)$anchor1_promoter_ids <- anchor1_promIDs
  
  n <- length(anchor2_prom_list)
  
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  
  for (idx in seq_along(anchor2_prom_list)) {
    
    i <- as.integer(names(anchor2_prom_list)[idx])
    
    anchor2_promIDs[[i]] <- anchor2_prom_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor2_prom: %d%% (%d/%d)", round(idx/n*100), idx, n))
  
  }
  
  mcols(loops_gi)$anchor2_promoter_ids <- anchor2_promIDs
  
  # ---- for super enhancer ----
  anchor1_se_list <- anchor_wise_ids(se_gr, hits1_se)
  anchor2_se_list <- anchor_wise_ids(se_gr, hits2_se)
  
  anchor1_seIDs <- vector("list", N)
  anchor2_seIDs <- vector("list", N)
  
  n <- length(anchor1_se_list)
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  for (idx in seq_along(anchor1_se_list)) {
    
    i <- as.integer(names(anchor1_se_list)[idx])
    
    anchor1_seIDs[[i]] <- anchor1_se_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor1_se: %d%% (%d/%d)", round(idx/n*100), idx, n))
  }
  
  mcols(loops_gi)$anchor1_se_ids <- anchor1_seIDs
  
  n <- length(anchor2_se_list)
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  for (idx in seq_along(anchor2_se_list)) {
    
    i <- as.integer(names(anchor2_se_list)[idx])
    
    anchor2_seIDs[[i]] <- anchor2_se_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor2_se: %d%% (%d/%d)", round(idx/n*100), idx, n))
  }
  
  mcols(loops_gi)$anchor2_se_ids <- anchor2_seIDs
  
  # ---- for typical enhancer ----
  anchor1_te_list <- anchor_wise_ids(te_gr, hits1_te)
  anchor2_te_list <- anchor_wise_ids(te_gr, hits2_te)
  
  anchor1_teIDs <- vector("list", N)
  anchor2_teIDs <- vector("list", N)
  
  n <- length(anchor1_te_list)
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  for (idx in seq_along(anchor1_te_list)) {
    
    i <- as.integer(names(anchor1_te_list)[idx])
    
    anchor1_teIDs[[i]] <- anchor1_te_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor1_te: %d%% (%d/%d)", round(idx/n*100), idx, n))
  }
  
  mcols(loops_gi)$anchor1_te_ids <- anchor1_teIDs
  
  n <- length(anchor2_te_list)
  checkpoints <- unique(floor(seq(0.1, 1.0, by = 0.1) * n))
  for (idx in seq_along(anchor2_te_list)) {
    
    i <- as.integer(names(anchor2_te_list)[idx])
    
    anchor2_teIDs[[i]] <- anchor2_te_list[[idx]]
    
    if (idx %in% checkpoints) message(sprintf("anchor2_te: %d%% (%d/%d)", round(idx/n*100), idx, n))
  }
  
  mcols(loops_gi)$anchor2_te_ids <- anchor2_teIDs
  
  return(loops_gi)
  
}


# ---- 4. Finding the interaction types of each anchor ----

# Function name: classify_interaction(), make_loop_interaction_categories()

# Input:
# loops_gi <- GInteraction object containing the ids of intersecting genomic elements of each loop anchor in metadata
# loops <- dataframe of loops

# What it does: 

# Output:
# loops_gi <- GInteraction object containing the interaction type of each loop observed

# priority of elements (if a loop anchor has hits against SE,TE and P, then this priority would be used to assign a loop anchor to either SE, TE or P) 
# priority chosen here (P>SE>TE)

classify_interaction<-function(p1,se1,te1){
  if(p1){
    return("P")
  }
  else if(se1){
    return("SE")
  }
  else if(te1){
    return("TE")
  }
  else{
    return("other")
  }
}


make_loop_interaction_categories<-function(loops_gi){
  
  # ---- apply the priority of elements to each anchor of the loop ---- 
  df<-as.data.frame(loops_gi)
  
  df$interaction_a1 <- mapply(
    classify_interaction,
    df$is_P_a1,  df$is_SE_a1,  df$is_TE_a1
  )
  
  df$interaction_a2 <- mapply(
    classify_interaction,
    df$is_P_a2,  df$is_SE_a2,  df$is_TE_a2
  )
  
  df$interaction_type <- paste0(df$interaction_a1,"-",df$interaction_a2)
  
  df$interaction_type_2 <- mapply(function(a, b) {
    paste(sort(c(a, b)), collapse = "-")
  }, df$interaction_a1, df$interaction_a2)
  
  mcols(loops_gi)$interaction_type <- df$interaction_type
  mcols(loops_gi)$interaction_type_2 <- df$interaction_type_2
  b<-as.data.frame(mcols(loops_gi)) %>% group_by(interaction_type_2) %>% summarise(count=n())
  
  library(ggplot2)
  
  p1<-ggplot(b, aes(x = interaction_type_2, y = count, fill =interaction_type_2)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = 0, size = 4) +  # Add count labels above bars
    labs(title = "Counts of the loop anchors",
         x = "",
         y = "Counts") +
    theme_minimal()
  
  print(p1)
  
  n_rel_interactions<-as.data.frame(mcols(loops_gi)) %>% select(interaction_type) %>% filter(interaction_type %in% c("P-P","P-SE","P-TE","SE-P","TE-P","TE-TE","SE-SE","SE-TE","TE-SE")) %>% nrow()
  
  cat("\n The number of total relevant loops is : ", n_rel_interactions)
  
  return(loops_gi)
}


# ---- 5. Make the edges and nodes ----

# Function name: give_total_edges_nodes()

# Input:
# loops_g1 <- GRanges object of Loops containing the interaction type
# loops <- loops data frame 

# What it does: for each of the loops that are observes, it removes redundancy in interaction type, 
#               makes it consistent and fetches the relevant ids based on the interaction type and 
#               generates pairwise combinations between the nodes to make the edge dataframe and nodes dataframe
#               generates some graphs

# Output:
# total_nodes <- df containing all the nodes and the type (SE, TE or P)
# total_edges <- df containing all the edges (ids), their interaction frequency and interaction type

give_total_edges_nodes <- function(loops_g1, loops){
  
  library(purrr)
  
  # split the GInteraction object by interaction type
  loops_by_type <- split(loops_g1, mcols(loops_g1)$interaction_type)
  
  # fetch the colname of an anchor based on the type of the anchor
  get_colname <- function(anchor, type) {
    switch(type,
           "P" = paste0("anchor", anchor, "_promoter_ids"),
           "SE" = paste0("anchor", anchor, "_se_ids"),
           "TE" = paste0("anchor", anchor, "_te_ids")
    )
  }
  
  # function to fetch relevant ids for each anchor based on its type 
  loops_by_type_processed <- lapply(names(loops_by_type), function(type_name) {
    
    loop_obj <- loops_by_type[[type_name]]
    
    # Split interaction type 
    
    types <- strsplit(type_name, "-")[[1]]
    type_a1 <- types[1]
    type_a2 <- types[2]
    
    # Get corresponding column names
    
    col_a1 <- get_colname(1, type_a1)
    col_a2 <- get_colname(2, type_a2)
    
    # Subset metadata to only relevant columns
    
    mcols(loop_obj) <- mcols(loop_obj)[, c(col_a1, col_a2, "interaction_type")]
    
    return(loop_obj)
    
  })
  
  # Keep names
  names(loops_by_type_processed) <- names(loops_by_type)
  
  # Determine the anchor types of A1 and A2, since they were not explicitly mentioned in the previous function 
  
  get_anchor_type <- function(isP, isSE, isTE) {
    if (isP) return("P")
    if (isSE) return("SE")
    if (isTE) return("TE")
    return("other")
  }
  
  priority <- c("P" = 1, "SE" = 2, "TE" = 3, "other" = 4)
  
  loops_df <- as.data.frame(loops_g1)
  
  # Get anchor types of anchor1 and 2 as T/F
  
  loops_df$type_a1 <- mapply(get_anchor_type,
                             loops_df$is_P_a1,
                             loops_df$is_SE_a1,
                             loops_df$is_TE_a1)
  
  loops_df$type_a2 <- mapply(get_anchor_type,
                             loops_df$is_P_a2,
                             loops_df$is_SE_a2,
                             loops_df$is_TE_a2)
  
  # Decide if swap needed (which loops require swaps) why are we swapping? - to make our interaction_type consistent
  
  loops_df$swap <- priority[loops_df$type_a1] > priority[loops_df$type_a2]
  
  swap_cols <- function(df, col1, col2, idx) {
    tmp <- df[[col1]][idx]
    df[[col1]][idx] <- df[[col2]][idx]
    df[[col2]][idx] <- tmp
    df
  }
  
  swap_idx <- which(loops_df$swap) # the rows which need to be swapped
  
  cols_to_swap <- list(
    c("is_P_a1", "is_P_a2"),
    c("is_SE_a1", "is_SE_a2"),
    c("is_TE_a1", "is_TE_a2"),
    c("anchor1_promoter_ids", "anchor2_promoter_ids"),
    c("anchor1_se_ids", "anchor2_se_ids"),
    c("anchor1_te_ids", "anchor2_te_ids")
  )
  
  # perform the swap
  for (pair in cols_to_swap) {
    loops_df <- swap_cols(loops_df, pair[1], pair[2], swap_idx)
  }
  
  # recompute the anchor types
  loops_df$type_a1 <- mapply(get_anchor_type,
                             loops_df$is_P_a1,
                             loops_df$is_SE_a1,
                             loops_df$is_TE_a1)
  
  loops_df$type_a2 <- mapply(get_anchor_type,
                             loops_df$is_P_a2,
                             loops_df$is_SE_a2,
                             loops_df$is_TE_a2)
  
  # rebuild interaction type
  loops_df$interaction_type <- paste(loops_df$type_a1,
                                     loops_df$type_a2,
                                     sep = "-")
  
  unique(loops_df$interaction_type)
  
  # function to extract correct IDs for each loop based on its interaction type
  get_ids <- function(row, anchor, type) {
    col <- switch(type,
                  "P" = paste0("anchor", anchor, "_promoter_ids"),
                  "SE" = paste0("anchor", anchor, "_se_ids"),
                  "TE" = paste0("anchor", anchor, "_te_ids"),
                  "other" = NULL
    )
    
    if (is.null(col)) return(NA)
    
    ids <- row[[col]]
    if (length(ids[[1]]) == 0) return(NA)
    
    return(ids[[1]])
  }
  
  # extract ids for each loop
  loops_df <- loops_df %>%
    rowwise() %>%
    mutate(
      ids_a1 = list(get_ids(cur_data(), 1, type_a1)),
      ids_a2 = list(get_ids(cur_data(), 2, type_a2))
    ) %>%
    ungroup()
  
  library(dplyr)
  
  n <- nrow(loops_df)
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  ids_a1 <- vector("list", n)
  ids_a2 <- vector("list", n)
  
  for (i in seq_len(n)) {
    ids_a1[[i]] <- get_ids(loops_df[i, ], 1, loops_df$type_a1[i])
    ids_a2[[i]] <- get_ids(loops_df[i, ], 2, loops_df$type_a2[i])
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  loops_df$ids_a1 <- ids_a1
  loops_df$ids_a2 <- ids_a2
  
  # create the edges
  loops_df %>% head()
  loops_df %>% select(interaction_type_2) %>% unique()
  loops_df %>% select(interaction_type) %>% unique()
  
  # loops_df <- loops_df %>% filter(interaction_type %in% c("P-TE", "P-SE", "SE-TE", "P-P", "SE-TE", "TE-TE", "SE-SE" ))
  loops_df <- loops_df %>% filter(interaction_type %in% c("P-TE","TE-P", 
                                                          "P-SE","SE-P",
                                                          "SE-TE","TE-SE",
                                                          "P-P", "TE-TE", "SE-SE" ))
  
  library(data.table)
  
  df <- loops_df %>%
    filter(!is.na(ids_a1), !is.na(ids_a2)) %>%
    select(interaction_type, ids_a1, ids_a2, value)
  
  chunk_size <- 5000
  
  n <- nrow(df)
  
  edges_list <- list()
  
  k <- 1
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  for (start in seq(1, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n)
    
    chunk <- df[start:end, ]
    
    chunk_edges <- lapply(seq_len(nrow(chunk)), function(i) {
      if (length(chunk$ids_a1[[i]]) == 0 || length(chunk$ids_a2[[i]]) == 0) return(NULL)
      
      data.frame(
        node1 = rep(chunk$ids_a1[[i]], each = length(chunk$ids_a2[[i]])),
        node2 = rep(chunk$ids_a2[[i]], times = length(chunk$ids_a1[[i]])),
        interaction_type = chunk$interaction_type[i],
        value = chunk$value[i]
      )
    })
    
    edges_list[[k]] <- data.table::rbindlist(chunk_edges)
    k <- k + 1
    
    setTxtProgressBar(pb, end)
  }
  
  close(pb)
  
  edges <- data.table::rbindlist(edges_list)
  edges %>% nrow()
  edges %>% distinct %>% nrow()
  
  unique(edges$interaction_type)
  edges <- edges %>% distinct()
  edges <- edges %>% select(node1, node2, value, interaction_type)
  colnames(edges) <- c("LA1", "LA2","value", "interaction_type")
  
  edges_unique <- edges
  edges_unique <- edges %>%
    group_by(
      LA1, LA2, interaction_type
    ) %>%
    summarise(
      value_mean = mean(value, na.rm = TRUE),
      n_loops = n(),
      .groups = "drop"
    )
  
  edges_unique %>% head()
  edges_unique <- edges_unique %>% dplyr::rename(value=value_mean) %>% select(1,2,4,3)
  edges <- edges_unique
  
  edges %>%
    group_by(
      LA1, LA2, interaction_type
    ) %>%
    summarise(
      value_mean = max(value, na.rm = TRUE),
      n_loops = n(),
      .groups = "drop"
    ) %>% nrow()
  
  
  edges %>% nrow()
  edges %>% filter(LA1!=LA2) %>% nrow()
  
  edges %>%
    group_by(
      LA1, LA2, interaction_type
    ) %>%
    summarise(
      value_mean = mean(value, na.rm = TRUE),
      n_loops = n(),
      .groups = "drop"
    ) %>% filter(LA1!=LA2) %>% nrow()
  
  edges <- edges %>% filter(LA1!=LA2)
  
  # nodes
  a <- edges %>% separate(interaction_type, into = c("LA1_type", "LA2_type"), sep = "-")
  a1<-a %>% select(LA1_type, LA1)
  a2<-a %>% select(LA2_type, LA2)
  colnames(a2)<-colnames(a1)
  nodes <- rbind(a1,a2) %>% unique()
  colnames(nodes) <- c("LA1_type", "LA1")
  nodes <- nodes %>% select(LA1, LA1_type)
  table(nodes$LA1_type)
  
  total_edges <- edges
  total_nodes <- nodes
  
  total_edges %>% group_by(interaction_type) %>% summarise(count=n())
  
  m<-total_edges %>% group_by(interaction_type) %>% summarise(count=n()) 
  
  p1<-ggplot(m, aes(x =interaction_type, y = count, fill =interaction_type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = 0, size = 4) +  # Add count labels above bars
    labs(title = "Counts of interacting nodes",
         x = "",
         y = "Counts") +
    theme_minimal()
  
  print(p1)
  
  return(list(total_nodes=total_nodes, total_edges=total_edges))
}

# ---- 6. Construct the network  ----

# Function name: construct_nw()

# Input:
# total_nodes <- data frame containing the nodes
# total_edges <- data frame containing the edges
# path <- output path

# What it does: generates a weighted or unweighted graph, using the input provided
#               identifies the connected components in the graph
#               sorts the clusters based on the kind of nodes it has into cluster categories
#               identifies which node belong to which cluster & cluster type
#               prints couple of graphs


# Output:
# graph <- (igraph object) enhancer promoter graph
# df <- number of P, SE and TE in the network
# df1 <- how many clusters are identifies in each cluster category
# df2 <- sizes of clusters belonging to each cluster category
# cluster_classification <- which cluster belongs to which category
# components_g <- output of components(g)

# construct_nw<-function(total_nodes,total_edges, path){
#   
#   library(igraph)
#   library(tidyverse)
#   
#   g<-graph_from_data_frame(d = total_edges, directed = F, vertices = total_nodes)
#   
#   if("value" %in% edge_attr_names(g)){
#     E(g)$weight<-E(g)$value
#   }else{
#   }
#   
#   
#   # graph attributes
#   # number of nodes of each type
#   a<-total_nodes %>% filter(LA1_type=="P") %>% unique() %>% nrow() # 2153
#   b<-total_nodes %>% filter(LA1_type=="SE") %>% unique() %>% nrow() # 3406
#   c<-total_nodes %>% filter(LA1_type=="TE") %>% unique() %>% nrow() # 6221
#   
#   df<-data.frame(char=c("Promoter", "SE constituent","Typical Enhancer"), Count=c(a,b,c))
#   
#   # calculating the components 
#   components_g <- components(g)
#   
#   total_cwise <- data.frame(Node=names(components_g$membership),
#                             Cluster=components_g$membership)
#   
#   clu_wise_enhancers<-total_cwise %>% left_join(total_nodes,by=c("Node"="LA1")) %>% filter(!LA1_type=="promoter")
#   
#   clu_wise_nodes<-total_cwise %>% left_join(total_nodes,by=c("Node"="LA1"))
#   
#   # classifying the clusters based on their constituents
#   cluster_classification<-clu_wise_nodes %>%
#     group_by(Cluster) %>%
#     summarise(
#       Types = list(unique(LA1_type))        # all LA1_type values in this cluster
#     ) %>%  mutate(type_set = map(Types, ~ sort(unique(.x)))) %>%
#     mutate(
#       cluster_type = map_chr(type_set, function(ts) {
#         # ts is the sorted unique vector for one cluster
#         if (identical(ts, "P")) {
#           "P"                # only P
#         } else if (identical(ts, "TE")) {
#           "TE"               # only TE
#         } else if (identical(ts, "SE")) {
#           "SE"               # only SE
#         } else if (identical(ts, c("P","SE"))) {
#           "P-SE"
#         } else if (identical(ts, c("P","TE"))) {
#           "P-TE"
#         } else if (identical(ts, c("SE","TE"))) {
#           "SE-TE"
#         } else if (identical(ts, c("P","SE","TE"))) {
#           "P-SE-TE"
#         } else {
#           paste(ts, collapse="-")
#         }
#       })
#     )
#   
#   
#   df1<-data.frame(Cluster_Type=character(), Count=integer(), stringsAsFactors = F)
#   
#   for(i in unique(cluster_classification$cluster_type)){
#     # print(i)
#     m<-cluster_classification %>% filter(cluster_type==i) %>% unique() %>% nrow()
#     # print(m)
#     df1<-rbind(df1, data.frame(Cluster_Type=i, Count=m))
#   }
#   
#   clusters<-components(g)
#   cluster_sizes<-clusters$csize
#   df3 <- data.frame(cluster_size = cluster_sizes) %>% mutate(Cluster=row_number())
#   
#   cluster_classification %>% head()
#   df2<-df3 %>% left_join(cluster_classification,by="Cluster",relationship = "many-to-many")
#   df2<-df2 %>% select(2,1,5)
#   df2 %>% head()
#   
#   p1 <- ggplot(df, aes(x = char, y = Count, fill = char)) +
#     geom_bar(stat = "identity") +
#     geom_text(aes(label = Count), vjust = -0.3, size = 4) +  # Slightly lift labels above bars
#     labs(
#       title = "Counts of the nodes",
#       x = "",
#       y = "Count"
#     ) +
#     scale_fill_manual(values = c(
#       "SE constituent" = "#af8dc3",
#       "Typical Enhancer" = "#5ab4ac",
#       "Promoter" = "#C3AF8D"
#     )) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Bars touch x-axis; space above
#     theme_minimal(base_size = 14) +
#     theme(
#       panel.background = element_rect(fill = "white", colour = NA),
#       plot.background = element_rect(fill = "white", colour = NA),
#       legend.background = element_rect(fill = "white", colour = NA),
#       legend.box.background = element_rect(fill = "white", colour = NA),
#       panel.grid = element_blank(),
#       legend.position = "none",
#       axis.line = element_line(color = "black", linewidth = 0.6),
#       axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
#       axis.title.x = element_text(face = "bold", size = 16),
#       axis.title.y = element_text(face = "bold", size = 16),
#       plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
#     )
#   print(p1)
#   
#   
#   
#   p1<- ggplot(df1, aes(x = Cluster_Type, y = Count, fill = Cluster_Type)) +
#     geom_bar(stat = "identity") +
#     geom_text(aes(label = Count), vjust = -0.3, size = 4) +
#     labs(title = "Number of clusters in each category",
#          x = "",
#          y = "Number of clusters") +
#     scale_fill_manual(
#       values = c(
#         "P-SE" = "#af8dc3",
#         "P-TE" = "#5ab4ac",
#         "P-SE-TE" = "#84a0b8"
#       ),
#       na.value = "lightgray"
#     ) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#     theme_minimal(base_size = 14) +
#     theme(
#       panel.background = element_rect(fill = "white", colour = NA),
#       plot.background = element_rect(fill = "white", colour = NA),
#       legend.background = element_rect(fill = "white", colour = NA),
#       legend.box.background = element_rect(fill = "white", colour = NA),
#       panel.grid = element_blank(),
#       legend.position = "none",
#       axis.line = element_line(color = "black", linewidth = 0.6),
#       axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
#       axis.title.x = element_text(face = "bold", size = 16),
#       axis.title.y = element_text(face = "bold", size = 16),
#       plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
#     )
#   
#   
#   print(p1)
#   
#   
#   p1<-ggplot(df2, aes(x = cluster_size)) +
#     geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
#     scale_x_continuous(breaks = seq(min(df2$cluster_size), max(df2$cluster_size), by = 1)) +  # Set x labels at interval of 1
#     labs(title = "Cluster Size Distribution",
#          x = "Cluster Size",
#          y = "Frequency") +
#     theme_minimal()
#   
#   print(p1)
#   
#   
#   p1<-ggplot(df2, aes(x = cluster_type, y = cluster_size, fill = cluster_type)) +
#     geom_boxplot() +
#     labs(
#       title = "Cluster Size Distribution by Cluster Category",
#       x = "",
#       y = "Cluster Size"
#     ) +
#     scale_fill_manual(values = c(
#       "P-SE" = "#af8dc3",
#       "P-TE" = "#5ab4ac",
#       "P-SE-TE" = "#84a0b8"
#     ),
#     na.value = "lightgray"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#       panel.background = element_rect(fill = "white", colour = NA),
#       plot.background = element_rect(fill = "white", colour = NA),
#       legend.background = element_rect(fill = "white", colour = NA),
#       legend.box.background = element_rect(fill = "white", colour = NA),
#       panel.grid = element_blank(),
#       legend.position = "none",
#       axis.line = element_line(color = "black", linewidth  = 0.6),
#       axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
#       axis.title.x = element_text(face = "bold", size = 16),
#       axis.title.y = element_text(face = "bold", size = 16),
#       plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
#     )
#   print(p1)
#   
#   dir.create(path, showWarnings = FALSE)
#   
#   invisible(saveRDS(object = g, file = paste0(path,"Graph.rds")))
#   
#   return(list(graph=g, df=df, df1=df1, df2=df2,  cluster_classification=cluster_classification, components_g=components_g))
# }

construct_nw<-function(total_nodes,total_edges, path, bipartite=F){
  
  # library(igraph)
  # library(tidyverse)
  
  if(bipartite==T){
    total_edges <- total_edges %>% filter(interaction_type%in%c("P-SE", "P-TE"))
    a <- total_edges$LA1
    b <- total_edges$LA2
    nodes <- c(a,b)
    total_nodes <- total_nodes %>% filter(LA1 %in% nodes)
  }
  
  
  g<-graph_from_data_frame(d = total_edges, directed = F, vertices = total_nodes)
  
  if("value" %in% edge_attr_names(g)){
    E(g)$weight<-E(g)$value
  }else{
  }
  
  
  # graph attributes
  # number of nodes of each type
  a<-total_nodes %>% filter(LA1_type=="P") %>% unique() %>% nrow() # 2153
  b<-total_nodes %>% filter(LA1_type=="SE") %>% unique() %>% nrow() # 3406
  c<-total_nodes %>% filter(LA1_type=="TE") %>% unique() %>% nrow() # 6221
  
  df<-data.frame(char=c("Promoter", "SE constituent","Typical Enhancer"), Count=c(a,b,c))
  
  # calculating the components 
  components_g <- components(g)
  
  total_cwise <- data.frame(Node=names(components_g$membership),
                            Cluster=components_g$membership)
  
  clu_wise_enhancers<-total_cwise %>% left_join(total_nodes,by=c("Node"="LA1")) %>% filter(!LA1_type=="promoter")
  
  clu_wise_nodes<-total_cwise %>% left_join(total_nodes,by=c("Node"="LA1"))
  
  # classifying the clusters based on their constituents
  cluster_classification<-clu_wise_nodes %>%
    group_by(Cluster) %>%
    summarise(
      Types = list(unique(LA1_type))        # all LA1_type values in this cluster
    ) %>%  mutate(type_set = map(Types, ~ sort(unique(.x)))) %>%
    mutate(
      cluster_type = map_chr(type_set, function(ts) {
        # ts is the sorted unique vector for one cluster
        if (identical(ts, "P")) {
          "P"                # only P
        } else if (identical(ts, "TE")) {
          "TE"               # only TE
        } else if (identical(ts, "SE")) {
          "SE"               # only SE
        } else if (identical(ts, c("P","SE"))) {
          "P-SE"
        } else if (identical(ts, c("P","TE"))) {
          "P-TE"
        } else if (identical(ts, c("SE","TE"))) {
          "SE-TE"
        } else if (identical(ts, c("P","SE","TE"))) {
          "P-SE-TE"
        } else {
          paste(ts, collapse="-")
        }
      })
    )
  
  
  df1<-data.frame(Cluster_Type=character(), Count=integer(), stringsAsFactors = F)
  
  for(i in unique(cluster_classification$cluster_type)){
    
    m<-cluster_classification %>% filter(cluster_type==i) %>% unique() %>% nrow()
    
    df1<-rbind(df1, data.frame(Cluster_Type=i, Count=m))
  }
  
  clusters<-components(g)
  cluster_sizes<-clusters$csize
  df3 <- data.frame(cluster_size = cluster_sizes) %>% mutate(Cluster=row_number())
  
  cluster_classification %>% head()
  df2<-df3 %>% left_join(cluster_classification,by="Cluster",relationship = "many-to-many")
  df2<-df2 %>% select(2,1,5)
  df2 %>% head()
  
  p1 <- ggplot(df, aes(x = char, y = Count, fill = char)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.3, size = 4) +  # Slightly lift labels above bars
    labs(
      title = "Counts of the nodes",
      x = "",
      y = "Count"
    ) +
    scale_fill_manual(values = c(
      "SE constituent" = "#af8dc3",
      "Typical Enhancer" = "#5ab4ac",
      "Promoter" = "#C3AF8D"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Bars touch x-axis; space above
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  print(p1)
  
  
  
  p1<- ggplot(df1, aes(x = Cluster_Type, y = Count, fill = Cluster_Type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.3, size = 4) +
    labs(title = "Number of clusters in each category",
         x = "",
         y = "Number of clusters") +
    scale_fill_manual(
      values = c(
        "P-SE" = "#af8dc3",
        "P-TE" = "#5ab4ac",
        "P-SE-TE" = "#84a0b8"
      ),
      na.value = "lightgray"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  
  print(p1)
  
  
  p1<-ggplot(df2, aes(x = cluster_size)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
    scale_x_continuous(breaks = seq(min(df2$cluster_size), max(df2$cluster_size), by = 1)) +  # Set x labels at interval of 1
    labs(title = "Cluster Size Distribution",
         x = "Cluster Size",
         y = "Frequency") +
    theme_minimal()
  
  print(p1)
  
  
  p1<-ggplot(df2, aes(x = cluster_type, y = cluster_size, fill = cluster_type)) +
    geom_boxplot() +
    labs(
      title = "Cluster Size Distribution by Cluster Category",
      x = "",
      y = "Cluster Size"
    ) +
    scale_fill_manual(values = c(
      "P-SE" = "#af8dc3",
      "P-TE" = "#5ab4ac",
      "P-SE-TE" = "#84a0b8"
    ),
    na.value = "lightgray"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth  = 0.6),
      axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(face = "bold"),     # <- Bold x-axis labels
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  print(p1)
  
  # dir.create(path, showWarnings = FALSE)
  
  invisible(saveRDS(object = g, file = paste0(path,"Graph.rds")))
  
  return(list(graph=g, df=df, df1=df1, df2=df2,  cluster_classification=cluster_classification, components_g=components_g))
}

# ---- 7. Calculate the network properties and plot ----

# Function name: calculate_nw_ppts()

# Input:
# g <- graph (igraph object)
# components_g <- output of components(g)
# cluster_classification <-  which cluster belongs to which category
# path <- output path

# What it does: for the nodes in each cluster, 
#               it computes following network properties: degree, strength, betweenness, closeness, harmonic & eigenvector centrality 
#               filters only the enhancers

# Output:
# enhancer_props_df <- network properties of each enhancer in the network


# calculate_nw_ppts<-function(g, components_g, cluster_classification, path){
#   
#   library(igraph)
#   library(dplyr)
#   library(tibble)
#   library(purrr)  # for bind_rows/map if needed
#   
#   cat("\nThe total number of components in the network is:", components_g$no)
#   cat("\nThe largest cluster has ", components_g$csize %>% max(), " nodes ")
#   enhancer_props_list <- list()
#   
#   for (i in unique(components_g$membership)) {
#     
#     # 1. Get nodes in component i and build the subgraph
#     subgraph_nodes <- which(components_g$membership == i)
#     subgraph <- induced_subgraph(g, vids = subgraph_nodes)
#     
#     # 2. Compute network properties on this subgraph
#     deg   <- degree(subgraph, normalized = FALSE)
#     norm_deg <- degree(subgraph, normalized = TRUE)
#     
#     if ("value" %in% edge_attr_names(subgraph)){
#       E(subgraph)$weight<-E(subgraph)$value
#       strength <- igraph::strength(subgraph)
#       E(subgraph)$inv_weight<-1/(E(subgraph)$value)
#       
#       clo <- closeness(subgraph,weights = E(subgraph)$inv_weight, normalized = F)
#       clo_norm <- closeness(subgraph,weights = E(subgraph)$inv_weight, normalized = T)
#       
#       bet <- betweenness(subgraph,weights = E(subgraph)$inv_weight,normalized = F)
#       bet_norm <- betweenness(subgraph,weights = E(subgraph)$inv_weight,normalized = T)
#       
#       harmonic <- harmonic_centrality(subgraph,weights = E(subgraph)$inv_weight,normalized = F)
#       harmonic_norm <- harmonic_centrality(subgraph,weights = E(subgraph)$inv_weight,normalized = T)
#       
#       eigen <- eigen_centrality(subgraph,weights = E(subgraph)$weight)$vector
#       
#     } else{
#       
#       strength <- strength(subgraph)
#       
#       
#       clo <- closeness(subgraph, normalized = F)
#       clo_norm <- closeness(subgraph, normalized = T)
#       
#       bet <- betweenness(subgraph,normalized = F)
#       bet_norm <- betweenness(subgraph,normalized = T)
#       
#       harmonic <- harmonic_centrality(subgraph,normalized = F)
#       harmonic_norm <- harmonic_centrality(subgraph,normalized = T)
#       
#       eigen <- eigen_centrality(subgraph)$vector
#     }
#     
#     # 3. Turn vertex attributes into a tibble
#     subgraph_vertices <- as_tibble(
#       igraph::as_data_frame(subgraph, what = "vertices")
#     )
#     
#     # 4. Filter only enhancers, then add centrality measures
#     enh_sub <- subgraph_vertices %>%
#       filter(LA1_type %in% c("SE", "TE")) %>%
#       mutate(
#         Cluster     = i,
#         degree      = deg[name],
#         norm_degree = norm_deg[name],
#         strength    = strength[name],
#         betweenness = bet[name],
#         betweenness_norm = bet_norm[name],
#         closeness   = clo[name],
#         closeness_norm = clo_norm[name], 
#         eigenvector = eigen[name],
#         harmonic =harmonic[name],
#         harmonic_norm=harmonic_norm[name]
#       )
#     
#     # 5. Store tibble for this component
#     enhancer_props_list[[as.character(i)]] <- enh_sub
#   }
#   
#   # 6. Combine all components into one dataframe
#   enhancer_props_df <- bind_rows(enhancer_props_list)
#   
#   # View first few rows
#   enhancer_props_df %>% head()
#   
#   
#   # for each of the network property plot the graph and perform the wilcox test
#   library(dplyr)
#   
#   properties <- enhancer_props_df %>%
#     select(-name, -LA1_type, -Cluster) %>%
#     select(where(is.numeric)) %>%
#     names()
#   
#   write.csv(enhancer_props_df,paste0(path,"enhancer_properties_df.csv"),quote = F, row.names = F )
#   return(enhancer_props_df)
# }

calculate_nw_ppts<-function(g, components_g, cluster_classification, path){

  cat("\nThe total number of components in the network is:", components_g$no)
  cat("\nThe largest cluster has ", components_g$csize %>% max(), " nodes ")
  enhancer_props_list <- list()
  
  # small helper to sanitize a weight vector to be strictly positive & finite
  sanitize_weights <- function(w, min_pos = .Machine$double.eps) {
    w[is.na(w)]        <- min_pos
    w[is.infinite(w)]  <- max(w[is.finite(w)], na.rm = TRUE)  # cap Inf at max finite
    w[w <= 0]          <- min_pos
    w
  }
  
  enhancer_props_list <- list()
  
  for (i in unique(components_g$membership)) {
    
    # 1. Subgraph for component i
    subgraph_nodes <- which(components_g$membership == i)
    subgraph      <- induced_subgraph(g, vids = subgraph_nodes)
    
    # 2. Degree (always works)
    deg      <- degree(subgraph, normalized = FALSE)
    # norm_deg <- degree(subgraph, normalized = TRUE)
    
    if ("value" %in% edge_attr_names(subgraph) && ecount(subgraph) > 0) {
      
      raw_val <- E(subgraph)$value
      
      # weight for strength / eigen: positive, finite
      w_pos <- sanitize_weights(raw_val)
      E(subgraph)$weight <- w_pos
      
      # inverse weight for distance-based centralities (closeness, betweenness, harmonic)
      inv_w <- 1 / w_pos
      inv_w <- sanitize_weights(inv_w)
      E(subgraph)$inv_weight <- inv_w
      
      strength <- igraph::strength(subgraph, weights = w_pos)
      
      clo           <- tryCatch(closeness(subgraph, weights = inv_w, normalized = FALSE),
                                error = function(e) rep(NA_real_, vcount(subgraph)))
      # clo_norm      <- tryCatch(closeness(subgraph, weights = inv_w, normalized = TRUE),
                                # error = function(e) rep(NA_real_, vcount(subgraph)))
      
      bet           <- tryCatch(betweenness(subgraph, weights = inv_w, normalized = FALSE),
                                error = function(e) rep(NA_real_, vcount(subgraph)))
      # bet_norm      <- tryCatch(betweenness(subgraph, weights = inv_w, normalized = TRUE),
                                # error = function(e) rep(NA_real_, vcount(subgraph)))
      
      harmonic      <- tryCatch(harmonic_centrality(subgraph, weights = inv_w, normalized = FALSE),
                                error = function(e) rep(NA_real_, vcount(subgraph)))
      # harmonic_norm <- tryCatch(harmonic_centrality(subgraph, weights = inv_w, normalized = TRUE),
                                # error = function(e) rep(NA_real_, vcount(subgraph)))
      
      eigen <- tryCatch(eigen_centrality(subgraph, weights = w_pos)$vector,
                        error = function(e) rep(NA_real_, vcount(subgraph)))
      
    } else {
      
      strength      <- igraph::strength(subgraph)
      clo           <- closeness(subgraph, normalized = FALSE)
      # clo_norm      <- closeness(subgraph, normalized = TRUE)
      bet           <- betweenness(subgraph, normalized = FALSE)
      # bet_norm      <- betweenness(subgraph, normalized = TRUE)
      harmonic      <- harmonic_centrality(subgraph, normalized = FALSE)
      # harmonic_norm <- harmonic_centrality(subgraph, normalized = TRUE)
      eigen         <- eigen_centrality(subgraph)$vector
    }
    
    # 3. Vertex attributes
    subgraph_vertices <- as_tibble(
      igraph::as_data_frame(subgraph, what = "vertices")
    )
    
    # 4. Filter enhancers + attach centralities
    enh_sub <- subgraph_vertices %>%
      filter(LA1_type %in% c("SE", "TE")) %>%
      mutate(
        Cluster          = i,
        degree           = deg[name],
        # norm_degree      = norm_deg[name],
        strength         = strength[name],
        betweenness      = bet[name],
        # betweenness_norm = bet_norm[name],
        closeness        = clo[name],
        # closeness_norm   = clo_norm[name],
        eigenvector      = eigen[name],
        harmonic         = harmonic[name]
        # harmonic_norm    = harmonic_norm[name]
      )
    
    enhancer_props_list[[as.character(i)]] <- enh_sub
  }
  
  # 6. Combine all components into one dataframe
  enhancer_props_df <- bind_rows(enhancer_props_list)
  
  # View first few rows
  enhancer_props_df %>% head()
  
  
  # for each of the network property plot the graph and perform the wilcox test
  library(dplyr)
  
  properties <- enhancer_props_df %>%
    select(-name, -LA1_type, -Cluster) %>%
    select(where(is.numeric)) %>%
    names()
  
  # write.csv(enhancer_props_df,paste0(path,"enhancer_properties_df.csv"),quote = F, row.names = F )
  return(enhancer_props_df)
}


# --- for pretty_plots ----

# plot_func <- function(data_frame, metric_name, y_axis_name, by_inp, file_prefix){
#   
#   library(dplyr)
#   library(ggplot2)
#   
#   #-----------------------------
#   # Helper: safe Wilcoxon test
#   #-----------------------------
#   safe_wilcox <- function(x, y, label="") {
#     x <- x[!is.na(x)]
#     y <- y[!is.na(y)]
#     
#     cat(label, "\n")
#     cat("n_x =", length(x), "| n_y =", length(y), "\n")
#     
#     if (length(x) > 1 && length(y) > 1) {
#       res <- wilcox.test(x, y)
#       cat("W statistic:", res$statistic, "\n")
#       cat("p-value:", res$p.value, "\n")
#     } else {
#       cat("Skipping test: not enough data\n")
#     }
#   }
#   
#   #-----------------------------
#   # Filter relevant data
#   #-----------------------------
#   rel_data <- data_frame %>%
#     filter(cluster_type %in% c("P-SE", "P-TE","P-SE-TE"),
#            LA1_type %in% c("SE", "TE"))
#   
#   #-----------------------------
#   # Print sample sizes (important)
#   #-----------------------------
#   cat("\n=== Sample sizes ===\n")
#   rel_data %>%
#     group_by(cluster_type, LA1_type) %>%
#     summarise(n = sum(!is.na(.data[[metric_name]])), .groups="drop") %>%
#     print()
#   
#   #-----------------------------
#   # Compute global whiskers
#   #-----------------------------
#   whiskers <- rel_data %>%
#     group_by(cluster_type, LA1_type) %>%
#     summarise(
#       Q1 = quantile(.data[[metric_name]], 0.25, na.rm = TRUE),
#       Q3 = quantile(.data[[metric_name]], 0.75, na.rm = TRUE),
#       IQR = Q3 - Q1,
#       lower_whisker = Q1 - 1.5 * IQR,
#       upper_whisker = Q3 + 1.5 * IQR,
#       .groups = "drop"
#     )
#   
#   global_min <- min(whiskers$lower_whisker, na.rm = TRUE)
#   global_max <- max(whiskers$upper_whisker, na.rm = TRUE)
#   range_padding <- 0.05 * (global_max - global_min)
#   
#   y_upper <- global_max + range_padding
#   
#   rel_data$cluster_type <- factor(rel_data$cluster_type, 
#                                   levels = c("P-SE", "P-TE", "P-SE-TE"))
#   
#   #-----------------------------
#   # Plot 1: cluster-wise
#   #-----------------------------
#   p1 <- ggplot(rel_data, aes(x = cluster_type, y = .data[[metric_name]], fill = LA1_type)) +
#     geom_boxplot(outlier.shape = NA) +
#     scale_fill_manual(values = c("SE" = "#af8dc3", "TE" = "#5ab4ac")) +
#     coord_cartesian(ylim = c(0, y_upper)) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
#     theme_minimal() +
#     theme(
#       legend.position = "none",
#       panel.background = element_rect(fill = "white"),
#       plot.background = element_rect(fill = "white"),
#       panel.grid = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 0.6),
#       axis.line.y = element_line(color = "black", linewidth = 0.6),
#       axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
#       axis.text.x = element_text(face = "bold", size = 12),
#       axis.title.y = element_text(face = "bold", size = 16)
#     ) +
#     labs(x = "", y = y_axis_name)
#   
#   ggsave(paste0(file_prefix, "_p1.svg"), plot = p1, width = 5, height = 4, dpi = 1200)
#   
#   #-----------------------------
#   # Plot 2: SE vs TE overall
#   #-----------------------------
#   p2 <- ggplot(rel_data, aes(x = LA1_type, y = .data[[metric_name]], fill = LA1_type)) +
#     geom_boxplot(outlier.shape = NA) +
#     scale_fill_manual(values = c("SE" = "#af8dc3", "TE" = "#5ab4ac")) +
#     coord_cartesian(ylim = c(0, y_upper)) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
#     theme_minimal() +
#     theme(
#       legend.position = "none",
#       panel.background = element_rect(fill = "white"),
#       plot.background = element_rect(fill = "white"),
#       panel.grid = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 0.6),
#       axis.line.y = element_line(color = "black", linewidth = 0.6),
#       axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.5),
#       axis.text.x = element_text(face = "bold", size = 12),
#       axis.title.y = element_text(face = "bold", size = 16)
#     ) +
#     labs(x = "", y = y_axis_name)
#   
#   ggsave(paste0(file_prefix, "_p2.svg"), plot = p2, width = 5, height = 4, dpi = 1200)
#   
#   print(p1)
#   print(p2)
#   
#   #-----------------------------
#   # Statistical tests
#   #-----------------------------
#   
#   cat("\n=== Overall Comparison: SE vs TE ===\n")
#   safe_wilcox(
#     rel_data %>% filter(LA1_type == "SE") %>% pull(.data[[metric_name]]),
#     rel_data %>% filter(LA1_type == "TE") %>% pull(.data[[metric_name]])
#   )
#   
#   cat("\n--- P-SE vs P-TE ---\n")
#   safe_wilcox(
#     rel_data %>% filter(cluster_type == "P-SE") %>% pull(.data[[metric_name]]),
#     rel_data %>% filter(cluster_type == "P-TE") %>% pull(.data[[metric_name]])
#   )
#   
#   cat("\n--- Within P-SE-TE (SE vs TE) ---\n")
#   safe_wilcox(
#     rel_data %>% filter(cluster_type == "P-SE-TE", LA1_type == "SE") %>% pull(.data[[metric_name]]),
#     rel_data %>% filter(cluster_type == "P-SE-TE", LA1_type == "TE") %>% pull(.data[[metric_name]])
#   )
# }

