
calculate_pci <- function(cell_line, csv_path, ed1, ed2, ed3, condition_levels){
  
  ed1 <- read_csv(ed1) %>% mutate(condition=condition_levels[1], cell_line=cell_line)
  ed2 <- read_csv(ed2)%>% mutate(condition=condition_levels[2], cell_line=cell_line)
  ed3 <- read_csv(ed3)%>% mutate(condition=condition_levels[3], cell_line=cell_line)
  
  
  edges_raw=rbind(ed1,ed2, ed3)
  
  # edge_path <- "results/k562/wt0/edges_wt0.csv"        # <-- your edge/interaction file
  # edges_raw <- read_csv(edge_path, show_col_types = FALSE)
  
  edges <- edges_raw %>%
    filter(cell_line == cell_line) %>%
    mutate(condition = factor(condition, levels = condition_levels, ordered = TRUE))
  
  # treat network as undirected: a node's partners = everything it shares an edge with.
  # duplicate edges in both directions so each endpoint "sees" the other.
  partner_long <- bind_rows(
    edges %>% select(condition, node = LA1, partner = LA2),
    edges %>% select(condition, node = LA2, partner = LA1)
  ) %>% distinct()
  
  # partner SET per node per condition (as a list-column)
  partner_sets <- partner_long %>%
    group_by(condition, node) %>%
    summarise(partners = list(unique(partner)), .groups = "drop")
  
  jaccard <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0) return(NA_real_)          # isolated in both -> undefined
    length(intersect(a, b)) / u
  }
  
  # PCI between two conditions, for nodes present in BOTH (shared nodes)
  partner_conservation <- function(condA, condB) {
    sa <- partner_sets %>% filter(condition == condA) %>% select(node, pA = partners)
    sb <- partner_sets %>% filter(condition == condB) %>% select(node, pB = partners)
    inner_join(sa, sb, by = "node") %>%
      mutate(pci = purrr::map2_dbl(pA, pB, jaccard),
             transition = paste0(condA, " \u2192 ", condB)) %>%
      select(node, pci, transition)
  }
  
  # run over your condition pairs, attach SE/TE label
  
  df_raw <- read_csv(csv_path, show_col_types = FALSE)
  stopifnot(all(c("name", "LA1_type", "cell_line", "condition") %in% names(df_raw)))
  metrics <- intersect(metrics, names(df_raw))
  
  df_raw <- df_raw %>%
    filter(cell_line == cell_line, LA1_type %in% c("SE", "TE")) %>%
    dplyr::rename(node_name = name, node_type = LA1_type) %>%
    mutate(node_type = as.character(node_type),
           condition = factor(condition, levels = condition_levels, ordered = TRUE))
  
  df_raw
  node_types <- df_raw %>% distinct(node_name, node_type) %>% dplyr::rename(node = node_name)
  
  cond_pairs <- combn(condition_levels, 2, simplify = FALSE)
  pci_all <- bind_rows(lapply(cond_pairs, function(p) partner_conservation(p[1], p[2]))) %>%
    inner_join(node_types, by = "node") %>%            # keep only SE/TE nodes
    mutate(transition = factor(transition,
                               levels = sapply(cond_pairs, function(p) paste0(p[1], " \u2192 ", p[2]))))
  
  
  star <- function(p) cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                          labels = c("***", "**", "*", "ns"))
  
  # SE vs TE test per transition
  pci_test <- pci_all %>%
    group_split(transition) %>%
    lapply(function(d) {
      d <- filter(d, !is.na(pci))
      w <- wilcox.test(pci ~ node_type, data = d)
      n1 <- sum(d$node_type == "SE"); n2 <- sum(d$node_type == "TE")
      data.frame(transition = d$transition[1], p = w$p.value,
                 r_rb = 1 - (2 * unname(w$statistic)) / (n1 * n2),
                 SE_med = median(d$pci[d$node_type == "SE"]),
                 TE_med = median(d$pci[d$node_type == "TE"]))
    }) %>% bind_rows() %>%
    mutate(p_adj = p.adjust(p, method = "BH"), sig = star(p))
  print(pci_test)
  
  ymax <- pci_all %>% group_by(transition) %>%
    summarise(y = max(pci, na.rm = TRUE), .groups = "drop") %>%
    left_join(pci_test, by = "transition")
  
  pci_all <- pci_all %>% mutate(across( c(node_type),
                                        ~ gsub("SE", "SEc", .x)))
  
  return(pci_all)
}