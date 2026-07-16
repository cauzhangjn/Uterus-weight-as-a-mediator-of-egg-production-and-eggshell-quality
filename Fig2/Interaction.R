library(tidyverse)
library(igraph)
library(ggraph)


df <- read.table(
  "互作.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 统一 RNA-seq 列名（防止有点号）
colnames(df)[colnames(df) == "RNA.seq"] <- "RNA_seq"



edges <- bind_rows(
  
  # GWAS → Gene
  df %>%
    filter(!is.na(GWAS)) %>%
    transmute(
      from = GWAS,
      to   = Gene,
      source = "GWAS"
    ),
  
  # TWAS_SMR → Gene（NA 不画）
  df %>%
    filter(!is.na(TWAS_SMR)) %>%
    transmute(
      from = TWAS_SMR,
      to   = Gene,
      source = "TWAS_SMR"
    )
)



nodes <- tibble(
  name = unique(c(edges$from, edges$to))
) %>%
  mutate(
    
    # 节点类型：只区分 SNP 来源
    node_class = case_when(
      name %in% df$GWAS ~ "GWAS",
      name %in% df$TWAS_SMR ~ "TWAS_SMR",
      TRUE ~ "Gene"
    ),
    
    # RNA-seq 信息（只对 Gene 有）
    RNA_seq = df$RNA_seq[match(name, df$Gene)],
    
    # Gene 填充色分组：调控方向
    gene_fill = case_when(
      RNA_seq == "Up" ~ "Up",
      RNA_seq == "Down" ~ "Down",
      TRUE ~ "No_Significant"
    ),
    
    # 标签规则
    label = case_when(
      node_class %in% c("GWAS", "TWAS_SMR") ~ name,
      RNA_seq %in% c("Up", "Down") ~ name,
      TRUE ~ NA_character_
    )
  )



g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = FALSE
)



ggraph(g, layout = "stress") +
  
  # 边：GWAS 实线，TWAS_SMR 虚线
  geom_edge_link(
    aes(linetype = source),
    color = "grey75",
    width = 0.6,
    alpha = 0.8
  ) +
  
  # SNP 节点（空心）
  geom_node_point(
    data = function(x) subset(x, node_class != "Gene"),
    aes(color = node_class),
    size = 6.5,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  
  # Gene 节点（实心，按调控方向）
  geom_node_point(
    data = function(x) subset(x, node_class == "Gene"),
    aes(fill = gene_fill),
    size = 3.8,
    shape = 21,
    color = "grey40",
    stroke = 0.5
  ) +
  
  # 标签（全部黑色，字体偏大）
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 4,
    color = "black",
    na.rm = TRUE
  ) +
  
  # SNP 颜色
  scale_color_manual(
    values = c(
      "GWAS" = "#7F0000",
      "TWAS_SMR" = "#08306B"
    )
  ) +
  
  # Gene 填充色（调控方向）
  scale_fill_manual(
    values = c(
      "Up" = "#F4A582",
      "Down" = "#92C5DE",
      "No_Significant" = "#D9D9D9"
    )
  ) +
  
  # 线型
  scale_linetype_manual(
    values = c(
      "GWAS" = "solid",
      "TWAS_SMR" = "dashed"
    )
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    text = element_text(family = "Arial"),
    plot.margin = margin(10, 10, 10, 10)
  )


ggsave(
  "GWAS_TWAS_RNAseq_network.pdf",
  width = 18,
  height = 6,
  device = cairo_pdf
)
