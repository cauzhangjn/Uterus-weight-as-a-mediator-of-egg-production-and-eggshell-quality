library(tidyverse)

df <- read.table("mediation_results_UW100_.txt",
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

df <- df %>% 
  filter(!grepl("_YP100_", SNP_PHENO))

df <- df %>%
  mutate(
    INT = sub(".*_(INT)$", "\\1", SNP_PHENO),
    Mediator = sub(".*_([^_]+)_INT$", "\\1", SNP_PHENO),
    SNP = sub("_(?:[^_]+)_INT$", "", SNP_PHENO)
  )

# 去掉 .100（如果你后面想统一命名）
df$Mediator <- gsub("\\.100", "", df$Mediator)


df$Mediator <- gsub("^EN85$", "EN85-100", df$Mediator)

df$SNP <- sub("_[^_]+$", "", df$SNP)




df_long <- df %>%
  select(SNP, Mediator, ACME, ADE, ACME_q_FDR, ADE_q_FDR) %>%
  pivot_longer(cols = c(ACME, ADE),
               names_to = "Effect_type",
               values_to = "Effect") %>%
  mutate(
    Direction = ifelse(Effect > 0, "Positive", "Negative"),
    Sig_class = case_when(
      Effect_type == "ACME" & ACME_q_FDR < 0.05 ~ "ACME_SIG",
      Effect_type == "ACME" & ACME_q_FDR >= 0.05 ~ "ACME_NS",
      Effect_type == "ADE" & ADE_q_FDR < 0.05 ~ "ADE_SIG",
      Effect_type == "ADE" & ADE_q_FDR >= 0.05 ~ "ADE_NS",
    )
  )

ggplot(df_long,
       aes(x = Effect,
           y = SNP,
           color = Direction,
           shape = Sig_class,
           fill  = Sig_class)) +
  geom_point(size = 3, stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Mediator, scales = "free_y") +
  scale_shape_manual(
    values = c(
      "ACME_SIG" = 21,
      "ACME_NS"  = 21,
      "ADE_SIG"  = 24,
      "ADE_NS"  =24
    )
  ) +
  scale_fill_manual(
    values = c(
      "ACME_SIG" = "black",
      "ACME_NS"  = "white",
      "ADE_SIG"  = "black",
      "ADE_NS"  = "white"
    )
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#D55E00",
      "Negative" = "#0072B5"
    )
  ) +
  labs(
    x = "Effect size",
    y = "SNP",
    shape = "Effect category",
    fill  = "Effect category",
    color = "Effect direction"
  ) +
  facet_wrap(~ Mediator, ncol = 4, scales = "free_y")+

  theme_classic(base_size = 12)


# 1️⃣ 创建 PDF 输出文件
pdf("mediation_plot_.pdf", width = 20, height = 10)  # 宽高可根据需要调整

# 2️⃣ 打印你的 ggplot 图
ggplot(df_long,
       aes(x = Effect,
           y = SNP,
           color = Direction,
           shape = Sig_class,
           fill  = Sig_class)) +
  geom_point(size = 3, stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Mediator, ncol = 4, scales = "free_y") +
  scale_shape_manual(
    values = c(
      "ACME_SIG" = 21,
      "ACME_NS"  = 21,
      "ADE_SIG"  = 24,
      "ADE_NS"   = 24
    )
  ) +
  scale_fill_manual(
    values = c(
      "ACME_SIG" = "black",
      "ACME_NS"  = "white",
      "ADE_SIG"  = "black",
      "ADE_NS"   = "white"
    )
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#D55E00",
      "Negative" = "#0072B5"
    )
  ) +
  labs(
    x = "Effect size",
    y = "SNP",
    shape = "Effect category",
    fill  = "Effect category",
    color = "Effect direction"
  ) +
  theme_classic(base_size = 12)

# 3️⃣ 关闭 PDF 设备
dev.off()
