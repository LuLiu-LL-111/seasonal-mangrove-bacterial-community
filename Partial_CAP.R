setwd("your_path")
library(vegan)
library(ape)
library(tidyverse)
library(ggplot2)

otu <- read.delim(
  "asv_table.txt",
  row.names = 1,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map <- read.delim(
  "group.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

otu <- t(otu)

otu_hell <- decostand(otu, method = "hellinger")
dis_bray <- vegdist(otu_hell, method = "bray")
p <- pcoa(dis_bray, correction = "lingoes")

site <- as.data.frame(p$vectors[, 1:2])
colnames(site) <- c("X1", "X2")

exp1 <- round(100 * p$values$Relative_eig[1], 2)
exp2 <- round(100 * p$values$Relative_eig[2], 2)

otu_cnt <- otu
otu_cnt <- otu_cnt[rownames(otu_hell), , drop = FALSE]

keep_rows <- rowSums(otu_cnt, na.rm = TRUE) > 0
keep_cols <- colSums(otu_cnt, na.rm = TRUE) > 0
otu_cnt <- otu_cnt[keep_rows, keep_cols, drop = FALSE]
otu_hell <- otu_hell[keep_rows, keep_cols, drop = FALSE]

W <- as.data.frame(p$vectors[, 1:4])
W <- W[rownames(otu_hell), , drop = FALSE]

badW <- apply(W, 1, function(z) any(is.na(z)))
if (any(badW)) {
  message("Some samples contain NA values in the PCoA coordinates; these samples were removed before species projection.")
  W <- W[!badW, , drop = FALSE]
  otu_hell <- otu_hell[rownames(W), , drop = FALSE]
  otu_cnt <- otu_cnt[rownames(W), , drop = FALSE]
}

zero_cols <- colSums(otu_hell, na.rm = TRUE) == 0
if (any(zero_cols)) {
  message(sprintf(
    "%d OTUs have zero abundance across the current sample set and were removed.",
    sum(zero_cols)
  ))
  otu_hell <- otu_hell[, !zero_cols, drop = FALSE]
}

species <- NULL
try({
  species <- vegan::wascores(otu_hell, W)
}, silent = TRUE)

if (is.null(species)) {
  X <- as.matrix(otu_hell)
  S <- as.matrix(W)
  wsum <- colSums(X, na.rm = TRUE)

  keep_sp <- wsum > 0
  if (any(!keep_sp)) {
    message(sprintf(
      "%d OTUs have zero effective abundance in the current dataset and were skipped.",
      sum(!keep_sp)
    ))
    X <- X[, keep_sp, drop = FALSE]
    wsum <- wsum[keep_sp]
  }

  num <- t(X) %*% S
  species <- sweep(num, 1, wsum, "/")
}

write.table(
  species,
  "pcoa_species_lingoes.txt",
  sep = "\t",
  col.names = NA,
  quote = FALSE
)

site2 <- as.data.frame(p$vectors[, 1:2])
site2 <- site2[rownames(W), , drop = FALSE]
colnames(site2) <- c("X1", "X2")

exp1 <- round(100 * p$values$Relative_eig[1], 2)
exp2 <- round(100 * p$values$Relative_eig[2], 2)

plot(
  site2$X1,
  site2$X2,
  pch = 19,
  col = "grey30",
  xlab = paste0("PCoA1 (", exp1, "%)"),
  ylab = paste0("PCoA2 (", exp2, "%)"),
  main = "PCoA (Bray + Lingoes)"
)
abline(h = 0, v = 0, col = "grey80")

library(tidyverse)
library(ggplot2)

map <- read.delim(
  "group.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(map) <- names(map) |>
  trimws() |>
  tolower() |>
  gsub("[^a-z0-9_]+", "_", x = _)

if ("otus" %in% names(map)) {
  names(map)[names(map) == "otus"] <- "sample"
}

map <- map |>
  mutate(
    season = factor(season, levels = sort(unique(season))),
    group = factor(group, levels = c("r", "b"))
  )

site2$sample <- rownames(site2)

df_plot <- site2 %>%
  left_join(map, by = "sample") %>%
  filter(!is.na(season), !is.na(group))

pal_season <- setNames(
  c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
  levels(df_plot$season)
)

shape_group <- c(r = 16, b = 17)

gg <- ggplot(df_plot, aes(X1, X2)) +
  geom_point(
    aes(color = season, shape = group),
    size = 3,
    alpha = 0.95
  ) +
  scale_color_manual(values = pal_season, name = "Season") +
  scale_shape_manual(values = shape_group, name = "Group") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_hline(yintercept = 0, color = "grey80") +
  labs(
    x = paste0("PCoA1 (", exp1, "%)"),
    y = paste0("PCoA2 (", exp2, "%)"),
    title = "PCoA: color = Season, shape = Group"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

print(gg)

gg +
  stat_ellipse(
    aes(color = season),
    level = 0.95,
    linewidth = 0.8,
    fill = NA,
    show.legend = FALSE
  )

library(vegan)

cap_season_partial <- capscale(
  otu_hell ~ season + Condition(site + group),
  data = map,
  distance = "bray"
)

anova(cap_season_partial)
anova(cap_season_partial, by = "terms")
anova(cap_season_partial, by = "axis")

sc_sites_sp <- scores(
  cap_season_partial,
  display = "sites",
  choices = 1:2
)

df_sp <- cbind(
  as.data.frame(sc_sites_sp),
  map[, c("season", "site", "group")]
)

eig_sp <- eigenvals(cap_season_partial, constrained = TRUE)
exp_sp <- round(100 * eig_sp[1:2] / sum(eig_sp), 2)

library(ggplot2)

gg_sp <- ggplot(df_sp, aes(CAP1, CAP2)) +
  geom_point(
    aes(color = season, shape = group),
    size = 3,
    alpha = 0.95
  ) +
  stat_ellipse(
    aes(color = season),
    level = 0.95,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_hline(yintercept = 0, color = "grey80") +
  labs(
    title = "Partial CAP (Bray): controlling for site and group, showing seasonal separation",
    x = paste0("CAP1 (", exp_sp[1], "% of constrained)"),
    y = paste0("CAP2 (", exp_sp[2], "% of constrained)")
  ) +
  theme_classic(base_size = 14)

print(gg_sp)
