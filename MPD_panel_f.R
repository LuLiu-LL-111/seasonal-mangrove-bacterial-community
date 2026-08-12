library(ape)
library(data.table)

setwd("your_path")

# ============================================================
# 1. Read input data
# ============================================================

grp <- fread("datanet1.txt", header = TRUE)
tr <- read.tree("root.tree")

required_cols <- c("OTU", "Group")
missing_cols <- setdiff(required_cols, names(grp))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "datanet1.txt is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

clean_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^'+|'+$", "", x)
  x <- gsub('^"+|"+$', "", x)
  x <- gsub("\r", "", x)
  x
}

grp[, OTU := clean_name(OTU)]
tr$tip.label <- clean_name(tr$tip.label)

if (anyDuplicated(grp$OTU) > 0) {
  dup <- unique(grp$OTU[duplicated(grp$OTU)])
  stop(
    paste0(
      "Duplicated OTU names were found in datanet1.txt: ",
      paste(head(dup, 20), collapse = ", ")
    )
  )
}

if (anyDuplicated(tr$tip.label) > 0) {
  dup <- unique(tr$tip.label[duplicated(tr$tip.label)])
  stop(
    paste0(
      "Duplicated tip labels were found in root.tree: ",
      paste(head(dup, 20), collapse = ", ")
    )
  )
}

grp[, Group := as.numeric(Group)]

if (any(is.na(grp$Group))) {
  stop("The Group column contains non-numeric or missing values.")
}

if (!all(grp$Group %in% c(0, 1))) {
  stop("The Group column must contain only 0 and 1.")
}

# ============================================================
# 2. Match OTUs between the group table and phylogenetic tree
# ============================================================

common <- intersect(grp$OTU, tr$tip.label)

if (length(common) == 0) {
  stop("No matching OTUs were found between datanet1.txt and root.tree.")
}

cat("Number of OTUs in datanet1.txt =", nrow(grp), "\n")
cat("Number of tips in root.tree =", length(tr$tip.label), "\n")
cat("Number of matched OTUs =", length(common), "\n")

tr2 <- drop.tip(
  tr,
  setdiff(tr$tip.label, common)
)

tips <- tr2$tip.label

grp2 <- grp[match(tips, grp$OTU)]

if (any(is.na(grp2$OTU))) {
  stop("Some tree tips could not be matched to datanet1.txt.")
}

stopifnot(all(grp2$OTU == tips))

# Calculate pairwise phylogenetic distances
D <- cophenetic(tr2)
D <- D[tips, tips, drop = FALSE]

# ============================================================
# 3. Define Group 1 and Group 0
# ============================================================

otu_g1 <- grp2$OTU[grp2$Group == 1]
otu_g0 <- grp2$OTU[grp2$Group == 0]

cat("Number of OTUs in Group 1 =", length(otu_g1), "\n")
cat("Number of OTUs in Group 0 =", length(otu_g0), "\n")

if (length(otu_g1) < 2) {
  stop("Group 1 must contain at least two OTUs to calculate MPD.")
}

if (length(otu_g0) < length(otu_g1)) {
  stop(
    paste0(
      "Group 0 contains only ",
      length(otu_g0),
      " OTUs, but ",
      length(otu_g1),
      " OTUs are required for each random sample."
    )
  )
}

# ============================================================
# 4. Define the MPD function
# ============================================================

calc_mpd <- function(otu_vec, D) {

  if (length(otu_vec) < 2) {
    return(NA_real_)
  }

  subD <- D[
    otu_vec,
    otu_vec,
    drop = FALSE
  ]

  mean(
    subD[lower.tri(subD, diag = FALSE)],
    na.rm = TRUE
  )
}

# Observed MPD of Group 1
obs_mpd <- calc_mpd(
  otu_g1,
  D
)

if (!is.finite(obs_mpd)) {
  stop("Observed MPD could not be calculated.")
}

# ============================================================
# 5. Generate the null MPD distribution
# ============================================================

set.seed(123)

nperm <- 9999
k <- length(otu_g1)

rand_mpd <- replicate(
  nperm,
  {
    sampled_otus <- sample(
      otu_g0,
      size = k,
      replace = FALSE
    )

    calc_mpd(
      sampled_otus,
      D
    )
  }
)

rand_mpd <- rand_mpd[
  is.finite(rand_mpd)
]

if (length(rand_mpd) == 0) {
  stop("No valid random MPD values were generated.")
}

# ============================================================
# 6. Calculate permutation-based statistics
# ============================================================

null_mean <- mean(
  rand_mpd,
  na.rm = TRUE
)

null_sd <- sd(
  rand_mpd,
  na.rm = TRUE
)

# One-sided test: Group 1 has a smaller MPD than expected
p_lower <- (
  sum(
    rand_mpd <= obs_mpd,
    na.rm = TRUE
  ) + 1
) / (
  length(rand_mpd) + 1
)

# One-sided test: Group 1 has a larger MPD than expected
p_upper <- (
  sum(
    rand_mpd >= obs_mpd,
    na.rm = TRUE
  ) + 1
) / (
  length(rand_mpd) + 1
)

# Two-sided test based on distance from the null mean
p_two <- (
  sum(
    abs(rand_mpd - null_mean) >=
      abs(obs_mpd - null_mean),
    na.rm = TRUE
  ) + 1
) / (
  length(rand_mpd) + 1
)

# Standardized effect size
ses_mpd <- (
  obs_mpd - null_mean
) / null_sd

cat("\n")
cat("Observed MPD =", obs_mpd, "\n")
cat("Null mean =", null_mean, "\n")
cat("Null SD =", null_sd, "\n")
cat("SES MPD =", ses_mpd, "\n")
cat("P lower (Group 1 smaller) =", p_lower, "\n")
cat("P upper (Group 1 larger) =", p_upper, "\n")
cat("Two-sided P =", p_two, "\n")

# ============================================================
# 7. Save numerical results
# ============================================================

mpd_summary <- data.frame(
  Observed_MPD = obs_mpd,
  Null_mean = null_mean,
  Null_SD = null_sd,
  SES_MPD = ses_mpd,
  P_lower = p_lower,
  P_upper = p_upper,
  P_two_sided = p_two,
  N_permutations = length(rand_mpd),
  Group1_size = length(otu_g1),
  Group0_size = length(otu_g0)
)

write.table(
  mpd_summary,
  file = "MPD_permutation_test_summary.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    Random_MPD = rand_mpd
  ),
  file = "MPD_null_distribution.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# 8. Plot the null MPD distribution
# ============================================================

plot_mpd_null <- function(rand_mpd, obs_mpd) {

  x_values <- c(
    rand_mpd,
    obs_mpd
  )

  x_range <- range(
    x_values,
    na.rm = TRUE
  )

  x_width <- diff(x_range)

  if (!is.finite(x_width) || x_width == 0) {
    x_width <- max(
      abs(x_range),
      1
    ) * 0.1
  }

  x_pad <- x_width * 0.05

  x_lim <- c(
    x_range[1] - x_pad,
    x_range[2] + x_pad
  )

  hist(
    rand_mpd,
    breaks = 40,
    probability = TRUE,
    col = "lightblue",
    border = "grey40",
    xlab = "MPD",
    ylab = "Density",
    main = "",
    xlim = x_lim
  )

  density_obj <- density(
    rand_mpd,
    na.rm = TRUE
  )

  lines(
    density_obj$x,
    density_obj$y,
    lwd = 1.5,
    col = "grey20"
  )

  # Draw a continuous baseline over the complete x-axis range
  segments(
    x0 = x_lim[1],
    y0 = 0,
    x1 = x_lim[2],
    y1 = 0,
    col = "grey40",
    lwd = 1
  )

  # Draw the observed MPD
  abline(
    v = obs_mpd,
    col = "red",
    lty = 2,
    lwd = 2
  )
}

# Display the figure
plot_mpd_null(
  rand_mpd,
  obs_mpd
)

# ============================================================
# 9. Save the figure
# ============================================================

pdf(
  "MPD_null_distribution.pdf",
  width = 5,
  height = 4
)

plot_mpd_null(
  rand_mpd,
  obs_mpd
)

dev.off()

png(
  "MPD_null_distribution.png",
  width = 1600,
  height = 1200,
  res = 300
)

plot_mpd_null(
  rand_mpd,
  obs_mpd
)

dev.off()

# ============================================================
# 10. Finish
# ============================================================

cat("\nAnalysis completed successfully.\n")
cat("Output files:\n")
cat("1. MPD_permutation_test_summary.txt\n")
cat("2. MPD_null_distribution.txt\n")
cat("3. MPD_null_distribution.pdf\n")
cat("4. MPD_null_distribution.png\n")
