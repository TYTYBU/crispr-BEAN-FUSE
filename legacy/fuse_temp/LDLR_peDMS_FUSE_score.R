# =============================================================================
# LDLR Prime Editing DMS Dataset — FUSE Score Calculation
# =============================================================================
# Standalone script extracted from get_LDLR_DMS_FUSE_score.R
# All required functions from AA_proj_functions_101421.R are integrated below.
#
# Required input files (copy these into the folders indicated):
#   FUNSUM/funsum_combined_121324.rds
#   FUNSUM/funsum_SS_combined_121324.rds
#   FUNSUM/funsum_combined_noLOF_121324.rds      (only if include_LOF = FALSE)
#   FUNSUM/funsum_SS_combined_noLOF_121324.rds   (only if include_LOF = FALSE)
#   functional_datasets/0725_HCT116_LDLRPE_finalBEAN_FUSE_e.xlsx
#   dssp_out/LDLR.dss
#
# Output:
#   output/LDLR_peDMS_FUSE_score_071725.csv
# =============================================================================

# ── Working directory ─────────────────────────────────────────────────────────
# Set this to the folder containing the subfolders above.
# If you run the script from its own location, leave it as-is.
setwd("/Users/tianyu/Downloads/fuse_temp")

# ── Libraries ─────────────────────────────────────────────────────────────────
library(openxlsx)
library(tidyverse)
library(ptm)
library(ggridges)   # required by draw_score_distribution_ggplot


# =============================================================================
# Helper functions (sourced from AA_proj_functions_101421.R)
# =============================================================================

# ── James-Stein shrinkage estimator ──────────────────────────────────────────
get_js <- function(m){
  mbar  <- mean(colMeans(m, na.rm = TRUE), na.rm = TRUE)   # global mean
  mu0   <- colMeans(m, na.rm = TRUE)                        # column means
  s2    <- var(as.vector(m), na.rm = TRUE) /
           (nrow(m) * ncol(m) - length(which(is.na(as.vector(m)))))
  cval  <- 1 - (ncol(m) - 2) * s2 / sum((mu0 - mbar)^2, na.rm = TRUE)
  adj_val <- cval * (mu0 - mbar)
  js_est  <- mbar + ifelse(is.na(adj_val), 0, adj_val)
  return(js_est)
}

# ── Convert FUNSUM matrix to substitution pair table ─────────────────────────
funsum_to_subTable <- function(df.funsum){
  df.sub_tb <- c()
  for (i in 1:nrow(df.funsum)){
    temp <- data.frame(
      aa_pair = paste0(rownames(df.funsum)[i], colnames(df.funsum)),
      score   = as.vector(df.funsum[i, ])
    )
    df.sub_tb <- rbind(df.sub_tb, temp)
  }
  return(df.sub_tb)
}

# ── Annotate variants with functional class (MIS / LOF / SYN) ────────────────
get_functional_class <- function(df){
  df$functional_class <- NA
  ind <- which(df$aaref != df$aaalt & df$aaalt != '*')
  df$functional_class[ind] <- "MIS"
  ind <- which(df$aaref != "*" & df$aaalt == "*")
  df$functional_class[ind] <- "LOF"
  ind <- which(df$aaref == df$aaalt)
  df$functional_class[ind] <- "SYN"
  return(df)
}

# ── Plot raw / normalised score distributions by functional class ─────────────
draw_score_distribution_ggplot <- function(df, score_type = "raw_score", bin_ct = 30){
  require(ggridges)
  plt <- df %>%
    ggplot(aes(x = .data[[score_type]], y = functional_class, fill = functional_class)) +
    geom_density_ridges(alpha = 0.6, stat = "binline", bins = bin_ct) +
    ylab("") + xlab(score_type) +
    theme_ridges()
  return(plt)
}

# ── Normalise a scoreset to a [0, 1]-anchored scale ──────────────────────────
normalize_scoreset <- function(df, lower_bound = 0.05, upper_bound = 0.95,
                               force_quantile = FALSE){
  if (!"functional_class" %in% colnames(df)){
    df$functional_class <- NA
    ind <- which(df$aaref != df$aaalt & df$aaalt != '*')
    df$functional_class[ind] <- "MIS"
    ind <- which(df$aaref != "*" & df$aaalt == "*")
    df$functional_class[ind] <- "LOF"
    ind <- which(df$aaref == df$aaalt)
    df$functional_class[ind] <- "SYN"
  }

  ind_lof <- which(df$functional_class == "LOF")
  median_lof <- if (length(ind_lof) > 10) median(df$raw_score[ind_lof], na.rm = TRUE) else NA

  ind_mis    <- which(df$functional_class == "MIS")
  median_mis <- median(df$raw_score[ind_mis], na.rm = TRUE)

  ind_syn    <- which(df$functional_class == "SYN")
  median_syn <- if (length(ind_syn) > 10) median(df$raw_score[ind_syn], na.rm = TRUE) else NA

  temp <- df$raw_score

  if (!is.na(median_syn) & !is.na(median_lof)){
    if (median_lof < median_syn){
      temp       <- -temp
      median_lof <- -median_lof
      median_syn <- -median_syn
    }
    df$norm_raw_score <- (temp - median_syn) / (median_lof - median_syn)

  } else {
    mis_scores <- df$raw_score[ind_mis]

    if (!is.na(median_syn) & is.na(median_lof)){
      if (median_mis < median_syn){
        temp       <- -temp
        median_syn <- -median_syn
        mis_scores <- -mis_scores
      }
      median_lof        <- quantile(mis_scores, upper_bound, na.rm = TRUE)
      df$norm_raw_score <- (temp - median_syn) / (median_lof - median_syn)

    } else if (is.na(median_syn) & !is.na(median_lof)){
      if (median_mis > median_lof){
        temp       <- -temp
        median_lof <- -median_lof
        mis_scores <- -mis_scores
      }
      median_syn        <- quantile(mis_scores, lower_bound, na.rm = TRUE)
      df$norm_raw_score <- (temp - median_syn) / (median_lof - median_syn)

    } else {
      ind_p    <- which(df$aaalt == "P")
      median_p <- median(df$raw_score[ind_p], na.rm = TRUE)
      if (median_p < median_mis){
        temp       <- -temp
        median_syn <- -median_syn
        mis_scores <- -mis_scores
      }
      median_syn        <- quantile(mis_scores, lower_bound, na.rm = TRUE)
      median_lof        <- quantile(mis_scores, upper_bound, na.rm = TRUE)
      df$norm_raw_score <- (temp - median_syn) / (median_lof - median_syn)
    }
  }

  if (force_quantile){
    ind_mis    <- which(df$functional_class == "MIS")
    temp       <- df$norm_raw_score
    mis_scores <- df$norm_raw_score[ind_mis]
    median_syn        <- quantile(mis_scores, lower_bound,  na.rm = TRUE)
    median_lof        <- quantile(mis_scores, upper_bound, na.rm = TRUE)
    df$norm_raw_score <- (temp - median_syn) / (median_lof - median_syn)
  }

  return(df)
}

# ── FUSE score calculation for a single gene (with per-SS FUNSUM) ─────────────
de_noise_ss_1gene <- function(df, pos_mean_method, df.funsum, ls.funsum_ss,
                               dss_path, include_LOF = TRUE, show_func_class = FALSE){
  # Filter rows without amino acid info
  ind <- which(!is.na(df$aaref) & !is.na(df$aaalt))
  df  <- df[ind, ]

  if (include_LOF){
    aa_list <- unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list <- unlist(strsplit("RHKDESTNQCGPAVILMFYW",  split = ""))
    ind <- which((df$aaref != "*") & (df$aaalt != "*"))
    df  <- df[ind, ]
  }

  df.sub_tb <- funsum_to_subTable(df.funsum)

  if (is.null(df[["gene"]])) df[["gene"]] <- "gene"
  df$gene_aa_str <- paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)

  # Collapse duplicate substitutions
  df <- df %>%
    group_by(gene_aa_str) %>%
    summarise(
      gene          = unique(gene),
      aapos         = unique(aapos),
      aaref         = unique(aaref),
      aaalt         = unique(aaalt),
      raw_score     = mean(raw_score),
      norm_raw_score = mean(norm_raw_score)
    )

  # ── Positional component ───────────────────────────────────────────────────
  df.pos_score <- df %>%
    group_by(gene, aapos) %>%
    summarise(aaref = unique(aaref), pos_mean = NA)
  df.pos_score$gene_aapos <- paste0(df.pos_score$gene, "---", df.pos_score$aapos)

  temp <- matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list

  if (pos_mean_method == "funsum"){
    temp2 <- matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
    colnames(temp2) <- aa_list
  }

  for (i in 1:nrow(df.pos_score)){
    ind  <- which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 <- df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] <- df$norm_raw_score[ind[ind2]]

    if (pos_mean_method == "funsum"){
      ind <- which(!is.na(temp[i, ]))
      temp2[i, ind] <- temp[i, ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }
  }

  if      (pos_mean_method == "mean")   df.pos_score$pos_mean <- rowMeans(temp, na.rm = TRUE)
  else if (pos_mean_method == "median") df.pos_score$pos_mean <- apply(temp, 1, FUN = median, na.rm = TRUE)
  else if (pos_mean_method == "js")     df.pos_score$pos_mean <- get_js(t(temp))
  else if (pos_mean_method == "funsum") df.pos_score$pos_mean <- get_js(t(temp2))

  # ── Build output table with all possible substitutions ────────────────────
  df.out <- df.pos_score %>%
    select(gene, aapos, aaref) %>%
    dplyr::slice(rep(1:n(), each = length(aa_list)))
  df.out$aaalt <- rep(aa_list, nrow(df.pos_score))

  if (show_func_class){
    df.out$functional_class <- NA
    ind <- which(df.out$aaref != df.out$aaalt & df.out$aaalt != '*')
    df.out$functional_class[ind] <- "MIS"
    ind <- which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] <- "LOF"
    ind <- which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] <- "SYN"
  }

  df.out$gene_aa_str <- paste0(df.out$gene, "---", df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$gene_aapos  <- paste0(df.out$gene, "---", df.out$aapos)
  df.out$aa_pair     <- paste0(df.out$aaref, df.out$aaalt)

  # ── DSSP secondary structure annotation ───────────────────────────────────
  require(ptm)
  df.out$ss  <- NA
  df.out$acc <- NA
  df.dssp    <- parse.dssp(file = dss_path, keepfiles = TRUE)
  ind        <- match(df.out$aapos, table = df.dssp$respdb)
  df.out$ss  <- df.dssp$ss[ind]
  df.out$acc <- df.dssp$sasa[ind]

  # ── Assign scores ──────────────────────────────────────────────────────────
  ind                <- match(df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$raw_score   <- df$raw_score[ind]
  df.out$norm_raw_score <- df$norm_raw_score[ind]

  ind                <- match(df.out$gene_aapos, table = df.pos_score$gene_aapos)
  df.out$pos_score   <- df.pos_score$pos_mean[ind]

  ind                <- match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score   <- df.sub_tb$score[ind]

  ss_list <- c("G" = "Helices", "H" = "Helices", "I" = "Helices",
               "E" = "Strands", "B" = "Strands",
               "T" = "Loops",   "S" = "Loops",   "C" = "Loops")
  df.out$sub_score_ss <- NA

  for (ss in names(ss_list)){
    ind <- which(df.out$ss == ss)
    if (length(ind)){
      df.sub_tb_ss <- funsum_to_subTable(ls.funsum_ss[[ss]])
      ind2 <- match(df.out$aa_pair[ind], table = df.sub_tb_ss$aa_pair)
      df.out$sub_score_ss[ind] <- df.sub_tb_ss$score[ind2]
    }
  }

  df.out$final_score    <- df.out$pos_score + df.out$sub_score
  df.out$final_score_ss <- df.out$pos_score + df.out$sub_score_ss

  ind <- which(is.na(df.out$final_score_ss))
  df.out$final_score_ss[ind] <- df.out$final_score[ind]
  df.out$sub_score_ss[ind]   <- df.out$sub_score[ind]

  return(df.out)
}


# =============================================================================
# LDLR prime editing DMS dataset
# =============================================================================

gene_id     <- "LDLR"
include_LOF <- TRUE

# ── Load FUNSUM matrices ───────────────────────────────────────────────────────
if (include_LOF){
  df.funsum_all <- readRDS("./FUNSUM/funsum_combined_121324.rds")
  ls.funsum_ss  <- readRDS("./FUNSUM/funsum_SS_combined_121324.rds")
} else {
  df.funsum_all <- readRDS("./FUNSUM/funsum_combined_noLOF_121324.rds")
  ls.funsum_ss  <- readRDS("./FUNSUM/funsum_SS_combined_noLOF_121324.rds")
}

# ── Load raw DMS data ──────────────────────────────────────────────────────────
df.raw <- read.xlsx(
  "./functional_datasets/0725_HCT116_LDLRPE_finalBEAN_FUSE_e.xlsx",
  sheet       = "data",
  check.names = TRUE
)

score_col      <- "bean_element_result.MixtureNormal_HCT116_rep7_normalization0702_muZ_adj_SD.1"
score_col_nice <- "MixtureNormal_HCT116_rep7_normalization0702_muZadj_SD_lt1"

df <- df.raw %>%
  filter(`Missense.` == TRUE | `Stop.` == TRUE | `Synonymous.` == TRUE) %>%
  mutate(
    gene      = gene_id,
    aaref     = Starting.AA,
    aaalt     = Final.AA,
    aapos     = as.numeric(LDLR.pos),
    raw_score = as.numeric(.[[score_col]])
  ) %>%
  select(gene, aapos, aaref, aaalt, raw_score) %>%
  filter(!is.na(raw_score))

df$aaalt[df$aaalt == "Z"] <- "*"

# Make observations unique (average duplicates)
df <- df %>%
  group_by(gene, aapos, aaref, aaalt) %>%
  summarise(raw_score = mean(raw_score))

# ── Annotate functional class ─────────────────────────────────────────────────
df <- get_functional_class(df)

# ── Normalise scores ──────────────────────────────────────────────────────────
df <- normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = TRUE)
draw_score_distribution_ggplot(df, score_type = "norm_raw_score")

# ── Calculate FUSE score ──────────────────────────────────────────────────────
df.out <- de_noise_ss_1gene(
  df             = df,
  pos_mean_method = "js",
  df.funsum      = df.funsum_all,
  ls.funsum_ss   = ls.funsum_ss,
  dss_path       = "./dssp_out/LDLR.dss",
  include_LOF    = include_LOF,
  show_func_class = TRUE
)

df.fuse <- df.out %>%
  mutate(
    FUSE_score    = final_score,
    FUSE_SS_score = final_score_ss
  ) %>%
  select(gene, aapos, aaref, aaalt, functional_class,
         raw_score, norm_raw_score, FUSE_score, FUSE_SS_score) %>%
  mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt))

# ── Write output ──────────────────────────────────────────────────────────────
write_csv(df.fuse, file = "./output/LDLR_peDMS_FUSE_score_071725.csv")
message("Done! Output saved to: ./output/LDLR_peDMS_FUSE_score_071725.csv")
