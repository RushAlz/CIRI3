library(tidyverse)
library(edgeR)  # for cpm() normalization

setwd("/pastel/projects/PARDoS/circrna_pardos")

# ---------------------------------------------------------------
# Load decoupled BSJ/FSJ results
# ---------------------------------------------------------------
decoupled_res_df  <- read_tsv("FINALIZE/test_10")
decoupled_BSJ_res <- read_tsv("FINALIZE/test_10.BSJ_Matrix")
decoupled_FSJ_res <- read_tsv("FINALIZE/test_10.FSJ_Matrix")

# Sanity check: same circRNA_ID order in both matrices
stopifnot(identical(decoupled_BSJ_res$circRNA_ID, decoupled_FSJ_res$circRNA_ID))

# Convert to numeric matrices with rownames
BSJ_mat <- as.matrix(decoupled_BSJ_res[,-1])
FSJ_mat <- as.matrix(decoupled_FSJ_res[,-1])
rownames(BSJ_mat) <- decoupled_BSJ_res$circRNA_ID
rownames(FSJ_mat) <- decoupled_FSJ_res$circRNA_ID

# ---------------------------------------------------------------
# Filter low-count circRNAs
# Keep circRNAs with BSJ >= 2 in at least 50% of samples
# (matches edgeR::filterByExpr philosophy but simpler/explicit)
# ---------------------------------------------------------------
n_samples    <- ncol(BSJ_mat)
min_samples  <- ceiling(n_samples * 0.5)
keep_circRNA <- rowSums(BSJ_mat >= 2) >= min_samples

message(sprintf("Keeping %d / %d circRNAs (>=4 BSJ in >=%d/%d samples)",
                sum(keep_circRNA), length(keep_circRNA), min_samples, n_samples))

BSJ_filt <- BSJ_mat[keep_circRNA, ]
FSJ_filt <- FSJ_mat[keep_circRNA, ]

# ---------------------------------------------------------------
# Normalization
#   - BSJ and FSJ counts: CPM with TMM library-size normalization,
#     then log2(CPM + 1) for distribution plotting
#   - Ratio: log2((BSJ + 1) / (FSJ + 1))  -> avoids Inf/NaN from zeros
# ---------------------------------------------------------------
# FSJ: TMM-normalized log2 CPM
dge_FSJ <- DGEList(counts = FSJ_filt)
dge_FSJ <- calcNormFactors(dge_FSJ, method = "TMM")
FSJ_logcpm <- cpm(dge_FSJ, log = TRUE, prior.count = 1)

# BSJ: TMM-normalized log2 CPM
dge_BSJ <- DGEList(counts = BSJ_filt)
dge_BSJ <- calcNormFactors(dge_BSJ, method = "TMM")
BSJ_logcpm <- cpm(dge_BSJ, log = TRUE, prior.count = 1)

# Override the auto-computed lib.size with the linear-derived one
lib_size_from_linear <- dge_FSJ$samples$lib.size * dge_FSJ$samples$norm.factors
names(lib_size_from_linear) <- rownames(dge_FSJ$samples)

dge_BSJ <- DGEList(counts = BSJ_filt)
dge_BSJ$samples$lib.size    <- lib_size_from_linear
dge_BSJ$samples$norm.factors <- 1   # already baked into lib.size above
BSJ_logcpm <- cpm(dge_BSJ, log = TRUE, prior.count = 1)

# Circularization ratio (log2): pseudocount avoids 0/0 and log(0)
ratio_log2 <- log2((BSJ_filt + 1) / (FSJ_filt + 1))

# ---------------------------------------------------------------
# Helper to convert matrix -> long tibble for ggplot
# ---------------------------------------------------------------
mat_to_long <- function(mat, value_name) {
  as.data.frame(mat) %>%
    rownames_to_column("circRNA_ID") %>%
    pivot_longer(cols = -circRNA_ID, names_to = "Sample", values_to = value_name)
}

# Hide legend if too many samples (typical for PARDoS)
sample_legend <- if (n_samples > 20) "none" else "right"

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
p_bsj <- mat_to_long(BSJ_logcpm, "logCPM") %>%
  ggplot(aes(x = logCPM, group = Sample, color = Sample)) +
  geom_density(alpha = 0.4) +
  labs(x = expression(log[2]*"(CPM)"),
       title = "BSJ distribution (TMM-normalized, filtered)") +
  theme_bw() + theme(legend.position = sample_legend)

p_fsj <- mat_to_long(FSJ_logcpm, "logCPM") %>%
  ggplot(aes(x = logCPM, group = Sample, color = Sample)) +
  geom_density(alpha = 0.4) +
  labs(x = expression(log[2]*"(CPM)"),
       title = "FSJ distribution (TMM-normalized, filtered)") +
  theme_bw() + theme(legend.position = sample_legend)

p_ratio <- mat_to_long(ratio_log2, "log2_ratio") %>%
  ggplot(aes(x = log2_ratio, group = Sample, color = Sample)) +
  geom_density(alpha = 0.4) +
  labs(x = expression(log[2]*"((BSJ+1)/(FSJ+1))"),
       title = "Circularization ratio distribution (filtered)") +
  theme_bw() + theme(legend.position = sample_legend)

print(p_bsj)
print(p_fsj)
print(p_ratio)

# ---------------------------------------------------------------
# Compare with original single-BAM results
# ---------------------------------------------------------------
decoupled_BSJ_res <- read_tsv("FINALIZE/test_10.BSJ_Matrix")
decoupled_FSJ_res <- read_tsv("FINALIZE/test_10.FSJ_Matrix")
res_BSJ_merged_long = decoupled_BSJ_res %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample")) 
res_FSJ_merged_long = decoupled_FSJ_res %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample"))

# Compare PARDOS results with Jishu results
pardos_res = list.files("OLD", pattern = ".CIRI3", full.names = T)
names(pardos_res) = basename(pardos_res)
pardos_res_df = purrr::map_df(pardos_res, ~ read_tsv(.x, show_col_types = FALSE), .id = "file") %>%
  mutate(Sample = gsub("(.*)\\.CIRI3", "\\1", basename(file))) %>%
  select(-file)

table(pardos_res_df$Sample)

pardos_BSJ_check = pardos_res_df %>% select(circRNA_ID, Sample, value = `#junction_reads`) 
pardos_FSJ_check = pardos_res_df %>% select(circRNA_ID, Sample, value = `#non_junction_reads`)

res_BSJ_merged_long %>% 
  full_join(pardos_BSJ_check, by = c("circRNA_ID", "Sample"), suffix = c(".decoupled", ".old")) %>%
  ggplot(aes(x = value.old, y = value.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Single-BAM CIRI3 (original)", y = "Decoupled Analysis BSJ Count") +
  theme_bw() +
  facet_wrap(~ Sample, scales = "free", nrow = 2)

res_FSJ_merged_long %>% 
  full_join(pardos_FSJ_check, by = c("circRNA_ID", "Sample"), suffix = c(".decoupled", ".old")) %>%
  ggplot(aes(x = value.old, y = value.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Single-BAM CIRI3 (original)", y = "Decoupled Analysis FSJ Count") +
  theme_bw() +
  facet_wrap(~ Sample, scales = "free", nrow = 2)

# ---------------------------------------------------------------
# Compare with original multi-BAM results
# ---------------------------------------------------------------
joint_res_df = read_tsv("/pastel/Github_scripts/CIRI3/test_decoupled_manual/multi_bam/CIRI3_joint_results.txt")
joint_BSJ_res = read_tsv("/pastel/Github_scripts/CIRI3/test_decoupled_manual/multi_bam/CIRI3_joint_results.txt.BSJ_Matrix")
joint_FSJ_res = read_tsv("/pastel/Github_scripts/CIRI3/test_decoupled_manual/multi_bam/CIRI3_joint_results.txt.FSJ_Matrix")
joint_BSJ_res_long = joint_BSJ_res %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample"))
joint_BSJ_res_long$Sample = gsub("(.*)_S(.*)","\\1",
                                 gsub("STAR_output_","",basename(dirname(joint_BSJ_res_long$Sample))))
joint_FSJ_res_long = joint_FSJ_res %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample"))
joint_FSJ_res_long$Sample = gsub("(.*)_S(.*)","\\1",
                                 gsub("STAR_output_","",basename(dirname(joint_FSJ_res_long$Sample))))

res_BSJ_merged_long %>%
  full_join(joint_BSJ_res_long, by = c("circRNA_ID", "Sample"), suffix = c(".decoupled", ".joint")) %>%
  ggplot(aes(x = value.joint, y = value.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Joint Analysis BSJ Count", y = "Decoupled Analysis BSJ Count") +
  theme_bw() +
  facet_wrap(~ Sample, scales = "free", nrow = 2)

res_FSJ_merged_long %>%
  full_join(joint_FSJ_res_long, by = c("circRNA_ID", "Sample"), suffix = c(".decoupled", ".joint")) %>%
  ggplot(aes(x = value.joint, y = value.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Joint Analysis FSJ Count", y = "Decoupled Analysis FSJ Count") +
  theme_bw() +
  facet_wrap(~ Sample, scales = "free", nrow = 2)