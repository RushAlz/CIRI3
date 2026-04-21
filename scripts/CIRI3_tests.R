# ============================================================================
# CIRI3 Test Data Analysis (Shell Script)
# ============================================================================
cd /pastel/tools/circRNA_tools/test_data
THREADS=64
REF_FASTA=GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=gencode.v32.primary_assembly.annotation.gtf
GENOMEDIR=${PWD}

# git clone https://github.com/gyjames/CIRI3.git
CIRI3_JAR_PATH=/pastel/tools/circRNA_tools/CIRI3/CIRI3_Java_18.0.1.jar

# List of sample prefixes
SAMPLES=(
  "Div_100_S91"
  "Div_101_S92"
  "PARDOS_1_S1"
  "PARDOS_2_S2"
)

# BWA
# conda create -n bwa_env bwa
conda activate bwa_env

# Index reference
bwa index -a bwtsw ${REF_FASTA}

for SAMPLE in "${SAMPLES[@]}"; do
  R1="${SAMPLE}_R1_001.fastq.gz"
  R2="${SAMPLE}_R2_001.fastq.gz"
  bwa mem -t ${THREADS} -T 19 ${REF_FASTA} ${R1} ${R2} > ${SAMPLE}.sam
done

# git clone https://github.com/gyjames/CIRI3.git
# cd CIRI3
# conda env create -n CIRI3 -f ./environment.yaml
conda activate CIRI3

# STAR
# Index reference
mkdir -p ${GENOMEDIR}/STAR
STAR --runThreadN ${THREADS} \
--runMode genomeGenerate \
--genomeDir ${GENOMEDIR}/STAR \
--genomeFastaFiles \
${REF_FASTA} \
--sjdbGTFfile ${GTF_FILE}

for SAMPLE in "${SAMPLES[@]}"; do
  R1="${SAMPLE}_R1_001.fastq.gz"
  R2="${SAMPLE}_R2_001.fastq.gz"
  
  STAR --runThreadN ${THREADS} \
  --genomeDir ${GENOMEDIR}/STAR \
  --outSAMtype SAM \
  --readFilesIn ${R1} ${R2} \
  --readFilesCommand zcat \
  --outFileNamePrefix STAR_output_${SAMPLE}/ \
  --outReadsUnmapped Fastx \
  --outSJfilterOverhangMin 15 12 12 12 \
  --alignSJoverhangMin 15 \
  --alignSJDBoverhangMin 15 \
  --outFilterMultimapNmax 20 \
  --outFilterScoreMin 1 \
  --outFilterMatchNmin 1 \
  --outFilterMismatchNmax 2 \
  --chimSegmentMin 15 \
  --chimScoreMin 15 \
  --chimJunctionOverhangMin 15
done

conda deactivate
conda activate bwa_env

for SAMPLE in "${SAMPLES[@]}"; do
  bwa mem -t ${THREADS} \
  -T 19 \
  ${REF_FASTA} \
  STAR_output_${SAMPLE}/Unmapped.out.mate1 \
  STAR_output_${SAMPLE}/Unmapped.out.mate2 > STAR_output_${SAMPLE}/bwa.sam
done

conda deactivate
conda activate CIRI3

# Processing CIRI3 by sample
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${PWD}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${PWD}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
  -O ${SAMPLE}.CIRI3results.txt \
  -F ${REF_FASTA} \
  -A ${GTF_FILE} \
  -Ma 1 \
  -W 0 \
  -T ${THREADS}
done

# Joint processing
cat /dev/null > my_samples.tsv
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=STAR_output_${SAMPLE}/bwa.sam
  
  echo "${ChimericOutJunction},${AlignedOutSam},${BWASam}" >> my_samples.tsv
done

$CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} -I my_samples.tsv \
-O CIRI3_joint_results.txt \
-F ${REF_FASTA} \
-W 1 \
-Ma 1 \
-T ${THREADS}

# Decoupled processing
DECOUPLED_JAR="/pastel/Github_scripts/CIRI3/CIRI3_decoupled.jar"

# SCAN1
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${PWD}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${PWD}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} SCAN1 \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
  -O ${SAMPLE}.SCAN1 \
  -F ${REF_FASTA} \
  -A ${GTF_FILE} \
  -Ma 1 \
  -W 0 \
  -T ${THREADS}
done

cat /dev/null > samples_scan1.tsv
for SAMPLE in "${SAMPLES[@]}"; do
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  echo "${BWASam}\t${SAMPLE}.SCAN1.scan1_meta" >> samples_scan1.tsv
done

# BUILD_UNIVERSE
$CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} BUILD_UNIVERSE \
-I  samples_scan1.tsv \
-F  ${REF_FASTA} \
-O  decoupled_test

# SCAN2
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${PWD}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${PWD}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} SCAN2 \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
  -CU decoupled_test.universe \
  -O ${SAMPLE}.SCAN2 \
  -F ${REF_FASTA} \
  -Ma 1 \
  -T ${THREADS}
done

# FINALIZE
cat /dev/null > finalize_samples.tsv
# /abs/path/sample1.sam   /abs/path/sample1.fsj_counts   8   sample1
# /abs/path/sample2.sam   /abs/path/sample2.fsj_counts   8   sample2

for SAMPLE in "${SAMPLES[@]}"; do
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  FSJ_counts=${PWD}/${SAMPLE}.SCAN2.fsj_counts
  echo "${BWASam}\t${FSJ_counts}\t8\t${SAMPLE}" >> finalize_samples.tsv
done


$CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} FINALIZE \
-I  samples_scan1.tsv \
-F  ${REF_FASTA} \
-O  decoupled_test \
-A ${GTF_FILE}

# ============================================================================
# CIRI3 Test Data Analysis (R Script - Downstream)
# ============================================================================

library(tidyverse)

setwd("/pastel/tools/circRNA_tools/test_data")
res1 = read_tsv("Div_100_S91.CIRI3results.txt", show_col_types = FALSE)
res2 = read_tsv("Div_101_S92.CIRI3results.txt", show_col_types = FALSE)
res = bind_rows(res1 %>% mutate(Sample = "Div_100_S91"),
                res2 %>% mutate(Sample = "Div_101_S92")) 

# Make Score Matrix; rownames are circRNA_ID, columns are Score (sample names)
BSJ_matrix = res %>%
  select(circRNA_ID, Sample, Score) %>%
  pivot_wider(names_from = Sample, values_from = Score) %>%
  column_to_rownames(var = "circRNA_ID") %>%
  as.matrix() %>% as.data.frame()
cor(BSJ_matrix, use = "complete.obs")

# Make #non_junction_reads matrix (rename as FSJ)
FSJ_matrix = res %>%
  select(circRNA_ID, Sample, `#non_junction_reads`) %>%
  pivot_wider(names_from = Sample, values_from = `#non_junction_reads`) %>%
  column_to_rownames(var = "circRNA_ID") %>%
  as.matrix() %>% as.data.frame()
cor(FSJ_matrix, use = "complete.obs")

# Make ration matrix (junction_reads_ratio)
ratio_matrix = res %>%
  select(circRNA_ID, Sample, junction_reads_ratio) %>%
  pivot_wider(names_from = Sample, values_from = junction_reads_ratio) %>%
  column_to_rownames(var = "circRNA_ID") %>%
  as.matrix() %>% as.data.frame()
cor(ratio_matrix, use = "complete.obs")

# Compare results from joint analysis 
joint_BSJ_res = read_tsv("CIRI3_joint_results.txt.BSJ_Matrix")
joint_FSJ_res = read_tsv("CIRI3_joint_results.txt.FSJ_Matrix")

joint_BSJ_matrix = joint_BSJ_res %>%
  column_to_rownames(var = "circRNA_ID") %>%
  as.matrix() %>% as.data.frame()
joint_FSJ_matrix = joint_FSJ_res %>%
  column_to_rownames(var = "circRNA_ID") %>%
  as.matrix() %>% as.data.frame()

colnames(joint_BSJ_matrix) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_BSJ_matrix))
colnames(joint_FSJ_matrix) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_FSJ_matrix))
colnames(BSJ_matrix) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(BSJ_matrix))
colnames(FSJ_matrix) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(FSJ_matrix))

cat("Matching BSJ counts:\n")
for (sample in colnames(BSJ_matrix)) {
  matched_circRNAs = intersect(rownames(BSJ_matrix), rownames(joint_BSJ_matrix))
  cat("Number of matched circRNAs:", length(matched_circRNAs), " of ", nrow(BSJ_matrix),"\n")
  original_counts = BSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_BSJ_matrix[matched_circRNAs, sample]
  matching = sum(original_counts == joint_counts, na.rm = TRUE)
  total = length(matched_circRNAs)
  cat(sprintf("Sample %s: %d/%d matching BSJ counts\n", sample, matching, total))
}
cat("Matching FSJ counts:\n")
for (sample in colnames(FSJ_matrix)) {
  matched_circRNAs = intersect(rownames(FSJ_matrix), rownames(joint_FSJ_matrix))
  cat("Number of matched circRNAs:", length(matched_circRNAs), " of ", nrow(FSJ_matrix),"\n")
  original_counts = FSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_FSJ_matrix[matched_circRNAs, sample]
  matching = sum(original_counts == joint_counts, na.rm = TRUE)
  total = length(matched_circRNAs)
  cat(sprintf("Sample %s: %d/%d matching FSJ counts\n", sample, matching, total))
}

# ============================================================================
# EXTENDED SIMILARITY ANALYSIS
# ============================================================================

# --- 1. Correlation Analysis ---
cat("\n=== CORRELATION ANALYSIS ===\n")

for (sample in colnames(BSJ_matrix)) {
  matched_circRNAs = intersect(rownames(BSJ_matrix), rownames(joint_BSJ_matrix))
  original_counts = BSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_BSJ_matrix[matched_circRNAs, sample]
  
  # Pearson correlation
  cor_pearson = cor(original_counts, joint_counts, use = "complete.obs")
  # Spearman correlation (rank-based, better for count data)
  cor_spearman = cor(original_counts, joint_counts, use = "complete.obs", method = "spearman")
  
  cat(sprintf("Sample %s BSJ - Pearson: %.4f, Spearman: %.4f\n", 
              sample, cor_pearson, cor_spearman))
}

for (sample in colnames(FSJ_matrix)) {
  matched_circRNAs = intersect(rownames(FSJ_matrix), rownames(joint_FSJ_matrix))
  original_counts = FSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_FSJ_matrix[matched_circRNAs, sample]
  
  cor_pearson = cor(original_counts, joint_counts, use = "complete.obs")
  cor_spearman = cor(original_counts, joint_counts, use = "complete.obs", method = "spearman")
  
  cat(sprintf("Sample %s FSJ - Pearson: %.4f, Spearman: %.4f\n", 
              sample, cor_pearson, cor_spearman))
}

# --- 2. Agreement Statistics ---
cat("\n=== AGREEMENT STATISTICS ===\n")

compute_agreement_stats = function(original, joint, name) {
  # Mean Absolute Difference
  mad = mean(abs(original - joint), na.rm = TRUE)
  # Root Mean Square Difference
  rmsd = sqrt(mean((original - joint)^2, na.rm = TRUE))
  # Percent within ±1 count
  within_1 = mean(abs(original - joint) <= 1, na.rm = TRUE) * 100
  # Percent within ±10%
  within_10pct = mean(abs(original - joint) <= pmax(original, joint) * 0.1, na.rm = TRUE) * 100
  
  cat(sprintf("%s - MAD: %.2f, RMSD: %.2f, Within ±1: %.1f%%, Within ±10%%: %.1f%%\n",
              name, mad, rmsd, within_1, within_10pct))
  
  return(data.frame(name = name, MAD = mad, RMSD = rmsd, 
                    within_1 = within_1, within_10pct = within_10pct))
}

agreement_results = list()
for (sample in colnames(BSJ_matrix)) {
  matched_circRNAs = intersect(rownames(BSJ_matrix), rownames(joint_BSJ_matrix))
  original_counts = BSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_BSJ_matrix[matched_circRNAs, sample]
  agreement_results[[paste0(sample, "_BSJ")]] = 
    compute_agreement_stats(original_counts, joint_counts, paste(sample, "BSJ"))
}

for (sample in colnames(FSJ_matrix)) {
  matched_circRNAs = intersect(rownames(FSJ_matrix), rownames(joint_FSJ_matrix))
  original_counts = FSJ_matrix[matched_circRNAs, sample]
  joint_counts = joint_FSJ_matrix[matched_circRNAs, sample]
  agreement_results[[paste0(sample, "_FSJ")]] = 
    compute_agreement_stats(original_counts, joint_counts, paste(sample, "FSJ"))
}

agreement_df = bind_rows(agreement_results)

# --- 3. Visualization ---
cat("\n=== GENERATING PLOTS ===\n")

# Create comparison data frame for plotting
comparison_data = list()

for (sample in colnames(BSJ_matrix)) {
  matched_circRNAs = intersect(rownames(BSJ_matrix), rownames(joint_BSJ_matrix))
  comparison_data[[paste0(sample, "_BSJ")]] = data.frame(
    circRNA_ID = matched_circRNAs,
    original = BSJ_matrix[matched_circRNAs, sample],
    joint = joint_BSJ_matrix[matched_circRNAs, sample],
    sample = sample,
    type = "BSJ"
  )
}

for (sample in colnames(FSJ_matrix)) {
  matched_circRNAs = intersect(rownames(FSJ_matrix), rownames(joint_FSJ_matrix))
  comparison_data[[paste0(sample, "_FSJ")]] = data.frame(
    circRNA_ID = matched_circRNAs,
    original = FSJ_matrix[matched_circRNAs, sample],
    joint = joint_FSJ_matrix[matched_circRNAs, sample],
    sample = sample,
    type = "FSJ"
  )
}

comparison_df = bind_rows(comparison_data)

# Scatter plot: Original vs Joint
p1 = ggplot(comparison_df %>% filter(original > 0 & joint > 0), 
            aes(x = original, y = joint)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(type ~ sample, scales = "free") +
  labs(title = "Comparison: Individual vs Joint Analysis",
       x = "Individual Analysis Counts",
       y = "Joint Analysis Counts") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "lightblue"))

p1

# Log-scale scatter plot (better for count data with wide range)
p2 = ggplot(comparison_df %>% filter(original > 0 & joint > 0), 
            aes(x = original, y = joint)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(type ~ sample, scales = "free") +
  labs(title = "Comparison: Individual vs Joint Analysis (Log Scale)",
       x = "Individual Analysis Counts (log10)",
       y = "Joint Analysis Counts (log10)") +
  theme_bw()

p2

# Bland-Altman plot (difference vs mean)
comparison_df = comparison_df %>% filter(original > 0 & joint > 0) %>%
  mutate(mean_count = (original + joint) / 2,
         difference = original - joint)

p3 = ggplot(comparison_df, aes(x = mean_count, y = difference)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(type ~ sample, scales = "free") +
  labs(title = "Bland-Altman Plot: Agreement Between Methods",
       x = "Mean Count",
       y = "Difference (Individual - Joint)") +
  theme_bw()

p3

# --- 4. Venn Diagram / Overlap Analysis ---
cat("\n=== OVERLAP ANALYSIS ===\n")

# CircRNAs detected in each method
individual_circRNAs = unique(res$circRNA_ID)
joint_circRNAs = rownames(joint_BSJ_matrix)

only_individual = setdiff(individual_circRNAs, joint_circRNAs)
only_joint = setdiff(joint_circRNAs, individual_circRNAs)
both = intersect(individual_circRNAs, joint_circRNAs)

cat(sprintf("CircRNAs in individual analysis only: %d\n", length(only_individual)))
cat(sprintf("CircRNAs in joint analysis only: %d\n", length(only_joint)))
cat(sprintf("CircRNAs in both analyses: %d\n", length(both)))
cat(sprintf("Jaccard similarity: %.4f\n", 
            length(both) / length(union(individual_circRNAs, joint_circRNAs))))

# --- 5. Summary Statistics Table ---
cat("\n=== SUMMARY TABLE ===\n")

summary_table = comparison_df %>%
  group_by(sample, type) %>%
  summarise(
    n_circRNAs = n(),
    exact_match = sum(original == joint, na.rm = TRUE),
    exact_match_pct = mean(original == joint, na.rm = TRUE) * 100,
    pearson_cor = cor(original, joint, use = "complete.obs"),
    spearman_cor = cor(original, joint, use = "complete.obs", method = "spearman"),
    mean_diff = mean(original - joint, na.rm = TRUE),
    sd_diff = sd(original - joint, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# --- 6. Heatmap of discrepancies ---
cat("\n=== IDENTIFYING LARGEST DISCREPANCIES ===\n")

top_discrepancies = comparison_df %>%
  mutate(abs_diff = abs(difference),
         pct_diff = ifelse(mean_count > 0, abs_diff / mean_count * 100, NA)) %>%
  group_by(sample, type) %>%
  slice_max(order_by = abs_diff, n = 10) %>%
  select(circRNA_ID, sample, type, original, joint, difference, pct_diff)

print(top_discrepancies)


################################################################################
################################################################################
################################################################################

library(tidyverse)

setwd("/pastel/tools/circRNA_tools/test_data")
res1 = read_tsv("Div_100_S91.CIRI3results.txt", show_col_types = FALSE)
res2 = read_tsv("Div_101_S92.CIRI3results.txt", show_col_types = FALSE)
res = bind_rows(res1 %>% mutate(Sample = "Div_100_S91"),
                res2 %>% mutate(Sample = "Div_101_S92")) 

pardos_res = list.files("/pastel/projects/PARDoS/circrna_pardos", pattern = ".CIRI3", full.names = T)
names(pardos_res) = basename(pardos_res)
pardos_res_df = purrr::map_df(pardos_res, ~ read_tsv(.x, show_col_types = FALSE), .id = "file") %>%
  mutate(Sample = gsub("(.*)\\.CIRI3", "\\1", basename(file))) %>%
  select(-file)

res = bind_rows(res, pardos_res_df) %>%
  filter(chr %in% paste0("chr", c(1:22, "X", "Y")))
table(res$chr)

# Check the overlap of circRNAs detected in each sample (upset)
cat("CircRNAs detection overlap\n")
library(UpSetR)
circ_list <- res %>%
  # filter(Sample %in% c("PARDOS_2","PARDOS_1")) %>%
  select(circRNA_ID, Sample) %>%
  distinct() %>%
  group_by(Sample) %>%
  summarise(circRNAs = list(circRNA_ID), .groups = "drop") %>%
  { setNames(.$circRNAs, .$Sample) }
UpSetR::upset(
  fromList(circ_list),
  order.by = "freq",
  nsets    = length(circ_list)
)

circ_counts <- res %>%
  select(circRNA_ID, Sample) %>%
  distinct() %>%
  group_by(circRNA_ID) %>%
  summarise(n_samples = n_distinct(Sample), .groups = "drop")

pct_shared <- mean(circ_counts$n_samples >= 10) * 100
cat(sprintf("%.1f%% of circRNAs detected in >= 10 samples (%d / %d)\n",
            pct_shared,
            sum(circ_counts$n_samples >= 10),
            nrow(circ_counts)))

##### Testing a merging strategy 

res1 = read_tsv("Div_100_S91.CIRI3results.txt", show_col_types = FALSE)
res2 = read_tsv("Div_101_S92.CIRI3results.txt", show_col_types = FALSE)
single_results = bind_rows(res1 %>% mutate(Sample = "Div_100_S91"),
                res2 %>% mutate(Sample = "Div_101_S92")) 
length(unique(single_results$circRNA_ID))

joint_results = read_tsv("CIRI3_joint_results.txt")
length(unique(joint_results$circRNA_ID))

# Compare results from joint analysis 
joint_BSJ_res = read_tsv("CIRI3_joint_results.txt.BSJ_Matrix")
joint_FSJ_res = read_tsv("CIRI3_joint_results.txt.FSJ_Matrix")

joint_FSJ_res_MERGE = data.table::fread("final_output.FSJ_Matrix", fill = TRUE)
empty_cols <- names(joint_FSJ_res_MERGE)[colSums(!is.na(joint_FSJ_res_MERGE)) == 0]
joint_FSJ_res_MERGE[, (empty_cols) := NULL]

joint_BSJ_res_MERGE = data.table::fread("final_output.BSJ_Matrix", fill = TRUE)
empty_cols_bsj <- names(joint_BSJ_res_MERGE)[colSums(!is.na(joint_BSJ_res_MERGE)) == 0]
joint_BSJ_res_MERGE[, (empty_cols_bsj) := NULL]

colnames(joint_FSJ_res) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_FSJ_res))
colnames(joint_FSJ_res_MERGE) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_FSJ_res_MERGE))
colnames(joint_BSJ_res) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_BSJ_res))
colnames(joint_BSJ_res_MERGE) = gsub("(.*)Div_(.*?)_(.*)", "Div_\\2", colnames(joint_BSJ_res_MERGE))

colnames(joint_FSJ_res_MERGE) = colnames(joint_FSJ_res)
colnames(joint_BSJ_res_MERGE) = colnames(joint_BSJ_res)

FSJ_combined = joint_FSJ_res %>%
  full_join(joint_FSJ_res_MERGE, by = "circRNA_ID", suffix = c("_joint", "_merge"))
head(FSJ_combined)

BSJ_combined = joint_BSJ_res %>%
  full_join(joint_BSJ_res_MERGE, by = "circRNA_ID", suffix = c("_joint", "_merge"))
head(BSJ_combined)

# plot scatter plot of FSJ counts from joint vs merge
FSJ_combined %>%
  pivot_longer(cols = -circRNA_ID, names_to = "sample", values_to = "FSJ_count") %>%
  mutate(method = gsub("(.*)_(joint|merge)", "\\2", sample),
         sample = gsub("(.*)_(joint|merge)", "\\1", sample)) %>%
  pivot_wider(names_from = method, values_from = FSJ_count) %>%
  ggplot(aes(x = joint, y = merge)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "FSJ Counts: Joint vs Merge", x = "Joint Analysis FSJ Count", y = "Merge Analysis FSJ Count") +
  theme_bw() +
  facet_wrap(~ sample, scales = "free")

# plot scatter plot of BSJ counts from joint vs merge
BSJ_combined %>%
  pivot_longer(cols = -circRNA_ID, names_to = "sample", values_to = "BSJ_count") %>%
  mutate(method = gsub("(.*)_(joint|merge)", "\\2", sample),
         sample = gsub("(.*)_(joint|merge)", "\\1", sample)) %>%
  pivot_wider(names_from = method, values_from = BSJ_count) %>%
  ggplot(aes(x = joint, y = merge)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "BSJ Counts: Joint vs Merge", x = "Joint Analysis BSJ Count", y = "Merge Analysis BSJ Count") +
  theme_bw() +
  facet_wrap(~ sample, scales = "free")



cd /pastel/tools/circRNA_tools/test_data
THREADS=8
REF_FASTA=GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=gencode.v32.primary_assembly.annotation.gtf
GENOMEDIR=${PWD}
CIRI3_JAR_PATH=/pastel/tools/circRNA_tools/CIRI3/CIRI3_Java_18.0.1.jar

# List of sample prefixes
SAMPLES=(
  "Div_100_S91"
  "Div_101_S92"
)

# Processing CIRI3 by sample
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${PWD}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${PWD}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
  -O ${SAMPLE}.CIRI3results_TEST.txt \
  -F ${REF_FASTA} \
  -A ${GTF_FILE} \
  -Ma 1 \
  -W 0 \
  --keep-bsj \
  -T ${THREADS}
done

echo -e "/path/sample1.bam\n/path/sample2.bam" > sam_list.txt
$CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} \
JoinSamples \
-I sam_list.txt \
-O joined \
-F ${REF_FASTA} \
-A ${GTF_FILE}
# For STAR: add -Ma 1



java -jar CIRI3.jar -W 0 --keep-bsj -I sample1.bam -O sample1.result -F ref.fa -A anno.gtf
java -jar CIRI3.jar -W 0 --keep-bsj -I sample2.bam -O sample2.result -F ref.fa -A anno.gtf
# For STAR: add -Ma 1
