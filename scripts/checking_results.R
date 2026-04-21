library(tidyverse)

setwd("/pastel/tools/circRNA_tools/test_data/decoupled_comparison")
# setwd("/pastel/tools/circRNA_tools/test_data/decoupled_bwa_comparison")

# Compare results from joint analysis 
joint_res_df = read_tsv("original/result")
head(joint_res_df)
joint_BSJ_res = read_tsv("original/result.BSJ_Matrix")
joint_FSJ_res = read_tsv("original/result.FSJ_Matrix")

# Compare results from decoupled analysis
decoupled_res_df = read_tsv("decoupled/finalize/result")
head(decoupled_res_df)
decoupled_BSJ_res = read_tsv("decoupled/finalize/result.BSJ_Matrix")
decoupled_FSJ_res = read_tsv("decoupled/finalize/result.FSJ_Matrix")

colnames(joint_BSJ_res) = colnames(decoupled_BSJ_res)
colnames(joint_FSJ_res) = colnames(decoupled_FSJ_res)

res_BSJ_merged = full_join(joint_BSJ_res, decoupled_BSJ_res, by = "circRNA_ID", suffix = c(".joint", ".decoupled"))
res_BSJ_merged_long = res_BSJ_merged %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample", "Run"), names_pattern = "(.*)\\.(.*)") %>%
  pivot_wider(names_from = Run, values_from = value)

ggplot(res_BSJ_merged_long, aes(x = joint, y = decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Sample) +
  labs(x = "Joint Analysis BSJ Count", y = "Decoupled Analysis BSJ Count") +
  theme_bw()

res_FSJ_merged = full_join(joint_FSJ_res, decoupled_FSJ_res, by = "circRNA_ID", suffix = c(".joint", ".decoupled"))
res_FSJ_merged_long = res_FSJ_merged %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample", "Run"), names_pattern = "(.*)\\.(.*)") %>%
  pivot_wider(names_from = Run, values_from = value)

ggplot(res_FSJ_merged_long, aes(x = joint, y = decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Sample) +
  labs(x = "Joint Analysis FSJ Count", y = "Decoupled Analysis FSJ Count") +
  theme_bw()

# Check the overlap of circRNAs detected in each sample (upset)
cat("CircRNAs detection overlap\n")
library(UpSetR)
circ_list <- res_BSJ_merged_long %>%
  filter(decoupled > 0) %>%
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

# Compare PARDOS results with Jishu results
pardos_res = list.files("/pastel/projects/PARDoS/circrna_pardos", pattern = ".CIRI3", full.names = T)
names(pardos_res) = basename(pardos_res)
pardos_res_df = purrr::map_df(pardos_res, ~ read_tsv(.x, show_col_types = FALSE), .id = "file") %>%
  mutate(Sample = gsub("(.*)\\.CIRI3", "\\1", basename(file))) %>%
  select(-file)

pardos_check = pardos_res_df %>% filter(Sample %in% c("PARDOS_1","PARDOS_2")) %>%
  select(Sample, circRNA_ID, `#junction_reads`) %>%
  pivot_wider(names_from = Sample, values_from = `#junction_reads`) %>%
  full_join(
    res_BSJ_merged_long %>% 
      mutate(Sample = gsub("(.*)_(.*)","\\1",Sample)) %>%
      filter(Sample %in% c("PARDOS_1","PARDOS_2")) %>% 
      pivot_wider(names_from = Sample, values_from = decoupled) %>%
      select(circRNA_ID, PARDOS_1, PARDOS_2), 
    by = c("circRNA_ID"), suffix = c(".jishu",".decoupled")) 

ggplot(pardos_check, aes(x = PARDOS_1.jishu, y = PARDOS_1.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "PARDOS_1 CIRI3 BSJ Count", y = "Decoupled Analysis BSJ Count") +
  theme_bw()
ggplot(pardos_check, aes(x = PARDOS_2.jishu, y = PARDOS_2.decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "PARDOS_2 CIRI3 BSJ Count", y = "Decoupled Analysis BSJ Count") +
  theme_bw()
    