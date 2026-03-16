library(tidyverse)

# 1. load data
res <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/TARGET_ALL_Bisque.rds')
ratio_matrix <- res$bulk.props

# 2. subset cells
target_cells <- c('High', 'Neutral', 'Low')
target_ratio <- ratio_matrix[target_cells, ]

# 3. check sum ratio
col_sums <- colSums(target_ratio)
print(summary(col_sums))

# 4. find samples which sum ratio equals 0
zero_sum_samples <- names(col_sums[col_sums == 0])
if (length(zero_sum_samples) > 0) {
  print(zero_sum_samples)
  print(ratio_matrix[, zero_sum_samples, drop = FALSE])
}

# 5. del sample  which sum ratio equals 0
target_ratio_noZero <- target_ratio[, !colnames(target_ratio) %in% zero_sum_samples, drop = FALSE]
target_ratio_norm <- apply(target_ratio_noZero, 2, function(x) x / sum(x))

# 6. K-means
set.seed(123)
kmeans_result <- kmeans(target_ratio_t, centers = 3, nstart = 25)
saveRDS(kmeans_result,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL.RDS')

# 7. change class number into subtype
kmeans_result <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL.RDS')
class <- kmeans_result$cluster
subtype <- data.frame(class = class, row.names = names(class))
subtype <- subtype %>%
  mutate(subtype = case_when(
    class == 1 ~ "Low",
    class == 2 ~ "Middle",
    class == 3 ~ "High"
  ))
saveRDS(subtype,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')


