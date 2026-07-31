# ============================================================
# Temporal Trajectory Clustering
# ============================================================
#
# Cluster subjects based on temporal changes between consecutive
# timepoints using hierarchical clustering.
#
# Author: Lu Sun
# License: UCLA
#
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cluster_temporal_trajectories <- function(
    data,
    patient_col = "Patient",
    timepoint_col = "Timepoint",
    value_col = "Value",
    timepoint_order,
    cluster_k = 2,
    hclust_method = "ward.D2",
    dist_method = "euclidean",
    na_imputation = c("mean", "none")
) {
  
  na_imputation <- match.arg(na_imputation)
  
  # ----------------------------------------------------------
  # Standardize column names
  # ----------------------------------------------------------
  
  df <- data %>%
    dplyr::select(
      Patient = all_of(patient_col),
      Timepoint = all_of(timepoint_col),
      Value = all_of(value_col)
    )
  
  # ----------------------------------------------------------
  # Wide format
  # ----------------------------------------------------------
  
  df_wide <- df %>%
    pivot_wider(
      names_from = Timepoint,
      values_from = Value
    ) %>%
    group_by(Patient) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  
  # Keep specified timepoints
  available_tp <- intersect(timepoint_order, colnames(df_wide))
  
  df_matrix <- df_wide %>%
    dplyr::select(Patient, all_of(available_tp))
  
  patient_ids <- df_matrix$Patient
  
  trajectory_mat <- as.data.frame(df_matrix[, -1])
  
  rownames(trajectory_mat) <- patient_ids
  
  # ----------------------------------------------------------
  # Compute temporal deltas
  # ----------------------------------------------------------
  
  delta_mat <- data.frame(Patient = patient_ids)
  
  for (i in seq_len(length(available_tp) - 1)) {
    
    tp1 <- available_tp[i]
    tp2 <- available_tp[i + 1]
    
    delta_name <- paste0("d", i)
    
    delta_mat[[delta_name]] <-
      trajectory_mat[[tp2]] - trajectory_mat[[tp1]]
  }
  
  # ----------------------------------------------------------
  # Missing value handling
  # ----------------------------------------------------------
  
  if (na_imputation == "mean") {
    
    delta_cols <- setdiff(colnames(delta_mat), "Patient")
    
    for (col in delta_cols) {
      delta_mat[[col]][is.na(delta_mat[[col]])] <-
        mean(delta_mat[[col]], na.rm = TRUE)
    }
  }
  
  # ----------------------------------------------------------
  # Clustering
  # ----------------------------------------------------------
  
  dist_matrix <- dist(
    delta_mat[, -1],
    method = dist_method
  )
  
  hc <- hclust(
    dist_matrix,
    method = hclust_method
  )
  
  clusters <- cutree(hc, k = cluster_k)
  
  # ----------------------------------------------------------
  # Attach cluster membership
  # ----------------------------------------------------------
  
  df_clustered <- trajectory_mat
  
  df_clustered$Patient <- patient_ids
  df_clustered$Cluster <- factor(clusters)
  
  # ----------------------------------------------------------
  # Long format for plotting
  # ----------------------------------------------------------
  
  df_long <- df_clustered %>%
    pivot_longer(
      cols = all_of(available_tp),
      names_to = "Timepoint",
      values_to = "Value"
    )
  
  df_long$Timepoint <-
    factor(df_long$Timepoint,
           levels = timepoint_order)
  
  # ----------------------------------------------------------
  # Summary statistics
  # ----------------------------------------------------------
  
  summary_df <- df_long %>%
    group_by(Cluster, Timepoint) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
  # Trajectory plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    summary_df,
    aes(
      x = Timepoint,
      y = mean,
      color = Cluster,
      fill = Cluster,
      group = Cluster
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(
      aes(
        ymin = mean - sd,
        ymax = mean + sd
      ),
      alpha = 0.2,
      color = NA
    ) +
    theme_minimal() +
    labs(
      x = "Timepoint",
      y = "Value",
      title = "Trajectory Clusters"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  return(list(
    hclust = hc,
    clusters = clusters,
    patient_clusters = df_clustered,
    delta_matrix = delta_mat,
    summary = summary_df,
    plot = p
  ))
}