# Example data for heatmap
heatmap_data <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, byrow = TRUE)
heatmap_df <- as.data.frame(heatmap_data)

# Ensure that the dataframe contains only numeric values
heatmap_df <- sapply(heatmap_df, as.numeric)

# Example data for custom distance values
distance_data <- matrix(c(0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0), nrow = 4, byrow = TRUE)
distance_df <- as.data.frame(distance_data)

# Perform hierarchical clustering based on custom distance values
hclust_result <- hclust(dist(distance_df), method = "average")

# Convert the dataframe to a numeric matrix
heatmap_matrix <- as.matrix(heatmap_df)

# Create a heatmap with dendrogram
heatmap(heatmap_matrix, Rowv = hclust_result, Colv = FALSE, col = heat.colors(256), scale = "row", main = "Heatmap with Dendrogram")
