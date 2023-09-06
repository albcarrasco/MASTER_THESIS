#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(tidyr)

# Read the TSV file
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

data <- read.table(input_file, header = TRUE, sep = "\t")

# Define the order of parameters
parameter_order <- c("mapped", "on.target", "detected_targets", "cov_uniformity")

# Reshape the data for plotting
data_long <- pivot_longer(data, cols = c(`mapped`, `on.target`, `detected_targets`, `cov_uniformity`),
                          names_to = "parameter", values_to = "percentage")

# Convert the parameter factor to the defined order
data_long$parameter <- factor(data_long$parameter, levels = parameter_order)

# Define custom colors
custom_colors <- c("#336699", "#86BBD8", "#2F4858", "#9EE493")

# Define custom labels for the legend
legend_labels <- c("Mapped reads", "On-Target reads", "Detected targets", "Coverage uniformity")

# Create a grouped bar plot with modified settings
p <- ggplot(data_long, aes(x = parameter, y = percentage, fill = parameter)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +  # Adjust position and width
  scale_fill_manual(values = custom_colors, labels = legend_labels) +  # Apply custom colors and labels
  labs(x = NULL,
       y = "Percentage",
       fill = NULL) +  # Remove legend title
  facet_wrap(~ sample, ncol = 1, scales = "free_x") +  # Create separate plots for each sample
  theme_classic() +
  theme(legend.position = "top",
        # legend.box.background = element_rect(color = "black"),  # Add legend box border color
        legend.box.margin = margin(6, 6, 6, 6),  # Add margin to the legend box
        plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and make title bold
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.title = element_blank(),  # Remove legend title
        strip.text = element_text(size = 13),  # Adjust facet title font size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_line(color = "lightgray"),  # Add minor grid lines
        panel.background = element_blank()) +  # Remove plot background
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") +  # Add dashed line at 90%
  scale_x_discrete(labels = NULL) +  # Remove labels on the x-axis
  labs(title = "Flomics nf-core/sarek mapping QC metrics",  # Set static title
       x = "Metrics",  # Set x-axis label
       y = "Percentage",
       fill = NULL)

# Save the plot as a png file
output_file <- gsub(".tsv$", ".png", input_file)
ggsave(output_file, plot = p, width = 10, height = 6)


