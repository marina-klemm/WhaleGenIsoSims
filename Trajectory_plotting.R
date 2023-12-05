#November 16th, 2023
#Plotting trajectory from the whale_anewtry_restored.py script
#I want to plot the allele frequencies from the aforementioned script, to
#first visualize them before the script is complete.

#I exported the data by pasting them to excel, for ease of handling at first.


# Load libraries ----------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(reshape)
library(patchwork)

# Load dataset ------------------------------------------------------------
subpop0_model3 <- read.table("C:/Users/marin/anaconda3/envs/simuPOP-env/subpop0_model3.txt",
                    sep = "\t", header = TRUE)

subpop1_model3 <- read.table("C:/Users/marin/anaconda3/envs/simuPOP-env/subpop1_model3.txt",
                      sep = "\t", header = TRUE)

# Reshape the dataset -----------------------------------------------------

long_subpop0_model3 <- data.frame(x = seq_along(subpop0_model3[,1]),
                           subpop0_model3)
long_subpop0_model3 <- melt(long_subpop0_model3, id.vars = "x")

long_subpop1_model3 <- data.frame(x = seq_along(subpop1_model3[,1]),
                           subpop1_model3)
long_subpop1_model3 <- melt(long_subpop1_model3, id.vars = "x")


# Plot --------------------------------------------------------------------

sub0plot_model3 <- ggplot(long_subpop0_model3, aes(x = x, y = value, colour = variable)) +
  geom_line() +
  labs(
    title = "Simulated trajectory for subpopulation 0 (Model 3)",
    x = "Generation",
    y = "Allele frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title with increased font size
    axis.title = element_text(size = 14)  # Increased font size for axis labels
  ) +
  scale_x_continuous(breaks = 1:10) 


sub1plot_model3 <- ggplot(long_subpop1_model3, aes(x = x, y = value, colour = variable)) +
  geom_line() +
  labs(
    title = "Simulated trajectory for subpopulation 1 (Model 3)",
    x = "Generation",
    y = "Allele frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title with increased font size
    axis.title = element_text(size = 14)  # Increased font size for axis labels
  ) +
  scale_x_continuous(breaks = 1:10) 
  
combined_plot_model3 <- sub0plot_model3 / sub1plot_model3
#combined_plot_vertical <- combined_plot / NULL
print(combined_plot_model3) #horizontal view





# Now, using model4 instead: ----------------------------------------------


# Load dataset ------------------------------------------------------------
subpop0_model4 <- read.table("C:/Users/marin/anaconda3/envs/simuPOP-env/subpop0_model4.txt",
                      sep = "\t", header = TRUE)

subpop1_model4 <- read.table("C:/Users/marin/anaconda3/envs/simuPOP-env/subpop1_model4.txt",
                      sep = "\t", header = TRUE)

# Reshape the dataset -----------------------------------------------------

long_subpop0_model4 <- data.frame(x = seq_along(subpop0_model4[,1]),
                           subpop0_model4)
long_subpop0_model4 <- melt(long_subpop0_model4, id.vars = "x")

long_subpop1_model4 <- data.frame(x = seq_along(subpop1_model4[,1]),
                           subpop1_model4)
long_subpop1_model4 <- melt(long_subpop1_model4, id.vars = "x")


# Plot --------------------------------------------------------------------

sub0plot_model4 <- ggplot(long_subpop0_model4, aes(x = x, y = value, colour = variable)) +
  geom_line() +
  labs(
    title = "Simulated trajectory for subpopulation 0 (Model 4)",
    x = "Generation",
    y = "Allele frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title with increased font size
    axis.title = element_text(size = 14)  # Increased font size for axis labels
  ) +
  scale_x_continuous(breaks = 1:10) 


sub1plot_model4 <- ggplot(long_subpop1_model4, aes(x = x, y = value, colour = variable)) +
  geom_line() +
  labs(
    title = "Simulated trajectory for subpopulation 1 (Model 4)",
    x = "Generation",
    y = "Allele frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title with increased font size
    axis.title = element_text(size = 14)  # Increased font size for axis labels
  ) +
  scale_x_continuous(breaks = 1:10) 

combined_plot_model4 <- sub0plot_model4 / sub1plot_model4
#combined_plot_vertical <- combined_plot / NULL
print(combined_plot_model4) #horizontal view























