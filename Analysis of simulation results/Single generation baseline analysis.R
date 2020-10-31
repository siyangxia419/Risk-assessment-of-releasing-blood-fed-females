### Visualize the dengue cases transmitted by Aedes aegypti in the first generations under four different releasing scenario
### 2020.9.26
### Input data: results from the python simulation

library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)



# Data management ---------------------------------------------------------

# Paths to the simulation output:
bloodfed.dir <- "../Bloodfed_specific_parameter_values/"

# Baseline results
single.gen.results <- list(
  results1 = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_MN_Gens_2001.csv"),
                      header = FALSE, stringsAsFactors = FALSE),
  results2 = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_MN_Gens_2002.csv"),
                      header = FALSE, stringsAsFactors = FALSE),
  results3 = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_MN_Gens_2003.csv"),
                      header = FALSE, stringsAsFactors = FALSE),
  results4 = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_MN_Gens_2004.csv"),
                      header = FALSE, stringsAsFactors = FALSE)
)
n <- nrow(single.gen.results$results1)

# the 1st column contains the number of infections from the wild mosquitoes in the first generation
# the 11th column contains the number of infections from the released mosquitoes in the first generation

d <- lapply(X = single.gen.results, 
            FUN = function(x) {
              y <- x %>% 
                select(wild = V1, released = V11) %>% 
                mutate(iteration = 1:nrow(x)) %>%
                melt(id.var = "iteration")
            }) %>% 
  bind_rows(.id = "id") %>%
  mutate(id = ordered(id, levels = paste0("results", 1:4),
                      labels = c("VC of release = 0.2\nsize of release = 0.1",
                                 "VC of release = 0.05\nsize of release = 0.1",
                                 "VC of release = 0.2\nsize of release = 0.2",
                                 "VC of release = 0.05\nsize of release = 0.2"))) %>%
  rename(Attribution = variable)


# The mean of each group
d.summary <- d %>% group_by(id, Attribution) %>%
  summarise(mean = mean(value)) %>%
  mutate(y = 25, label = paste("mean =", round(mean, digits = 0)))
d.summary$label[d.summary$Attribution == "released"] <- 
  paste0(d.summary$label[d.summary$Attribution == "released"], 
         " (", 
         round(d.summary$mean[d.summary$Attribution == "released"] / 
                 d.summary$mean[d.summary$Attribution == "wild"], digits = 2),
         ")")
  


# Visualization -----------------------------------------------------------

# text size:
et <- element_text(size = 14)

# A colorblind-friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sel.col <- cbPalette[c(2, 3)]

# Second axis to show the proportion relative to the wild mosquitoes
sec.axis.scale <- d.summary$mean[1]

# The plot:
single.gen.plot <- ggplot(data = d, 
                          aes(y = value, fill = Attribution)) + 
  geom_hline(data = d.summary,
             aes(yintercept = mean, color = Attribution),
             linetype = "dashed", size = 0.5) +
  geom_text(data = d.summary, aes(y = mean, x = y, label = label),
            size = 4, nudge_y = 100) + 
  geom_histogram(binwidth = 10) + 
  facet_wrap(vars(id), nrow = 1, ncol = 4, strip.position = "top") + 
  scale_y_continuous(name = "Number of human infections in the first generation",
                     sec.axis = sec_axis(trans=~./sec.axis.scale, 
                                         name="Porpotion relative to the wild mosquitoes",
                                         breaks = seq(0, 1, 0.2))) + 
  scale_x_continuous(name = "Frequency among all iterations") + 
  scale_fill_manual(values = sel.col, 
                    labels = c("Wild mosquitoes", "Released mosquitoes")) + 
  scale_color_manual(values = sel.col, 
                     labels = c("Wild mosquitoes", "Released mosquitoes"),
                     guide = "none") + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.minor.y = element_line(color = "gray95"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.y = et,
    axis.text.y = et,
    axis.title.x = et,
    axis.text.x = et,
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )
print(single.gen.plot)
ggsave(filename = "Single_generation_baseline_histogram.pdf", 
       plot = single.gen.plot,
       width = 11, height = 8, units = "in", dpi = 300)
ggsave(filename = "Single_generation_baseline_histogram.png", 
       plot = single.gen.plot,
       width = 11, height = 8, units = "in", dpi = 300)





# Numerical summary -------------------------------------------------------

d.summary2 <- d %>% group_by(id, Attribution) %>%
  summarise(mean = mean(value), 
            quantile_low = quantile(value, probs = 0.025), 
            quantile_high = quantile(value, probs = 0.975))
baseline_wild <- rep(d.summary2$mean[d.summary2$Attribution =="wild"], each = 2)
d.summary2 <- d.summary2 %>% ungroup() %>%
  mutate(mean_prop = mean / baseline_wild,
         quantile_low_prop = quantile_low / baseline_wild,
         quantile_high_prop = quantile_high / baseline_wild)
write.csv(x = d.summary2, 
          file = "Single generation baseline analysis_numerical summary.csv", 
          row.names = FALSE, quote = FALSE)


save.image("Single generation baseline analysis_20201005.RData")
