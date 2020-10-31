### Visualize the dengue cases transmitted by Aedes aegypti across 10 generations under four different releasing scenario
### 2020.9.20
### Input data: results from the python simulation

library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)



# Data management ---------------------------------------------------------

# Paths to the simulation output:
bloodfed.dir <- "../Bloodfed_specific_parameter_values/"
bothsex.dir <- "../Bothsex_specific_parameter_values/"
maleonly.dir <- "../Maleonly_specific_parameter_values/"
norelease.dir <- "../Norelease_specific_parameter_values/"

# Choose which parameter combination to visualize:
index <- 4
name_index <- as.character(index + 2000)


# Baseline results
baseline.results <- list(
  bloodfed = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_MN_Gens_", name_index, ".csv"),
                      header = FALSE, stringsAsFactors = FALSE),
  bothsex = read.csv(paste0(bothsex.dir, "simulation_Bothsex_MN_Gens_", name_index, ".csv"),
                     header = FALSE, stringsAsFactors = FALSE),
  maleonly = read.csv(paste0(maleonly.dir, "simulation_Maleonly_MN_Gens_", name_index, ".csv"),
                      header = FALSE, stringsAsFactors = FALSE),
  norelease = read.csv(paste0(norelease.dir, "simulation_NoRelease_MN_Gens_", name_index, ".csv"),
                      header = FALSE, stringsAsFactors = FALSE)
)

n <- nrow(baseline.results$bloodfed)


# Some classification labels:
strategies <- names(baseline.results)
outcomes <- c("case.wild", "case.release", "mosq.wild", "mosq.release")
metrics <- c("cases", "cases", "mosquitoes", "mosquitoes")
populations <- c("wild", "release", "wild", "release")


# Separating the different output metrics:
d <- vector(mode = "list", length = 4 * 4)
names(d) <- paste(rep(strategies, each = 4), rep(outcomes, times = 4), sep = ".")

# The first 10 columns: the number of new dengue cases contributed by the wild mosquitoes
# The second 10 columns: the number of new dengue cases contributed by the released mosquitoes
# The third 10 columns: the number of infectious mosquitoes in the wild population
# The last 10 columns: the number of infectious mosquitoes in the releases population
# The 10 columns in each group represents the corresponding values for the 10 generations of release

for(i in strategies){
  for(j in 1:4){
    temp.name <- paste(i, outcomes[j], sep = ".")
    d[[temp.name]] <- baseline.results[[i]][, (1:10 + 10*(j-1))] %>%
      mutate(iteration = 1:n, 
             strategy = i, 
             metric = metrics[j],
             population = populations[j]) 
    names(d[[temp.name]])[1:10] <- paste0("Gen.", 1:10)
  }
}


# Combine all dengue case output:
d.case <- d[grepl(pattern = "case", x = names(d))] %>% 
  bind_rows() %>%
  melt(id.var = c("iteration", "strategy", "metric", "population")) %>%
  mutate(generation = as.numeric(gsub(pattern = "Gen.", replacement = "", x = variable))) %>%
  select(iteration, strategy, population, generation, value) %>%
  dcast(formula = iteration + strategy + generation ~ population, value.var = "value") %>%
  arrange(iteration, generation) %>%
  mutate(total = wild + release, increase = ifelse(wild > 0, release / wild, 0))


# summarize mean, SD and 95% interval:
d.summary <- d.case %>% group_by(strategy, generation) %>%
  summarise(wild_mean = mean(wild), 
            wild_median = median(wild), 
            wild_sd = sd(wild), 
            wild_IQR = IQR(wild),
            wild_quantile_low = quantile(wild, 0.025), 
            wild_quantile_high = quantile(wild, 0.975),
            release_mean = mean(release), 
            release_median = median(release), 
            release_sd = sd(release), 
            release_IQR = IQR(release),
            release_quantile_low = quantile(release, 0.025), 
            release_quantile_high = quantile(release, 0.975),
            total_mean = mean(total), 
            total_median = median(total), 
            total_sd = sd(total), 
            total_IQR = IQR(total),
            total_quantile_low = quantile(total, 0.025), 
            total_quantile_high = quantile(total, 0.975),
            increase_mean = mean(increase), 
            increase_median = median(increase), 
            increase_sd = sd(increase), 
            increase_IQR = IQR(increase),
            increase_quantile_low = quantile(increase, 0.025), 
            increase_quantile_high = quantile(increase, 0.975)) 
d.summary <- d.summary %>%
  mutate(strategy = ordered(strategy, levels = c("bloodfed", "bothsex", "maleonly", "norelease")))




# Visualization -----------------------------------------------------------

# shared options:
col.scenarios <- c("gray70", brewer.pal(n = 3, name = "Set2")[c(2, 3, 1)])
pd <- position_dodge(0.5)
et <- element_text(size = 14)



# Figure 1: the change of new dengue infections through generations -------

d.plot1 <- d.summary %>%
  select(generation, strategy, total_mean, total_quantile_low, total_quantile_high) %>%
  rename(mean = total_mean, quantile_low = total_quantile_low, quantile_high = total_quantile_high) %>%
  mutate(attribution = "total")

d.plot2 <- d.summary %>%
  select(generation, strategy, wild_mean, wild_quantile_low, wild_quantile_high) %>%
  rename(mean = wild_mean, quantile_low = wild_quantile_low, quantile_high = wild_quantile_high) %>%
  mutate(quantile_low = NA, quantile_high = NA) %>%
  filter(strategy %in% c("bloodfed", "bothsex")) %>%
  mutate(attribution = "wild\nmosquitoes")
# d.plot2$mean[d.plot2$strategy %in% c("maleonly", "norelease")] <- NA

d.plot3 <- d.summary %>%
  select(generation, strategy, release_mean, release_quantile_low, release_quantile_high) %>%
  rename(mean = release_mean, quantile_low = release_quantile_low, quantile_high = release_quantile_high) %>%
  # mutate(quantile_low = NA, quantile_high = NA) %>%
  filter(strategy %in% c("bloodfed", "bothsex")) %>%
  mutate(attribution = "released\nmosquitoes")

d.plot <- bind_rows(d.plot1, d.plot2, d.plot3) %>% 
  mutate(strategy = ordered(strategy, levels = c("norelease", "bothsex", "maleonly", "bloodfed")),
         attribution = ordered(attribution, levels = c("total", "wild\nmosquitoes", "released\nmosquitoes")))

f1.1 <- ggplot(data = d.plot, aes(x = generation, y = mean, 
                                  color = strategy, 
                                  linetype = attribution,
                                  shape = attribution, 
                                  group = interaction(attribution, strategy))) + 
  # Total cases
  geom_errorbar(aes(ymin = quantile_low, ymax = quantile_high), 
                width = 0.3, position = pd) + 
  geom_line(size = 0.6, position = pd) +
  geom_point(size = 2, position = pd) +
  
  scale_x_continuous(name="Generation of releases", breaks = 1:10) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("no release", "un-fed", "only male", "blood-fed")) + 
  scale_shape_manual(name = "Attribution", values = c(19, 1, 8)) + 
  scale_linetype_manual(name = "Attribution", values = c(1, 2, 3)) + 
  ylab("Number of new dengue infections") + 
  
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
    legend.text = element_text(size = 10)
    # legend.key.height = unit(1.8, "line")
  )
print(f1.1)
ggsave(filename = paste0("Baseline_", index, "_generation_mean and 95 percent CI_v1.png"),
       plot = f1.1, 
       width = 9, height = 5, units = "in", dpi = 320)



# Alternative Figure 1:
d.plot <- d.summary
d.plot$release_IQR[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release_mean[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release_sd[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release_median[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release_quantile_low[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release_quantile_high[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot <- d.plot %>%
  mutate(strategy = ordered(strategy, levels = c("norelease", "bothsex", "maleonly", "bloodfed")))

f1.2 <- ggplot(data = d.plot) +
  # Cases attributed to wild mosquitoes (overlay on top of the total cases)
  geom_line(aes(x = generation, y = wild_mean, color = strategy),
            linetype = "dashed", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = wild_mean, color = strategy),
             size = 2, shape = 1, position = pd) +

  # Total cases
  geom_errorbar(aes(x = generation, ymin = total_quantile_low, ymax = total_quantile_high,
                    color = strategy),
                width = 0.3, position = pd) +
  geom_line(aes(x = generation, y = total_mean, color = strategy),
            linetype = "solid", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = total_mean, color = strategy),
             size = 2, shape = 19, position = pd) +

  scale_x_continuous(name="Generation of releases", breaks = 1:10) +
  scale_color_manual(name = "Release strategies", values = col.scenarios) +

  # Cases attributed to released mosquitoes
  geom_errorbar(aes(x = generation, ymin = release_quantile_low, ymax = release_quantile_high,
                    color = strategy),
                width = 0.3, position = pd) +
  geom_line(aes(x = generation, y = release_mean, color = strategy),
            linetype = "dotted", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = release_mean, color = strategy),
             size = 2, shape = 8, position = pd) +

  ylab("Number of new dengue infections") +
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
    legend.text = element_text(size = 10)
  )
print(f1.2)
ggsave(filename = paste0("Baseline_", index, "_generation_mean and 95 percent CI_v2.png"), 
       plot = f1.2,
       width = 9, height = 5, units = "in", dpi = 320)




# Figure 2: total number of cases by releasing strategies -----------------

# Total number of cases for the 10 generations:
d.sum <- d.case %>% 
  # select(iteration, strategy, generation, release, wild, total) %>%
  # melt(id.var = c("iteration", "strategy", "generation")) %>% 
  group_by(strategy, iteration) %>% 
  summarise(release_sum = sum(release),
            wild_sum = sum(wild),
            total_sum = sum(total)) %>%
  mutate(sum_increase = ifelse(wild_sum > 0, release_sum / wild_sum, 0))

# Summarize over the 100 iterations
d.sum.summary <- d.sum %>% 
  melt(id.var = c("strategy", "iteration")) %>%
  group_by(strategy, variable) %>%
  summarise(mean = mean(value), median = median(value), 
            sd = sd(value), IQR = IQR(value),
            quantile.low = quantile(value, 0.025), 
            quantile.high = quantile(value, 0.975))

# Make a plot
d.sum.plot <- d.sum.summary[!grepl(pattern = "sum_", x = d.sum.summary$variable), ]
d.sum.plot$strategy <- ordered(d.sum.plot$strategy, 
                               levels = c("norelease", "bothsex", "maleonly", "bloodfed"))
d.sum.plot$variable <- ordered(d.sum.plot$variable, 
                               levels = c("total_sum", "wild_sum", "release_sum"))

sec.axis.scale <- d.sum.plot$mean[d.sum.plot$strategy == "norelease" & d.sum.plot$variable == "total_sum"]
pd <- position_dodge(0.8)

f2.1 <- ggplot(data = d.sum.plot, 
             aes(x = strategy, y = mean, 
                 color = strategy, 
                 linetype = variable,
                 shape = variable, 
                 group = interaction(variable, strategy))) + 
  geom_hline(yintercept = sec.axis.scale, color = "gray90", linetype = "longdash") + 
  geom_hline(yintercept = 0, color = "gray90", linetype = "longdash") + 
  geom_errorbar(aes(ymin = quantile.low, ymax = quantile.high), 
                width = 0.3, position = pd) + 
  geom_point(size = 3, position = pd) +
  scale_x_discrete(name="Release strategies", 
                   labels = c("no release", "un-fed", "only male", "blood-fed")) + 
  scale_y_continuous(name = "Total number of infections in ten generations",
                     sec.axis = sec_axis(~ . / sec.axis.scale, 
                                         name = "",
                                         breaks = seq(0, 1, 0.2),
                                         labels = scales::percent)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("no release", "un-fed", "only male", "blood-fed"),
                     guide = 'none') + 
  scale_shape_manual(name = "Attribution", values = c(19, 1, 8),
                     labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
  scale_linetype_manual(name = "Attribution", values = c(1, 2, 3),
                        labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
  ylab("Total number of infections in ten generations") + 
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
    axis.text.x = element_text(size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(2, "line")
  )
print(f2.1)
ggsave(filename = paste0("Baseline_", index, "_total_mean and 95 percent CI_v1.png"),
       plot = f2.1, 
       width = 7, height = 5, units = "in", dpi = 320)


# Only the proportion relative to the no-release
d.sum.plot.proportion <- d.sum.plot %>%
  mutate(mean = mean / sec.axis.scale,
         median = median / sec.axis.scale,
         sd = sd / sec.axis.scale,
         IQR = IQR / sec.axis.scale,
         quantile.low = quantile.low / sec.axis.scale,
         quantile.high = quantile.high / sec.axis.scale)

f2.2 <- ggplot(data = d.sum.plot.proportion, 
               aes(x = strategy, y = mean, 
                   color = strategy, 
                   linetype = variable,
                   shape = variable, 
                   group = interaction(variable, strategy))) + 
  geom_errorbar(aes(ymin = quantile.low, ymax = quantile.high), 
                width = 0.3, position = pd) + 
  geom_point(size = 3, position = pd) +
  scale_x_discrete(name="Release strategies", 
                   labels = c("no release", "un-fed", "only male", "blood-fed")) + 
  scale_y_continuous(name = "Total infections relative to no release",
                     breaks = seq(0, 1, 0.2)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("no release", "un-fed", "only male", "blood-fed"),
                     guide = 'none') + 
  scale_shape_manual(name = "Attribution", values = c(19, 1, 8),
                     labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
  scale_linetype_manual(name = "Attribution", values = c(1, 2, 3),
                        labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
  ylab("Total number of infections in ten generations") + 
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
    axis.text.x = element_text(size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(2, "line")
  )
print(f2.2)
ggsave(filename = paste0("Baseline_", index, "_total_mean and 95 percent CI_v2.png"),
       plot = f2.2, 
       width = 6, height = 5, units = "in", dpi = 320)


save.image(paste0("Baseline_", index, "_visualiation.RData"))
save(f1.1, f1.2, f2, file = paste0("Baseline_", index, "_figures.RData"))
