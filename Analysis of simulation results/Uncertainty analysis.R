### Uncertainty analysis
### 2020.9.8
### Input data: results from the python simulation

library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# load("Uncertainty_analysis.RData")



# Data management ---------------------------------------------------------

# Paths to the simulation output:
bloodfed.dir <- "../Bloodfed_uncertainty/"
bothsex.dir <- "../Bothsex_uncertainty/"
maleonly.dir <- "../Maleonly_uncertainty/"
norelease.dir <- "../Norelease_uncertainty/"

load("Uncertainty_parameters.RData")



# The global list for all runs:
results.list <- vector(mode = "list", length = 1001)

# Some classification labels:
strategies <- c("bloodfed", "bothsex", "maleonly", "norelease")
outcomes <- c("case.wild", "case.release", "mosq.wild", "mosq.release")
metrics <- c("cases", "cases", "mosquitoes", "mosquitoes")
populations <- c("wild", "release", "wild", "release")


# Read in the results for all runs (n = 1000):
for(s in 1:1001){
  ss <- sprintf("%03d", s)
  
  results <- list(
    bloodfed = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_Gens_", ss, ".csv"),
                        header = FALSE, stringsAsFactors = FALSE),
    bothsex = read.csv(paste0(bothsex.dir, "simulation_Bothsex_UA_Gens_", ss, ".csv"),
                       header = FALSE, stringsAsFactors = FALSE),
    maleonly = read.csv(paste0(maleonly.dir, "simulation_Maleonly_UA_Gens_", ss, ".csv"),
                        header = FALSE, stringsAsFactors = FALSE),
    norelease = read.csv(paste0(norelease.dir, "simulation_NoRelease_UA_Gens_", ss, ".csv"),
                         header = FALSE, stringsAsFactors = FALSE)
  )
  
  n <- nrow(results$bloodfed)
  

  
  # Separating the different output metrics:
  d <- vector(mode = "list", length = 4 * 4)
  names(d) <- paste(rep(strategies, each = 4), rep(outcomes, times = 4), sep = ".")
  
  for(i in strategies){
    for(j in 1:4){
      temp.name <- paste(i, outcomes[j], sep = ".")
      d[[temp.name]] <- results[[i]][, (1:10 + 10*(j-1))] %>%
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
    mutate(total = wild + release)
  
  
  # Summarize means:
  d.summary <- d.case %>% group_by(strategy, generation) %>%
    summarise(wild_mean = mean(wild), 
              release_mean = mean(release), 
              total_mean = mean(total)) 
  
  
  # Convert means as the proportions of the cases transmitted by the wild mosquitoes at the first generation 
  d.summary <- d.summary %>%
    mutate(wild_mean_prop = -1, 
           release_mean_prop = -1,
           total_mean_prop = -1)
  for(i in strategies){
    baseline_mean <- as.numeric(d.summary[d.summary$strategy == i & d.summary$generation == 1, "wild_mean"])
    d.summary[d.summary$strategy == i, c("wild_mean_prop", "release_mean_prop", "total_mean_prop")] <- 
      d.summary[d.summary$strategy == i, c("wild_mean", "release_mean", "total_mean")] / baseline_mean
  }
  d.summary$mean_increase <- d.summary$release_mean / d.summary$wild_mean
  d.summary$run <- s
  
  results.list[[s]] <- d.summary
  
}

# Combine all results:
d.ua <- results.list[2:1001] %>% bind_rows()


# Summarize over all 1000 iterations
d.ua.summary <- d.ua %>%
  group_by(strategy, generation) %>%
  summarise(wild.mean = mean(wild_mean), 
            wild.median = median(wild_mean), 
            wild.sd = sd(wild_mean), 
            wild.IQR = IQR(wild_mean),
            wild.quantile_low = quantile(wild_mean, 0.025), 
            wild.quantile_high = quantile(wild_mean, 0.975),
            release.mean = mean(release_mean), 
            release.median = median(release_mean), 
            release.sd = sd(release_mean), 
            release.IQR = IQR(release_mean),
            release.quantile_low = quantile(release_mean, 0.025), 
            release.quantile_high = quantile(release_mean, 0.975),
            total.mean = mean(total_mean), 
            total.median = median(total_mean), 
            total.sd = sd(total_mean), 
            total.IQR = IQR(total_mean),
            total.quantile_low = quantile(total_mean, 0.025), 
            total.quantile_high = quantile(total_mean, 0.975),
            
            wild.prop.mean = mean(wild_mean_prop), 
            wild.prop.median = median(wild_mean_prop), 
            wild.prop.sd = sd(wild_mean_prop), 
            wild.prop.IQR = IQR(wild_mean_prop),
            wild.prop.quantile_low = quantile(wild_mean_prop, 0.025), 
            wild.prop.quantile_high = quantile(wild_mean_prop, 0.975),
            release.prop.mean = mean(release_mean_prop), 
            release.prop.median = median(release_mean_prop), 
            release.prop.sd = sd(release_mean_prop), 
            release.prop.IQR = IQR(release_mean_prop),
            release.prop.quantile_low = quantile(release_mean_prop, 0.025), 
            release.prop.quantile_high = quantile(release_mean_prop, 0.975),
            total.prop.mean = mean(total_mean_prop), 
            total.prop.median = median(total_mean_prop), 
            total.prop.sd = sd(total_mean_prop), 
            total.prop.IQR = IQR(total_mean_prop),
            total.prop.quantile_low = quantile(total_mean_prop, 0.025), 
            total.prop.quantile_high = quantile(total_mean_prop, 0.975),
            
            increase.mean = mean(mean_increase), 
            increase.median = median(mean_increase), 
            increase.sd = sd(mean_increase), 
            increase.IQR = IQR(mean_increase),
            increase.quantile_low = quantile(mean_increase, 0.025), 
            increase.quantile_high = quantile(mean_increase, 0.975))



# Visualization -----------------------------------------------------------

# shared options:
col.scenarios <- c("gray70", brewer.pal(n = 3, name = "Set2")[c(2, 3, 1)])
pd <- position_dodge(0.5)
et <- element_text(size = 14)



# Figure 1: the change of new dengue infections through generations -------

d.plot1 <- d.ua.summary %>%
  select(generation, strategy, total.prop.mean, total.prop.quantile_low, total.prop.quantile_high) %>%
  rename(mean = total.prop.mean, quantile_low = total.prop.quantile_low, quantile_high = total.prop.quantile_high) %>%
  mutate(attribution = "total")

d.plot2 <- d.ua.summary %>%
  select(generation, strategy, wild.prop.mean, wild.prop.quantile_low, wild.prop.quantile_high) %>%
  rename(mean = wild.prop.mean, quantile_low = wild.prop.quantile_low, quantile_high = wild.prop.quantile_high) %>%
  mutate(quantile_low = NA, quantile_high = NA) %>%
  filter(strategy %in% c("bloodfed", "bothsex")) %>%
  mutate(attribution = "wild\nmosquitoes")
# d.plot2$mean[d.plot2$strategy %in% c("maleonly", "norelease")] <- NA

d.plot3 <- d.ua.summary %>%
  select(generation, strategy, release.prop.mean, release.prop.quantile_low, release.prop.quantile_high) %>%
  rename(mean = release.prop.mean, quantile_low = release.prop.quantile_low, quantile_high = release.prop.quantile_high) %>%
  # mutate(quantile_low = NA, quantile_high = NA) %>%
  filter(strategy %in% c("bloodfed", "bothsex")) %>%
  mutate(attribution = "released\nmosquitoes")

d.ua.plot <- bind_rows(d.plot1, d.plot2, d.plot3) %>% 
  mutate(strategy = ordered(strategy, 
                            levels = c("norelease", "bothsex", "maleonly", "bloodfed")),
         attribution = ordered(attribution, 
                               levels = c("total", "wild\nmosquitoes", "released\nmosquitoes")))

f1.1 <- ggplot(data = d.ua.plot, aes(x = generation, y = mean, 
                                  color = strategy, 
                                  linetype = attribution,
                                  shape = attribution, 
                                  group = interaction(attribution, strategy))) + 
  # Total cases
  geom_errorbar(aes(ymin = quantile_low, ymax = quantile_high), 
                width = 0.3, position = pd) + 
  geom_line(size = 0.6, position = pd) +
  geom_point(size = 2, position = pd) +
  
  scale_x_continuous(name = "Generation of releases", breaks = 1:10) + 
  scale_y_continuous(name = "Mean proportion of dengue infections", breaks = seq(0, 1, 0.2)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("no release", "un-fed", "only male", "blood-fed")) + 
  scale_shape_manual(name = "Attribution", values = c(19, 1, 8)) + 
  scale_linetype_manual(name = "Attribution", values = c(1, 2, 3)) + 
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
ggsave(filename = "UA_generation_proportion_mean and 95 percent CI_v1.png", plot = f1.1, 
       width = 9, height = 5, units = "in", dpi = 320)



# Alternative Figure 1:
d.plot <- d.ua.summary
d.plot$release.prop.mean[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release.prop.median[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release.prop.sd[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release.prop.IQR[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release.prop.quantile_low[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot$release.prop.quantile_high[d.plot$strategy %in% c("maleonly", "norelease")] <- NA
d.plot <- d.plot %>%
  mutate(strategy = ordered(strategy, levels = c("norelease", "bothsex", "maleonly", "bloodfed")))

f1.2 <- ggplot(data = d.plot) +
  # Cases attributed to wild mosquitoes (overlay on top of the total cases)
  geom_line(aes(x = generation, y = wild.prop.mean, color = strategy),
            linetype = "dashed", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = wild.prop.mean, color = strategy),
             size = 2, shape = 1, position = pd) +
  
  # Total cases
  geom_errorbar(aes(x = generation, ymin = total.prop.quantile_low, ymax = total.prop.quantile_high,
                    color = strategy),
                width = 0.3, position = pd) +
  geom_line(aes(x = generation, y = total.prop.mean, color = strategy),
            linetype = "solid", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = total.prop.mean, color = strategy),
             size = 2, shape = 19, position = pd) +
  
  scale_x_continuous(name = "Generation of releases", breaks = 1:10) +
  scale_y_continuous(name = "Mean proportion of dengue infections", breaks = seq(0, 1, 0.2)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios) +
  
  # Cases attributed to released mosquitoes
  geom_errorbar(aes(x = generation, ymin = release.prop.quantile_low, ymax = release.prop.quantile_high,
                    color = strategy),
                width = 0.3, position = pd) +
  geom_line(aes(x = generation, y = release.prop.mean, color = strategy),
            linetype = "dotted", size = 0.6, position = pd) +
  geom_point(aes(x = generation, y = release.prop.mean, color = strategy),
             size = 2, shape = 8, position = pd) +

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
ggsave(filename = "UA_generation_proportion_mean and 95 percent CI_v2.png", plot = f1.2,
       width = 9, height = 5, units = "in", dpi = 320)





# Figure 2: total number of cases by releasing strategies -----------------

# Total number of cases for the 10 generations:
d.ua.sum <- d.ua %>% 
  group_by(strategy, run) %>%
  summarise(wild_sum_mean = sum(wild_mean),
            release_sum_mean = sum(release_mean),
            total_sum_mean = sum(total_mean))

# Calculate the proportion relative to the no-release strategy:
ua.norelease.wild <- rep(d.ua.sum$wild_sum_mean[d.ua.sum$strategy == "norelease"], n = 4)

d.ua.sum <- d.ua.sum %>%
  mutate(wild_sum_mean_prop = wild_sum_mean / ua.norelease.wild, 
         release_sum_mean_prop = release_sum_mean / ua.norelease.wild,
         total_sum_mean_prop = total_sum_mean / ua.norelease.wild,
         sum_mean_increase = release_sum_mean / wild_sum_mean)
d.ua.sum$strategy <- ordered(d.ua.sum$strategy, 
                             levels = c("norelease", "bothsex", "maleonly", "bloodfed"))

rm(ua.norelease.wild)


# Summarize over the 1000 runs
d.ua.sum.summary <- d.ua.sum %>% 
  melt(id.var = c("strategy", "run")) %>%
  group_by(strategy, variable) %>%
  summarise(mean = mean(value), median = median(value), 
            sd = sd(value), IQR = IQR(value),
            quantile.low = quantile(value, 0.025), 
            quantile.high = quantile(value, 0.975))


### Make a plot
# Prepare the data:
d.ua.sum.plot1 <- d.ua.sum.summary %>%
  filter(grepl(pattern = "_prop", x = variable),
         strategy != "norelease"
         # , !(strategy == "maleonly" & variable %in% c("wild_sum_mean_prop", "release_sum_mean_prop"))
         ) %>%
  mutate(strategy = ordered(strategy, 
                            levels = c("bothsex", "maleonly", "bloodfed")),
         variable = ordered(variable, 
                            levels = c("total_sum_mean_prop", "wild_sum_mean_prop", "release_sum_mean_prop")))

d.ua.sum.plot2 <- d.ua.sum %>% 
  melt(id.var = c("strategy", "run")) %>%
  filter(grepl(pattern = "_prop", x = variable),
         strategy != "norelease"
         # , !(strategy == "maleonly" & variable %in% c("wild_sum_mean_prop", "release_sum_mean_prop"))
         ) %>%
  mutate(strategy = ordered(strategy, 
                            levels = c("bothsex", "maleonly", "bloodfed")),
         variable = ordered(variable, 
                            levels = c("total_sum_mean_prop", "wild_sum_mean_prop", "release_sum_mean_prop")))

# Making the figure
col.scenarios <- brewer.pal(n = 3, name = "Set2")[c(2, 3, 1)]
pd <- position_dodge(0.8)
et <- element_text(size = 14)

f2 <- ggplot() + 
  geom_violin(data = d.ua.sum.plot2,
              aes(x = strategy, y = value,
                  color = strategy, 
                  linetype = variable,
                  group = interaction(variable, strategy)),
              draw_quantiles = c(0.025, 0.5, 0.975),
              trim = TRUE, scale = "width",
              position = pd) + 
  # geom_point(data = d.ua.sum.plot2,
  #            aes(x = strategy, y = value,
  #                color = strategy, 
  #                shape = variable,
  #                group = interaction(variable, strategy)),
  #            size = 0.1, position = pd) + 
  # geom_errorbar(data = d.ua.sum.plot1,
  #               aes(x = strategy, 
  #                   ymin = quantile.low, 
  #                   ymax = quantile.high,
  #                   color = strategy, 
  #                   shape = variable, 
  #                   group = interaction(variable, strategy)), 
  #               width = 0.3, position = pd) + 
  geom_point(data = d.ua.sum.plot1, 
             aes(x = strategy, y = mean, 
                 color = strategy, 
                 shape = variable, 
                 group = interaction(variable, strategy)),
             size = 2, position = pd) + 
  scale_x_discrete(name="Release strategies", 
                   labels = c("un-fed", "only male", "blood-fed")) + 
  scale_y_continuous(name = "Total infections relative to no release",
                     breaks = seq(0, 1, 0.2)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("no release", "un-fed", "only male", "blood-fed"),
                     guide = 'none') + 
  scale_shape_manual(name = "Attribution", values = c(19, 1, 8),
                     labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
  scale_linetype_manual(name = "Attribution", values = c(1, 2, 3),
                        labels = c("total", "wild\nmosquitoes", "released\nmosquitoes")) + 
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
print(f2)
ggsave(filename = "UA_total_proportion_mean and 95 percent CI_v2.png",
       plot = f2, width = 6, height = 5, units = "in", dpi = 320)





# Effect of each parameters -----------------------------------------------

# Join the parameter matrix with the results matrix:
d.ua.sum.metric <- d.ua.sum %>%
  select(run, strategy, total_sum_mean_prop, release_sum_mean_prop) %>%
  # melt(id.var = c("name", "strategy")) %>%
  # dcast(name ~ strategy + variable, value.var = "value") %>%
  filter(strategy != "norelease")

pe.ua <- bind_rows(Uncertainty_parameters_bloodfed[-1,],
                   Uncertainty_parameters_bothsex[-1, ],
                   Uncertainty_parameters_maleonly[-1, ]) %>%
  bind_cols(d.ua.sum.metric)
print(identical(pe.ua$name, pe.ua$run))



# Figure 3.1 Parameter effect on total infection --------------------------
pe.ua.long <- pe.ua %>%
  select(ovi, bitephase, EIP, mubites, 
         a, b, s, c, 
         sdovi, sdhost, sdEIP,  
         run, strategy, 
         total_sum_mean_prop, release_sum_mean_prop) %>%
  melt(id.var = c("run", "strategy", 
                  "total_sum_mean_prop", "release_sum_mean_prop")) %>%
  mutate(variable = ordered(as.character(variable), 
                            levels = c("ovi", "sdovi", 
                                       "bitephase", "sdhost", 
                                       "mubites", 
                                       "EIP", "sdEIP", 
                                       "a", "b", "s", "c"),
                            labels = c("mean oviposition duration",
                                       "SD of oviposition duration",
                                       "mean blood-feeding duration",
                                       "SD of blood-feeding duration",
                                       "mean bites per cycle",
                                       "mean EIP", 
                                       "SD of EIP",
                                       "Hazard parameter: a",
                                       "Hazard parameter: b",
                                       "Hazard parameter: s",
                                       "Hazard parameter: c")))

f3.1 <- ggplot(data = pe.ua.long, 
             aes(x = value, y = total_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 4, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Total infections relative to no release",
                     breaks = seq(0, 1, 0.1)) + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("un-fed", "only male", "blood-fed")) + 
  scale_fill_manual(name = "Release strategies", 
                    values = alpha(col.scenarios, alpha = 0.5),
                    labels = c("un-fed", "only male", "blood-fed"),
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
    axis.title.x = element_blank(),
    axis.text.x = et,
    legend.title = et,
    legend.text = et,
    legend.position = "bottom", 
    strip.background = element_blank(),
    strip.text.x = et,
    strip.placement = "outside"
  )
print(f3.1)
ggsave(filename = "UA_parameter_effect_total.png", plot = f3.1,
       width = 10, height = 12, units = "in", dpi = 320)
ggsave(filename = "UA_parameter_effect_total.pdf", plot = f3.1,
       width = 10, height = 12, units = "in", dpi = 320)


# Figure: number of infections contributed by the released mosquitoes
f3.1.r <- ggplot(data = pe.ua.long, 
               aes(x = value, y = release_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 4, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Total infections relative to no release") + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("un-fed", "only male", "blood-fed")) + 
  scale_fill_manual(name = "Release strategies", 
                    values = alpha(col.scenarios, alpha = 0.5),
                    labels = c("un-fed", "only male", "blood-fed"),
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
    axis.title.x = element_blank(),
    axis.text.x = et,
    legend.title = et,
    legend.text = et,
    legend.position = "bottom", 
    strip.background = element_blank(),
    strip.text.x = et,
    strip.placement = "outside"
  )
print(f3.1.r)
ggsave(filename = "UA_parameter_effect_released.png", plot = f3.1.r,
       width = 10, height = 12, units = "in", dpi = 320)
ggsave(filename = "UA_parameter_effect_released.pdf", plot = f3.1.r,
       width = 10, height = 12, units = "in", dpi = 320)

save.image("Uncertainty_analysis.RData")
