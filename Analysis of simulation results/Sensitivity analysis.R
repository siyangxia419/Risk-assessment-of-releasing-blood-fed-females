### Sensitivity analysis
### 2020.10.14
### Input data: results from the python simulation

library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(randomForest)
library(mmpf)

# load("Sensitivity_analysis.RData")




# Data management ---------------------------------------------------------

# Paths to the simulation output:
bloodfed.dir <- "../Bloodfed_sensitivity_adjusted_feeding/"
bothsex.dir <- "../Bothsex_sensitivity/"
maleonly.dir <- "../Maleonly_sensitivity/"
norelease.dir <- "../Norelease_sensitivity/"

load("Sensitivity_parameters.RData")



# The global list for all runs:
N <- nrow(Sensitivity_parameters_bloodfed)

# Some classification labels:
strategies <- c("bloodfed", "bothsex", "maleonly", "norelease")
outcomes <- c("case.wild", "case.release", "mosq.wild", "mosq.release")
metrics <- c("cases", "cases", "mosquitoes", "mosquitoes")
populations <- c("wild", "release", "wild", "release")


# List all result files for each releasing strategy:
bloodfed.file <- list.files(path = bloodfed.dir, pattern = "simulation_BloodfedFemales_SA_Gens_*")
bothsex.file <- list.files(path = bothsex.dir, pattern = "simulation_Bothsex_SA_Gens_*")
maleonly.file <- list.files(path = maleonly.dir, pattern = "simulation_Maleonly_SA_Gens_*")
norelease.file <- list.files(path = norelease.dir, pattern = "simulation_NoRelease_SA_Gens_*")

# check which runs failed:
bloodfed.index <- as.numeric(gsub(pattern = "simulation_BloodfedFemales_SA_Gens_", replacement = "", 
                                  gsub(pattern = ".csv", replacement = "", x = bloodfed.file)))
bloodfed.index.fail <- which(!(1:N %in% bloodfed.index))

bothsex.index <- as.numeric(gsub(pattern = "simulation_Bothsex_SA_Gens_", replacement = "", 
                                  gsub(pattern = ".csv", replacement = "", x = bothsex.file)))
bothsex.index.fail <- which(!(1:N %in% bothsex.index))

maleonly.index <- as.numeric(gsub(pattern = "simulation_Maleonly_SA_Gens_", replacement = "", 
                                 gsub(pattern = ".csv", replacement = "", x = maleonly.file)))
maleonly.index.fail <- which(!(1:N %in% maleonly.index))

norelease.index <- as.numeric(gsub(pattern = "simulation_NoRelease_SA_Gens_", replacement = "", 
                                 gsub(pattern = ".csv", replacement = "", x = norelease.file)))
norelease.index.fail <- which(!(1:N %in% norelease.index))

index.fail <- unique(c(bloodfed.index.fail,
                       bothsex.index.fail,
                       maleonly.index.fail,
                       norelease.index.fail))
index.successful <- (1:N)[-(index.fail)]




# Read in the results for all runs (n = 1000):
results.list <- vector(mode = "list", length = length(index.successful))

for(s in index.successful){
  ss <- sprintf("%03d", s)
  
  results <- list(
    bloodfed = read.csv(paste0(bloodfed.dir, "simulation_BloodfedFemales_SA_Gens_", ss, ".csv"),
                        header = FALSE, stringsAsFactors = FALSE),
    bothsex = read.csv(paste0(bothsex.dir, "simulation_Bothsex_SA_Gens_", ss, ".csv"),
                       header = FALSE, stringsAsFactors = FALSE),
    maleonly = read.csv(paste0(maleonly.dir, "simulation_Maleonly_SA_Gens_", ss, ".csv"),
                        header = FALSE, stringsAsFactors = FALSE),
    norelease = read.csv(paste0(norelease.dir, "simulation_NoRelease_SA_Gens_", ss, ".csv"),
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
d.sa <- results.list %>% bind_rows()

# Separate the purpose of the sensitivity matrix:
sa.aim <- read.csv("GSM_parameters_for_female_release_model_full.csv",
                   header = TRUE, stringsAsFactors = FALSE)
identical(sa.aim$Index, Sensitivity_parameters_bloodfed$name)
Sensitivity_parameters_bloodfed$aim <- sa.aim$Aim
Sensitivity_parameters_bothsex$aim <- sa.aim$Aim
Sensitivity_parameters_maleonly$aim <- sa.aim$Aim
Sensitivity_parameters_norelease$aim <- sa.aim$Aim







# # Summarize over all 1000 iterations
# d.sa.summary <- d.sa %>%
#   group_by(strategy, generation) %>%
#   summarise(wild.mean = mean(wild_mean), 
#             wild.median = median(wild_mean), 
#             wild.sd = sd(wild_mean), 
#             wild.IQR = IQR(wild_mean),
#             wild.quantile_low = quantile(wild_mean, 0.025), 
#             wild.quantile_high = quantile(wild_mean, 0.975),
#             release.mean = mean(release_mean), 
#             release.median = median(release_mean), 
#             release.sd = sd(release_mean), 
#             release.IQR = IQR(release_mean),
#             release.quantile_low = quantile(release_mean, 0.025), 
#             release.quantile_high = quantile(release_mean, 0.975),
#             total.mean = mean(total_mean), 
#             total.median = median(total_mean), 
#             total.sd = sd(total_mean), 
#             total.IQR = IQR(total_mean),
#             total.quantile_low = quantile(total_mean, 0.025), 
#             total.quantile_high = quantile(total_mean, 0.975),
#             
#             wild.prop.mean = mean(wild_mean_prop), 
#             wild.prop.median = median(wild_mean_prop), 
#             wild.prop.sd = sd(wild_mean_prop), 
#             wild.prop.IQR = IQR(wild_mean_prop),
#             wild.prop.quantile_low = quantile(wild_mean_prop, 0.025), 
#             wild.prop.quantile_high = quantile(wild_mean_prop, 0.975),
#             release.prop.mean = mean(release_mean_prop), 
#             release.prop.median = median(release_mean_prop), 
#             release.prop.sd = sd(release_mean_prop), 
#             release.prop.IQR = IQR(release_mean_prop),
#             release.prop.quantile_low = quantile(release_mean_prop, 0.025), 
#             release.prop.quantile_high = quantile(release_mean_prop, 0.975),
#             total.prop.mean = mean(total_mean_prop), 
#             total.prop.median = median(total_mean_prop), 
#             total.prop.sd = sd(total_mean_prop), 
#             total.prop.IQR = IQR(total_mean_prop),
#             total.prop.quantile_low = quantile(total_mean_prop, 0.025), 
#             total.prop.quantile_high = quantile(total_mean_prop, 0.975),
#             
#             increase.mean = mean(mean_increase), 
#             increase.median = median(mean_increase), 
#             increase.sd = sd(mean_increase), 
#             increase.IQR = IQR(mean_increase),
#             increase.quantile_low = quantile(mean_increase, 0.025), 
#             increase.quantile_high = quantile(mean_increase, 0.975))






# Summerizing the results ---------------------------------------------

# Total number of cases for the 10 generations:
d.sa.sum <- d.sa %>% 
  group_by(strategy, run) %>%
  summarise(wild_sum_mean = sum(wild_mean),
            release_sum_mean = sum(release_mean),
            total_sum_mean = sum(total_mean))

# Calculate the proportion relative to the no-release strategy:
sa.norelease.wild <- rep(d.sa.sum$wild_sum_mean[d.sa.sum$strategy == "norelease"], 4)

d.sa.sum$wild_sum_mean_prop <- d.sa.sum$wild_sum_mean / sa.norelease.wild
d.sa.sum$release_sum_mean_prop <- d.sa.sum$release_sum_mean / sa.norelease.wild
d.sa.sum$total_sum_mean_prop <- d.sa.sum$total_sum_mean / sa.norelease.wild
d.sa.sum$sum_mean_increase = d.sa.sum$release_sum_mean / d.sa.sum$wild_sum_mean
d.sa.sum$strategy <- ordered(d.sa.sum$strategy, 
                             levels = c("norelease", "bothsex", "maleonly", "bloodfed"))

rm(sa.norelease.wild)


# Join the parameter matrix with the results matrix:
d.sa.sum.metric <- d.sa.sum %>%
  select(run, strategy, total_sum_mean_prop, release_sum_mean_prop) %>%
  filter(strategy != "norelease")


pe.sa <- bind_rows(Sensitivity_parameters_bloodfed[Sensitivity_parameters_bloodfed$name %in% index.successful, ],
                   Sensitivity_parameters_bothsex[Sensitivity_parameters_bothsex$name %in% index.successful, ],
                   Sensitivity_parameters_maleonly[Sensitivity_parameters_maleonly$name %in% index.successful, ]) %>%
  bind_cols(d.sa.sum.metric)
pe.sa$aim <- gsub(pattern = "LSA_", replacement = "", x = pe.sa$aim)
print(identical(pe.sa$name, pe.sa$run))





# Figure 2.1 One-way sensitivity analysis ---------------------------------


lsa.parameters <- c("hi", "FF", 
                    "vcw", "sdVCw", "vcr", "sdVCr", 
                    "releaseratio", "fed")
lsa.sa <- pe.sa[pe.sa$aim %in% lsa.parameters, ]

lsa.sa$releaseratio <- lsa.sa$releaseratio * 0.1

lsa.sa$value <- NA
for(i in 1:nrow(lsa.sa)) lsa.sa$value[i] <- lsa.sa[i, lsa.sa$aim[i]]

lsa.sa$variable <- ordered(as.character(lsa.sa$aim),
                           levels = lsa.parameters,
                           labels = c("Proportion of infectious human",
                                      "Age of mosquitoes at release",
                                      "mean VC of wild mosquitoes",
                                      "SD of VC of wild mosquitoes",
                                      "mean VC of released mosquitoes",
                                      "SD of VC of released mosquitoes",
                                      "Relative size of release",
                                      "Proportion of females fed"))


# Plotting:
# shared options:
col.scenarios <- brewer.pal(n = 3, name = "Set2")[c(2, 3, 1)]
pd <- position_dodge(0.8)
et <- element_text(size = 14)

# Figure: total number of infections
f2.1 <- ggplot(data = lsa.sa, 
               aes(x = value, y = total_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Relative total infections") + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("un-fed", "only male", "blood-fed")) + 
  scale_fill_manual(name = "Release strategies", 
                    values = alpha(col.scenarios, alpha = 0.5),
                    labels = c("un-fed", "only male", "blood-fed"),
                    guide = "none") + 
  # ylim(min(lsa.sa$total_sum_mean_prop), max(lsa.sa$total_sum_mean_prop)) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.minor.y = element_line(color = "gray95"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1.2, "lines"), 
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
print(f2.1)
ggsave(filename = "SA_one_way_parameter_effect_total.png", plot = f2.1,
       width = 11, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_one_way_parameter_effect_total.pdf", plot = f2.1,
       width = 11, height = 9, units = "in", dpi = 320)


# Figure: number of infections contributed by the releases mosquitoes
f2.1.r <- ggplot(data = lsa.sa, 
                 aes(x = value, y = release_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Relative total infections") + 
  scale_color_manual(name = "Release strategies", values = col.scenarios,
                     labels = c("un-fed", "only male", "blood-fed")) + 
  scale_fill_manual(name = "Release strategies", 
                    values = alpha(col.scenarios, alpha = 0.5),
                    labels = c("un-fed", "only male", "blood-fed"),
                    guide = "none") + 
  # ylim(min(lsa.sa$total_sum_mean_prop), max(lsa.sa$total_sum_mean_prop)) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.minor.y = element_line(color = "gray95"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1.2, "lines"), 
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
print(f2.1.r)
ggsave(filename = "SA_one_way_parameter_effect_released.png", plot = f2.1.r,
       width = 11, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_one_way_parameter_effect_released.pdf", plot = f2.1.r,
       width = 11, height = 9, units = "in", dpi = 320)






# Figure 2.2 GSA by parameter -------------------------------------------

gsa.sa <- pe.sa[pe.sa$aim == "GSA", ]

gsa.sa.long <- gsa.sa %>%
  select(hi, FF, vcw, vcr, releaseratio, sdVCw, sdVCr, fed,  
         run, strategy, 
         total_sum_mean_prop, release_sum_mean_prop) %>%
  mutate(releaseratio = releaseratio * 0.1) %>%
  melt(id.var = c("run", "strategy", 
                  "total_sum_mean_prop", "release_sum_mean_prop")) %>%
  mutate(variable = ordered(as.character(variable), 
                            levels = lsa.parameters,
                            labels = c("Proportion of infectious human",
                                       "Age of mosquitoes at release",
                                       "mean VC of wild mosquitoes",
                                       "SD of VC of wild mosquitoes",
                                       "mean VC of released mosquitoes",
                                       "SD of VC of released mosquitoes",
                                       "Relative size of release",
                                       "Proportion of females fed")))


# Figure: total number of infections
f2.2 <- ggplot(data = gsa.sa.long, 
               aes(x = value, y = total_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Relative total infections", 
                     breaks = seq(0, 3, 0.5)) + 
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
    panel.spacing = unit(1.2, "lines"), 
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
print(f2.2)
ggsave(filename = "SA_GSA_parameter_effect_total.png", plot = f2.2,
       width = 11, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_GSA_parameter_effect_total.pdf", plot = f2.2,
       width = 11, height = 9, units = "in", dpi = 320)


# Figure: total number of infections
f2.2.r <- ggplot(data = gsa.sa.long, 
               aes(x = value, y = release_sum_mean_prop, color = strategy)) + 
  facet_wrap(vars(variable), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_point(size = 0.05) +
  stat_smooth(aes(fill = strategy), inherit.aes = TRUE) + 
  scale_y_continuous(name = "Relative total infections", 
                     breaks = seq(0, 3, 0.5)) + 
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
    panel.spacing = unit(1.2, "lines"), 
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
print(f2.2.r)
ggsave(filename = "SA_GSA_parameter_effect_released.png", plot = f2.2.r,
       width = 11, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_GSA_parameter_effect_released.pdf", plot = f2.2.r,
       width = 11, height = 9, units = "in", dpi = 320)





# Figure 2.3 GSA by strategy ----------------------------------------------
# Summarize over all GSA runs
d.sa.sum.gsa <- d.sa.sum[(d.sa.sum$strategy != "norelease"), ]
identical(d.sa.sum.gsa$run, pe.sa$run)
d.sa.sum.gsa <- d.sa.sum.gsa[pe.sa$aim == "GSA", ]

d.sa.sum.summary <- d.sa.sum.gsa %>%
  melt(id.var = c("strategy", "run")) %>%
  group_by(strategy, variable) %>%
  summarise(mean = mean(value), median = median(value), 
            sd = sd(value), IQR = IQR(value),
            quantile.low = quantile(value, 0.025), 
            quantile.high = quantile(value, 0.975))


### Make a plot
# Prepare the data:
d.sa.sum.plot1 <- d.sa.sum.summary %>%
  filter(grepl(pattern = "_prop", x = variable),
         strategy != "norelease") %>%
  mutate(strategy = ordered(strategy, 
                            levels = c("bothsex", "maleonly", "bloodfed")),
         variable = ordered(variable, 
                            levels = c("total_sum_mean_prop", "wild_sum_mean_prop", "release_sum_mean_prop")))

d.sa.sum.plot2 <- d.sa.sum.gsa %>% 
  melt(id.var = c("strategy", "run")) %>%
  filter(grepl(pattern = "_prop", x = variable),
         strategy != "norelease") %>%
  mutate(strategy = ordered(strategy, 
                            levels = c("bothsex", "maleonly", "bloodfed")),
         variable = ordered(variable, 
                            levels = c("total_sum_mean_prop", "wild_sum_mean_prop", "release_sum_mean_prop")))


# Making the figure
f2.3 <- ggplot() + 
  geom_violin(data = d.sa.sum.plot2,
              aes(x = strategy, y = value,
                  color = strategy, 
                  linetype = variable,
                  group = interaction(variable, strategy)),
              draw_quantiles = c(0.025, 0.5, 0.975),
              trim = TRUE, scale = "width",
              position = pd) + 
  geom_point(data = d.sa.sum.plot1, 
             aes(x = strategy, y = mean, 
                 color = strategy, 
                 shape = variable, 
                 group = interaction(variable, strategy)),
             size = 2, position = pd) + 
  scale_x_discrete(name="Release strategies", 
                   labels = c("un-fed", "only male", "blood-fed")) + 
  scale_y_continuous(name = "Relative total infections",
                     breaks = seq(0, 3, 0.2)) + 
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
print(f2.3)
ggsave(filename = "SA_GSA_total_proportion_mean and 95 percent CI_v2.png",
       plot = f2.3, width = 6, height = 5, units = "in", dpi = 320)





# GSA analysis: compare different strategies ------------------------------

# Prepare the dataset: extract the key output ("total_sum_mean_prop") and key parameters:
compare.strategy.outcome <- d.sa.sum.gsa %>%
  select(strategy, run, total_sum_mean_prop) %>%
  dcast(formula = run ~ strategy, value.var = "total_sum_mean_prop")
compare.strategy.outcome$min_strategy <- apply(X = as.matrix(compare.strategy.outcome[, c("bothsex", "maleonly", "bloodfed")]), 
                                               MARGIN = 1,  
                                               FUN = function(x) c("bothsex", "maleonly", "bloodfed")[which.min(x)])
table(compare.strategy.outcome$min_strategy)

compare.strategy.parameter <- pe.sa %>%
  filter(aim == "GSA" & strategy == "bloodfed") %>%
  select(run, all_of(lsa.parameters))

identical(compare.strategy.outcome$run, compare.strategy.parameter$run)

compare.strategy <- bind_cols(compare.strategy.parameter, 
                              compare.strategy.outcome) %>%
  select(!starts_with("run"))



# a) Some graphical exploration:
compare.strategy.each.para <- compare.strategy %>% 
  mutate(releaseratio = releaseratio * 0.1) %>%
  melt(id.var = c("bothsex", "maleonly", "bloodfed", "min_strategy")) %>%
  mutate(variable = ordered(as.character(variable), 
                            levels = lsa.parameters,
                            labels = c("Proportion of infectious human",
                                       "Age of mosquitoes at release",
                                       "mean VC of wild mosquitoes",
                                       "SD of VC of wild mosquitoes",
                                       "mean VC of released mosquitoes",
                                       "SD of VC of released mosquitoes",
                                       "Relative size of release",
                                       "Proportion of females fed")))

f2.4 <- ggplot(data = compare.strategy.each.para, 
               aes(x = value, y = min_strategy, color = min_strategy)) + 
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975),
              trim = TRUE, scale = "count", adjust = 0.5) + 
  facet_wrap(vars(variable), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  scale_y_discrete(name="Optimal release strategy", 
                   labels = c("blood-fed", "only male")) + 
  xlab(NULL) + 
  scale_color_manual(name = "", values = col.scenarios[c(3, 2)],
                     labels = c("only male", "blood-fed"),
                     guide = 'none') + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    axis.line = element_line(colour = "black"),
    axis.title.y = et,
    axis.text.y = et,
    axis.title.x = et,
    axis.text.x = element_text(size = 14),
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = et
  )
print(f2.4)
ggsave(filename = "SA_GSA_compare_strategies_each_parameters.png",
       plot = f2.4,
       width = 13, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_GSA_compare_strategies_each_parameters.pdf", 
       plot = f2.4,
       width = 13, height = 9, units = "in", dpi = 320)




# Random forest analysis --------------------------------------------------

# Construct the independent and dependent variables:
compare.strategy.x <- as.matrix(compare.strategy[, lsa.parameters])
compare.strategy.y <- as.factor(compare.strategy$min_strategy)

# Tuning the random forest to find the optimal "mtry":
rf.tune <- tuneRF(x = compare.strategy.x, y = compare.strategy.y)
mtry.opt <- rf.tune[, "mtry"][which.min(rf.tune[, "OOBError"])]

# Random forest:
rf.strategy <- randomForest(x = compare.strategy.x, 
                            y = compare.strategy.y,
                            importance = TRUE, 
                            na.action = na.omit, 
                            mtry = mtry.opt)

# Parameter importance
vi.strategy <- importance(rf.strategy)
vi.strategy <- vi.strategy[order(-vi.strategy[, "MeanDecreaseAccuracy"]),]
print(vi.strategy)

# Further processing of the VI results:
vi.strategy <- as.data.frame(vi.strategy) %>%
  mutate(variable = as.character(row.names(vi.strategy)),
         label = ordered(variable, 
                         levels = lsa.parameters,
                         labels = c("Proportion of infectious human",
                                    "Age of mosquitoes at release",
                                    "mean VC of wild mosquitoes",
                                    "SD of VC of wild mosquitoes",
                                    "mean VC of released mosquitoes",
                                    "SD of VC of released mosquitoes",
                                    "Relative size of release",
                                    "Proportion of females fed"))) %>%
  mutate(label = as.character(label)) %>%
  arrange(MeanDecreaseAccuracy) %>%
  mutate(label = ordered(label, 
                         levels = label))

# Save the results:
write.csv(x = vi.strategy, file = "SA_variable important.csv", 
          quote = FALSE, row.names = TRUE)


# Make a bar plot:
f3 <- ggplot(data = vi.strategy) + 
  geom_col(aes(x = MeanDecreaseAccuracy, y = label)) + 
  xlab("Parameter importance metrics") + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
print(f3)

ggsave(filename = "SA_GSA_compare_strategies_parameter_importance.png", 
       plot = f3,
       width = 7, height = 4, units = "in", dpi = 320)
ggsave(filename = "SA_GSA_compare_strategies_parameter_importance.pdf", 
       plot = f3,
       width = 7, height = 4, units = "in", dpi = 320)



# Marginal effect of each parameter:
mp.strategy <- data.frame(parameter = character(),
                          par.value = numeric(),
                          bloodfed = numeric())
for(v in colnames(compare.strategy.x)){
  mp.var <- marginalPrediction(data = as.data.frame(compare.strategy.x),
                   vars = v, 
                   n = c(50, nrow(compare.strategy.x)), 
                   model = rf.strategy,
                   predict.fun = function(object, newdata) predict(object, newdata = newdata, type = "prob"))
  names(mp.var) <- c("par.value", "bloodfed", "maleonly")
  mp.var$parameter <- v
  mp.strategy <- rbind.data.frame(mp.strategy, as.data.frame(mp.var[, c("parameter", "par.value", "bloodfed")]))
}

mp.strategy$parameter <- ordered(as.character(mp.strategy$parameter),
                                 levels = lsa.parameters,
                                 labels = c("Proportion of infectious human",
                                            "Age of mosquitoes at release",
                                            "mean VC of wild mosquitoes",
                                            "SD of VC of wild mosquitoes",
                                            "mean VC of released mosquitoes",
                                            "SD of VC of released mosquitoes",
                                            "Relative size of release",
                                            "Proportion of females fed"))

# Plotting the marginal effect of each parameters:
f4 <- ggplot(data = mp.strategy, 
            aes(x = par.value, y = bloodfed)) + 
  facet_wrap(vars(parameter), nrow = 3, ncol = 3, 
             scales = "free", 
             strip.position = "bottom") + 
  geom_line(size = 1.2) +
  scale_y_continuous(name = "Predicted probability of favoring blood-fed females", 
                     breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.minor.y = element_line(color = "gray95"), 
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.minor.x = element_line(color = "gray95"),
    panel.spacing = unit(1.2, "lines"), 
    axis.line = element_line(colour = "black"),
    axis.title.y = et,
    axis.text.y = et,
    axis.title.x = element_blank(),
    axis.text.x = et,
    strip.background = element_blank(),
    strip.text.x = et,
    strip.placement = "outside"
  )
print(f4)
ggsave(filename = "SA_GSA_compare_strategies_marginal_effect.png", plot = f4,
       width = 11, height = 9, units = "in", dpi = 320)
ggsave(filename = "SA_GSA_compare_strategies_marginal_effect.pdf", plot = f4,
       width = 11, height = 9, units = "in", dpi = 320)


save.image("Sensitivity_analysis.RData")
