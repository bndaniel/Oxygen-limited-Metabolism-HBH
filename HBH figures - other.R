#R code for other figures in "Persistence of Oxygen Dependent Dissolved Organic Matter (ODDOM) induced by Oxygen-limited Metabolism in Hypoxic Marine Environments"

#### Lineweaver-Burk plot and Kms #####
setwd("~/Desktop/initial_resp_R/")
#Lineweaver-Burk plot
resp_data = read.csv("resp_rate_compounds.csv")


#new dataset with 1/V and 1/[s]: km using Lineweaver-Burk plot equation
inv_resp_data = resp_data
inv_resp_data$resp_rate = 1/inv_resp_data$resp_rate
inv_resp_data$average_O = 1/inv_resp_data$average_O
inv_resp_data

#subset to exclude controls and rsq<0.6 (need to go through and manually remove outliers?)
inv_resp_data = subset(inv_resp_data, bottle <= 10 & rsq >= 0.65)
inv_resp_data

# Define which compounds are ODDOM or LDOM
ldom_compounds <- c("glucose1", "glucose3","glucose4","fumaricacid1","fumaricacid2","tryptophan2", "pyruvate1")  # Add your ODDOM compounds here
oddom_compounds <- c("toluene1","toluene2","benzoate1","benzoate2","protochatechuate1","protochatechuate2","gentisate1","gentisate2")  # Add your LDOM compounds here

# Add a new column based on compound type
inv_resp_data <- inv_resp_data %>%
  mutate(DOM_Type = case_when(
    compound %in% oddom_compounds ~ "ODDOM",
    compound %in% ldom_compounds ~ "LDOM",
    TRUE ~ NA_character_  # Assigns NA if the compound is not listed
  ))
inv_resp_data

#took em out "tyrosine1","tyrosine2"
inv_resp_data = subset(inv_resp_data, inv_resp_data$DOM_Type != "NA")

#linear model of all 1/v vs 1/s
library(scales)

resp_plot <- ggplot(inv_resp_data, aes(x = average_O, y = resp_rate, color = DOM_Type)) + 
  geom_point(size=3) + 
  scale_x_continuous(name = "1/[O]") +
  scale_y_continuous(name = "1/V") +
  theme_classic() +
  labs(color = "DOM Type") +  # Change legend title
  scale_color_manual(values = c("LDOM" = "skyblue3", "ODDOM" = "firebrick3")) + # Define colors
  theme(legend.position = c(0.05, 0.95),  # Top-left position
        legend.justification = c(0, 1),  
        legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16),  
        legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )  

# Add both types of regression lines
resp_plot + 
  # Faded regression lines for each compound
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, 
              aes(group = compound, color = DOM_Type), 
              size = 0.4, linetype = "dashed", alpha = 0.2) +  # Dashed & transparent
  
  # Bold regression lines for each DOM_Type
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, 
              aes(group = DOM_Type, color = DOM_Type),
              size = 1.8)  # Thicker line for main trend

#6.5x5.5 preview

# stats for each compound
stats <- inv_resp_data %>%
  group_by(compound) %>%
  summarize(
    lm_model = list(lm(resp_rate ~ average_O, data = cur_data())),  # Fit linear model
    pearson_cor = cor(average_O, resp_rate, method = "pearson")  # Compute Pearson correlation
  ) %>%
  mutate(
    r_squared = sapply(lm_model, function(m) summary(m)$r.squared),  # Extract R²
    label = paste0("R² = ", round(r_squared, 3), "\n", "r = ", round(pearson_cor, 3))  # Format text
  )
stats

drop.cols <- c('lm_model')
stat = as.data.frame(stats %>% select(-one_of(drop.cols)))
stat
write.csv(stat, "lineweaver_stats.csv")


# calculating Kms
compounds = c("glucose1", "glucose3", "glucose4","fumaricacid1","fumaricacid2", "tryptophan2", "pyruvate1", "toluene1","toluene2","tyrosine1","tyrosine2", "benzoate1", "benzoate2","protochatechuate1", "protochatechuate2", "gentisate1","gentisate2")

compound_km = data.frame()

for (comp in compounds){
  data <- subset(inv_resp_data, compound == comp)
  model <- summary(lm(resp_rate~average_O, data = data))
  vmax = 1/(model$coefficients["(Intercept)","Estimate"])
  km = model$coefficients["average_O","Estimate"]*(1/(model$coefficients["(Intercept)","Estimate"]))
  temp = tibble(comp,km,vmax)
  compound_km <- rbind(compound_km, temp)
}

compound_km

write.csv(compound_km, "compound_kms.csv")


#### end ####

#### michealas menten curves ####
#used rates calculted from autobod in excel - file below
read.csv(resp_rate_compounds.csv)

# in file resp_rate_compounds_Mich-men.xlsx
#### end ####

#### boxplot comparing Km's ####
setwd("~/Desktop/initial_resp_R/")

km_data = read.csv("kms_boxplot.csv")
km_data

km_box = ggplot(km_data, aes(x=DOM.type, y=Km, fill=DOM.type)) +
  scale_fill_manual(values=c("skyblue3","firebrick3")) +
  geom_boxplot(width=0.4)+
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, shape = 21, color = "black", fill = "black") + 
  theme_minimal() + 
  scale_y_log10() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "",           
    y = "Km (μM)",     
    fill = "DOM Type"       
  ) +
  stat_compare_means(vjust = 1.2,label = "p.format")
km_box


#### end ####

#### Nitrogen cycling figure ####

# Load Dark2 color palette
library(RColorBrewer)
install.packages("forecast")
install.packages("colortools")
library(ggplot2)
library(forecast)
library(dplyr)
library(colortools)
library(patchwork)


# Filter the data to include only 'Hypoxic' and 'Oxic' treatments
setwd("~/Desktop/NOdata/")

nutrient_metadata_avg <- read.csv("nutrient_metadata_avg.csv")
nutrient_metadata_avg
filtered_data <- nutrient_metadata_avg %>%
  filter(treatment %in% c('Hypoxic', 'Oxic'))

# Create the plot
n_cycle <- ggplot(filtered_data, aes(x = day)) +
  geom_line(aes(y = ammonia, color = "Ammonia"), size = 1) +
  geom_point(aes(y = ammonia, color = "Ammonia"), size = 2) +
  geom_errorbar(aes(ymin = ammonia - sd_ammonia, ymax = ammonia + sd_ammonia), width = 0.2) +
  
  geom_line(aes(y = nitrite, color = "Nitrite"), size = 1) +
  geom_point(aes(y = nitrite, color = "Nitrite"), size = 2) +
  geom_errorbar(aes(ymin = nitrite - sd_nitrite, ymax = nitrite + sd_nitrite), width = 0.2) +
  
  geom_line(aes(y = nitrate, color = "Nitrate"), size = 1) +
  geom_point(aes(y = nitrate, color = "Nitrate"), size = 2) +
  geom_errorbar(aes(ymin = nitrate - sd_nitrate, ymax = nitrate + sd_nitrate), width = 0.2) +
  
  labs(
    x = NULL,
    y = "Concentration (µM/L)",
    color = "Variable"
  ) +
  scale_color_manual(
    values = brewer.pal(8, "Dark2"),  # Use the Dark2 palette
    # Specify colors
    name = "Analyte",
    breaks = c("Ammonia", "Nitrite", "Nitrate"),  # Specify the order
    labels = c("Ammonia                ", "Nitrite", "Nitrate")  # Specify the labels
  ) +
  facet_wrap(~ treatment) +
  theme_linedraw() +
  theme(
    legend.box = "horizontal",
    legend.box.background = element_rect(color = "white"),
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey", color = "black"),  # Customize facet labels
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center, make bold and larger title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

n_cycle

# oxygen utilization rates

oxy_data = read.csv("oxygen_ut_data.csv")
oxy_data
oxygen_ut <- ggplot(oxy_data, aes(x = day, group = treatment)) +
  geom_line(aes(y = Oxygen_ut_ave, color = treatment), size = 1) +
  geom_point(aes(y = Oxygen_ut_ave, color = treatment), size = 2) +
  geom_errorbar(aes(ymin = Oxygen_ut_ave - sd_ut, ymax = Oxygen_ut_ave + sd_ut, color = treatment), width = 0.2) +
  labs(
    x = "Time (Day)",
    y = "Rate (µM Oxygen/min)",
    color = "Variable"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("Oxic" = "skyblue3", "Hypoxic" = "firebrick3")  # You can choose your own hex colors
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),  # Customize facet labels
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center, make bold and larger title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
  
oxygen_ut

# statistical comparison of oxygen utilization boxplot

# remove outliers - timepoint 0-10, and timepoint 59.5
head(oxy_data2)
oxy_data2 = subset(oxy_data, oxy_data$day > 10)
oxy_data2 = oxy_data2 %>%
  filter(oxy_data2$day != 59.5)
oxy_data2


install.packages('ggpubr')
library(ggpubr)

plot2 <- ggplot(oxy_data2, aes(x = treatment, y = Oxygen_ut_ave, fill = treatment)) +
  scale_fill_manual(values = c("firebrick3", "skyblue3")) +
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, size = 2, shape = 21 ,fill="black", color = "black") + 
  ylab("Rate (µM Oxygen/min)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box = "horizontal",
    legend.box.background = element_rect(color = "white"),
    legend.position = "right"
  ) +
  stat_compare_means(method="t.test", label = "p.format", label.x = 1.5,  # x-position near the right
                     label.y = max(oxy_data2$Oxygen_ut_ave, na.rm = TRUE)*0.9)
  
plot2

#combine with oxygen rates
install.packages("cowplot")
library(cowplot)

oxyplot = plot_grid(oxygen_ut, plot2, 
          ncol = 2, 
          align = "v", 
          rel_widths = c(3, 1.1))

oxyplot


# sig taxa at day 96 in hypoxic for nitrogen cycle 
setwd("~/Desktop/NOdata/")

st <- read.csv("sig_taxa_96.csv")
st

# Create the plot
st_96 <- ggplot(st, aes(x = day)) +
  geom_line(aes(y = nitro_pa, color = "Nitrospinaceae"), size = 1) +
  geom_point(aes(y = nitro_pa, color = "Nitrospinaceae"), size = 2) +
  geom_errorbar(aes(ymin = nitro_pa - nitro_pa_sd, ymax = nitro_pa + nitro_pa_sd), width = 0.2) +
  
  geom_line(aes(y = methyl_pa, color = "Methyloligellaceae"), size = 1) +
  geom_point(aes(y = methyl_pa, color = "Methyloligellaceae"), size = 2) +
  geom_errorbar(aes(ymin = methyl_pa - methyl_pa_sd, ymax = methyl_pa + methyl_pa_sd), width = 0.2) +
  
  geom_line(aes(y = favo_pa, color = "Flavobacteriales NS9 Marine Group"), size = 1) +
  geom_point(aes(y = favo_pa, color = "Flavobacteriales NS9 Marine Group"), size = 2) +
  geom_errorbar(aes(ymin = favo_pa - favo_pa_sd, ymax = favo_pa + favo_pa_sd), width = 0.2) +
  

  labs(
    x = NULL,
    y = "Relative Percent Abundance",
    color = "Variable"
  ) +
  scale_color_manual(
    values = brewer.pal(8, "Accent"),
    # Specify colors
    name = "Taxa (Family)",
    breaks = c("Nitrospinaceae", "Methyloligellaceae", "Flavobacteriales NS9 Marine Group"),  # Specify the order
    labels = c("Nitrospinaceae", "Methyloligellaceae", "Flavobacteriales NS9 Marine Group")  # Specify the labels
  ) +
  facet_wrap(~ treatment) +
  theme_linedraw() +
  theme(
    legend.box = "horizontal",
    legend.box.background = element_rect(color = "white"),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey", color = "black"),  # Customize facet labels
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center, make bold and larger title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

st_96

#10x3




















# nitrogen cycling taxa abundance

n_cyclers <- read.csv("n_cyclers.csv")

n_cycler_plot <- ggplot(n_cyclers, aes(x = day)) +
  geom_line(aes(y = Nitrincolaceae, color = "Nitrincolaceae"), size = 1) +
  geom_point(aes(y = Nitrincolaceae, color = "Nitrincolaceae"), size = 2) +
  
  geom_line(aes(y = Nitrosopumilaceae, color = "Nitrosopumilaceae"), size = 1) +
  geom_point(aes(y = Nitrosopumilaceae, color = "Nitrosopumilaceae"), size = 2) +
  
  geom_line(aes(y = Nitrospinaceae, color = "Nitrospinaceae"), size = 1) +
  geom_point(aes(y = Nitrospinaceae, color = "Nitrospinaceae"), size = 2) +
  
  geom_line(aes(y = Nitrosomonadaceae, color = "Nitrosomonadaceae"), size = 1) +
  geom_point(aes(y = Nitrosomonadaceae, color = "Nitrosomonadaceae"), size = 2) +
  
  labs(
    x = "Time (Day)",
    y = "Average percent abundance",
    color = "Variable"
  ) +
  scale_color_manual(
    values = brewer.pal(10, "Paired"),  # Use the Dark2 palette
    # Specify colors
    name = "Taxa",
    breaks = c("Nitrincolaceae", "Nitrosopumilaceae", "SAR86", "SAR11_CladeI", "SAR11_CladeII", "Saprospiraceae", "Actinomarinaceae", "Nitrospinaceae", "Cryomorphaceae", "Nitrosomonadaceae"),  # Specify the order
    labels = c("Nitrincolaceae", "Nitrosopumilaceae", "SAR86", "SAR11_CladeI", "SAR11_CladeII", "Saprospiraceae", "Actinomarinaceae", "Nitrospinaceae", "Cryomorphaceae", "Nitrosomonadaceae")  # Specify the labels
  ) +
  facet_wrap(~ treatment) +
  theme_linedraw() +
  theme(
    legend.box = "horizontal",
    legend.box.background = element_rect(color = "white"),
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "white", color = "black"), # Customize facet labels
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center, make bold and larger title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

n_cycler_plot


#final plot
install.packages('Rmisc')
library(Rmisc)
multiplot(n_cycle, n_cycler_plot, oxyplot, cols = 1)

#### end ####

#### oxygen utilization model fitting ####
setwd("~/Desktop/final_resp_R/")
df_all_mods = read.csv("best_fit_models.csv")

# plot em all
# Plot with legend showing only "Model 1" and "Model 2" (no equations)
plot_mod <- ggplot(df_all_mods, aes(x = x, y = y, color = type, linetype = group, shape = compound)) +
  geom_line(size = 1) +
  scale_color_manual(values = color_vals, labels = c("Exponential", "Sigmoidal")) +
  labs(
    title = "Best fit models",
    x = "Time (days)", y = "[Oxygen] (uM/L)",
    color = "Model", linetype = "Compound"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2)),
    linetype = guide_legend(override.aes = list(size = 1))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(size = 9)
  )
plot_mod

# or 
plot_mod <- ggplot(df_all_mods, aes(x = x, y = y, color = group, linetype = type,shape=compound)) +
  geom_line(size = 1) +
  scale_color_manual(values=c("skyblue3","firebrick3"), labels = c("LDOM", "ODDOM")) +
  labs(
    title = "Best fit models",
    x = "Time (days)", y = "[Oxygen] (uM/L)",
    color = "Model", linetype = "Compound"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2)),
    linetype = guide_legend(override.aes = list(size = 1))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(size = 9)
  )
plot_mod

# 7x5?

#### end ####
