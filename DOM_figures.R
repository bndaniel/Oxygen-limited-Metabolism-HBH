# R code for DOM related figures in Persistence of Oxygen Dependent Dissolved Organic Matter (ODDOM) induced by Oxygen-limited Metabolism in Hypoxic Marine Environments"

#### Overall weighted DOM characteristics between treatments ####
#### characteristics through time plots ####
library(vegan)
library(ggplot2)
head(data)
#compiled averages of each timepoint
setwd("~/Desktop/HBH_DOM/DOM_R2")

#for figure 1
### only doing AI, avCox, DBE, and TOC and adding MLB and iSD
#data for through time ai, avcox, dbe, and TOC
library(dplyr)
data=read.csv("weighted_DOM_all_r_treat.csv",header=T)
data = subset(data, data$type != "RT")
data_sub = subset(data, type %in% c("AI", "avCox","DBE","TOC"))
data_sub <- data_sub %>% select(-group)
head(data_sub)

#data for boxplot ai, avcox, dbe, and TOC
data=read.csv("weighted_DOM_ave_r.csv",header=T)
data_sub0 = subset(data, treat %in% c("AI", "avCox","DBE","TOC"))
data_sub0 <- na.omit(data_sub0)
data2 = subset(data_sub0, data$Time!= 0)
data2 <- na.omit(data2)
data2 <- data2 %>% select(-group)
data2 <- data2 %>% rename(type = treat)
head(data2)

#data for through time mlb and isd
all = read.csv("DOM_lability_chemodiv.csv")
all_sub = subset(all, type %in% c("MLB", "iSD"))
all_sub <- all_sub %>% select(-X)
head(all_sub)

# data for boxplot mlb and isd
dv = read.csv("DOM_lability_chemodiv_ave.csv")
dv = subset(dv, dv$Time != 0)
all_sub_0 = subset(dv, type %in% c("MLB", "iSD"))
all_sub_0 <- all_sub_0 %>% select(-X)
head(all_sub_0)

# combine datasets
#through time
tt = rbind(data_sub, all_sub)
bp = rbind(data2,all_sub_0)
tt
bp
#through time plots
plot_fig1_a = ggplot(tt, aes(x = Time, y = value, colour = treatment, fill = treatment)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_colour_manual(values = c("firebrick3", "skyblue3")) +
  scale_fill_manual(values = c("firebrick3", "skyblue3")) +
  scale_x_continuous(breaks = tt$Time, labels = tt$Time) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = 'top',
    strip.text = element_blank()
  ) +
  facet_grid(vars(factor(type, levels = c("AI", "avCox", "DBE", "TOC","MLB", "iSD"))), scales = 'free')

plot_fig1_a

#boxplots
plot_fig1_b = ggplot(bp, aes(x=treatment, y=value, fill=treatment)) +
  scale_fill_manual(values=c("firebrick3", "skyblue3")) +
  geom_boxplot(width=0.4)+
  theme_bw()+theme(panel.grid = element_blank(),legend.position='top',axis.text.y = element_blank(),   # Remove y-axis text
                   axis.ticks.y = element_blank(),  # Remove y-axis ticks
                   axis.title.y = element_blank())+ # Remove y-axis title)+
  stat_compare_means(data, comparisons = list(c("Hypoxic", "Saturated")), 
                     label = "p.value", method = "t.test", paired=T, vjust = 1.2)+
  facet_grid(vars(factor(type, levels=c("AI", "avCox", "DBE", "TOC","MLB", "iSD"))), scales = 'free')

plot_fig1_b


#combine em
library(cowplot)

plot_grid(plot_fig1_a,plot_fig1_b, 
          ncol = 2, 
          align = "v", 
          rel_widths = c(2, 0.8))


# permanova statistical test

# load data
data=read.csv("weighted_DOM_ave_r.csv",header=T)
data2 = subset(data, data$Time!= 0)
data2 <- data2 %>% rename(type = treat)
head(data2)

# data for mlb and isd
dv = read.csv("DOM_lability_chemodiv_ave.csv")
dv = subset(dv, dv$Time != 0)
all_sub_0 = subset(dv, type %in% c("MLB", "iSD"))
all_sub_0 <- all_sub_0 %>% select(-X)
head(all_sub_0)

#combine them
per = rbind(data2,all_sub_0)

library(vegan)
per
# Ensure variables are correctly formatted
per$treatment <- factor(per$treatment)
per$type <- factor(per$type)

library(tidyr)
library(dplyr)

bp_wide <- per %>%
  pivot_wider(names_from = type, values_from = value)
bp_wide
treatment_vector <- bp_wide$treatment
data_matrix <- bp_wide %>%
  dplyr::select(where(is.numeric),-Time)
data_matrix
set.seed(123)  # for reproducibility
permanova_result <- adonis2(data_matrix ~ treatment_vector, method = "euclidean", permutations = 999)

print(permanova_result)

### 4/29/25 can we add a lability comparison of overal DOM between each sample through time? 
# https://www.sciencedirect.com/science/article/pii/S0146638023001134
# % MLB = 100 x (# of molecules with H/C >= 1.5 / total # of molecules)


#### PCA plot ####

library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
data=read.csv("weighted_DOM_ave_r_PCA.csv",header=T,row.names=1)
data<-data[,-c(1,2)]
data

#remove RT
data <- subset(data, select = -RT)
data <- subset(data, select = -mz)


group=read.csv("weighted_DOM_ave_r_PCA.csv",header=T,row.names=1)
#group<-data[,c(1,2)]
data.pca<-rda(data,scale=TRUE)
summary(data.pca)
ggplot()+
  geom_point(aes(x=data.pca$CA$u[,1],y=data.pca$CA$u[,2],color=group$site, shape=group$season),size=6)+
  geom_hline(yintercept=0,linetype=3,size=1)+geom_vline(xintercept=0,linetype=3,size=1)+
  geom_segment(aes(x=0,y=0,xend=data.pca$CA$v[,1],yend=data.pca$CA$v[,2]),arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),type="closed"),linetype=1,size=1.5,color="orange")+
  geom_text_repel(aes(x=data.pca$CA$v[,1],y=data.pca$CA$v[,2]),label=row.names(data.pca$CA$v),color="orange")+
  theme(legend.title=element_blank())+labs(x="PCA1(71.5%)",y="PCA2(22.8%)",shape="Timepoint",color="Treatment")+theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values=c("firebrick3", "skyblue3"))+
  scale_shape_manual(values=c(16,17,15,18,8))


#### PCA with everything #### not sure this is the best to visualize ####

data=read.csv("weighted_DOM_all_r_treat_PCA.csv",header=T,row.names=1)
data
#remove RT and TOC
data <- subset(data, select = -RT)
#data <- subset(data, select = -TOC)
data <- subset(data, select = -OC)

group = read.csv("weighted_DOM_all_r_treat_PCA.csv",header=T,row.names=1)
data<-data[,-c(1,2,3)]

data.pca<-rda(data,scale=TRUE)
summary(data.pca)

ggplot()+
  geom_point(aes(x=data.pca$CA$u[,1],y=data.pca$CA$u[,2],color=group$site, shape=group$season),size=6)+
  geom_hline(yintercept=0,linetype=3,size=1)+geom_vline(xintercept=0,linetype=3,size=1)+
  geom_segment(aes(x=0,y=0,xend=data.pca$CA$v[,1],yend=data.pca$CA$v[,2]),arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),type="closed"),linetype=1,size=1.5,color="orange")+
  geom_text_repel(aes(x=data.pca$CA$v[,1],y=data.pca$CA$v[,2]),label=row.names(data.pca$CA$v),color="orange")+
  theme(legend.title=element_blank())+labs(x="PCA1(63.6%)",y="PCA2(30.6%)",shape="Timepoint",color="Treatment")+theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values=c("firebrick3", "skyblue3"))+
  scale_shape_manual(values=c(16,17,15,18,8))
#### end ####


#### DOC through time plot ####
setwd("~/Desktop/HBH_DOM")

data=read.csv("TOC.csv",header=T)

ggplot(data,aes(x= data$Time, y = value, colour = treatment)) + geom_point() + geom_smooth(method="loess")+ scale_colour_manual(values=c("firebrick3", "skyblue3")) +
  scale_x_continuous(breaks = data$Time, labels = data$Time, xlab("Time")) +
  theme_bw()+theme(panel.grid = element_blank())+
  facet_grid(vars(factor(treat, levels=c('TOC'))), scales = 'free')

#### end #### 

#### sorted DOM plots ####
#### violin ####
df_vk <- read.csv('elemental_info_sorted_dom_ave.csv',check.names = T,header=T)

# add in NOSC and DBE/c
df_vk$NOSC = 4 - (((4*df_vk$C_Sirius) + df_vk$H_Sirius - (2*df_vk$O_Sirius) - (3*df_vk$N_Sirius) + (5*df_vk$P_Sirius) - (2*df_vk$S_Sirius))/df_vk$C_Sirius)
df_vk$DBE_C = df_vk$DBE/df_vk$C_Sirius

data = df_vk
head(data)

#remove OIDOM
df_vk1 = subset(df_vk, df_vk$DOM!="OIDOM")
df_vk1      

#order DOM into LDOM, RDOM, ODDOM so that they overlap nice
li = list("LDOM","ODDOM","RDOM")
df_vk1$DOM = factor(df_vk1$DOM, levels=c("LDOM","ODDOM","RDOM"))
df_vk1= df_vk1 %>%
  arrange(DOM) 

#groups to compare

library(ggplot2)
library(ggpubr)

my_comparisons <- list(c("ODDOM", "LDOM"), c("ODDOM", "RDOM"), c("LDOM","RDOM"))

data = df_vk1
head(data)

compare_means(OC~DOM, data=data)

p1<-ggplot(data, aes(x=DOM, y=OC, fill =DOM ,color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "O/C", scales = "free_y")

p1

#Comparison of HC
compare_means(HC~DOM, data=data)

p2<-ggplot(data, aes(x=DOM, y=HC, fill =DOM ,color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "H/C", scales = "free_y")


#comparison of average C oxidation
compare_means(avCox~DOM, data=data)

p3<-ggplot(data, aes(x=DOM, y=NOSC, fill =DOM ,color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "NOSC", scales = "free_y")


#comparison of DBE
compare_means(DBE~DOM, data=data)

p4<-ggplot(data, aes(x=DOM, y=DBE, fill =DOM ,color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "DBE", scales = "free_y")


#comparison of AI
compare_means(AI~DOM, data=data)

p5<-ggplot(data, aes(x=DOM, y=AI, fill =DOM ,color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "AImod", scales = "free_y")
p5

#comparison of mz
compare_means(DBE_C~DOM, data=data)

p7<-ggplot(data, aes(x=DOM, y=DBE_C, fill =DOM, color= DOM)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ "DBE/C", scales = "free_y")
  
p7

#can add this to shift astriks into view + scale_y_continuous(limit=c(0,1250))
#combine them all better

#just have to add legend by changing legend.position = 'right' for one and edit in illustrator
plot_grid(p5+ theme(legend.position = "none"),p4+ theme(legend.position = "none"),p7+ theme(legend.position = "none"),p2+ theme(legend.position = "none"),p3+ theme(legend.position = "none"),p1+ theme(legend.position = "none"), ncol = 3, align = "v")

#7x5 landscape

#### van krev ####

li = list("LDOM","RDOM","ODDOM")
df_vk1$DOM = factor(df_vk1$DOM, levels=c("RDOM","LDOM","ODDOM"))
df_vk1= df_vk1 %>%
  arrange(DOM) 

#vk of HC vs OC
VK3 <- ggplot(df_vk1, aes(x = OC, y = HC, fill=`DOM`,color = `DOM`)) +
  geom_point(shape = 21, size = 3.5, stroke = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),legend.position = c(0.85, 0.12)) +
  scale_fill_manual(values = c("#C2A5CF","#ADEDED","#FA9A98"),
                    labels = c( "RDOM (n=950)","LDOM (n=392)","ODDOM (n=214)")) +
  scale_color_manual(values = c("#8073AC","#2166AC","#FE3298"),
                     labels = c("RDOM (n=950)","LDOM (n=392)","ODDOM (n=214)")) +
  stat_ellipse() + xlim(0, 1) +
  labs(x = "O/C", y = "H/C") +
  scale_size(range=c(3,10))+scale_x_continuous(limits = c(0, 0.9), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) + stat_ellipse() +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))

VK3  


#### stack plot of superclass ####
setwd("~/Desktop/HBH_DOM/DOM_R2")
library(RColorBrewer)

#upload each sorted DOM superclass percents (see methods for how it was done)
H_LDOM = read.csv("percents_LDOM_H.csv")
H_ODDOM = read.csv("percents_ODDOM_H.csv")
H_RDOM = read.csv("percents_RDOM_H.csv")

# convert group into day value
H_LDOM <- H_LDOM %>%
  mutate(group = case_when(
    group == "Avg0_H" ~ "0", group == "Avg12_H" ~ "12",group == "Avg42_H" ~ "42" ,group == "Avg70_H" ~ "70", group == "Avg96_H" ~ "96"))  # Assign numeric value to each group
H_ODDOM <- H_ODDOM %>%
  mutate(group = case_when(
    group == "Avg0_H" ~ "0", group == "Avg12_H" ~ "12",group == "Avg42_H" ~ "42" ,group == "Avg70_H" ~ "70", group == "Avg96_H" ~ "96"))  # Assign numeric value to each group
H_RDOM <- H_RDOM %>%
  mutate(group = case_when(
    group == "Avg0_H" ~ "0", group == "Avg12_H" ~ "12",group == "Avg42_H" ~ "42" ,group == "Avg70_H" ~ "70", group == "Avg96_H" ~ "96"))  # Assign numeric value to each group

#colors
H_ODDOM$Group_name <- as.factor(H_ODDOM$Group_name)
H_LDOM$Group_name <- as.factor(H_LDOM$Group_name)
H_RDOM$Group_name <- as.factor(H_RDOM$Group_name)
color_levels <- levels(H_ODDOM$Group_name)
nb.cols <- length(color_levels)

#attach colors to group names
install.packages("pals")
library(pals)
base_colors <- tol(12)
myColors <- colorRampPalette(base_colors)(nb.cols)
names(myColors) <- color_levels
colScale <- scale_color_manual(name = "Group_name", values = myColors)


#making the stack plot!
#LDOM
stackPlot <- ggplot(H_LDOM, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme_minimal() +  
  ggtitle("LDOM") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = myColors) 
stackPlot

#ODDOM
stackPlot2 <- ggplot(H_ODDOM, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme_minimal() + 
  labs(y=NULL) +
  ggtitle("ODDOM") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(),  # Removes y-axis label
        axis.text.y = element_blank(),   # Removes y-axis text
        axis.ticks.y = element_blank())+
  scale_fill_manual(values = myColors) 
stackPlot2

#RDOM
stackPlot3 <- ggplot(H_RDOM, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme_minimal() +  
  ggtitle("RDOM") +
  labs(y=NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(),  # Removes y-axis label
        axis.text.y = element_blank(),   # Removes y-axis text
        axis.ticks.y = element_blank())+
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = myColors) 
stackPlot3

#plot it up
install.packages("patchwork")  
library(patchwork)

stackPlot+stackPlot2+stackPlot3+ theme(legend.position = "none")


#### end ####

#### Kegg identified ODDOM comparison to experimental ODDOM ####
#### superclass identification annotations ####
setwd("~/Desktop/HBH_DOM/DOM_R2")
library(RColorBrewer)

#upload each sorted DOM superclass percents (see methods for how it was done)
K_sup = read.csv("kegg_superclass.csv")
K_sup
#colors
K_sup$Superclass <- factor(
  K_sup$Superclass,
  levels = unique(K_sup$Superclass))
color_levels <- levels(K_sup$Superclass)

nb.cols <- length(color_levels)

#attach colors to group names
install.packages("pals")
library(pals)
myColors <- c("#332288FF", "#6699CCFF", "#88CCEEFF", "#44AA99FF", "#117733FF", "#999933FF", "#DDCC77FF", "#CC6677FF", "#AA4466FF", "#882255FF", "#AA4499FF","#B7957CFF", "#A68A7BFF", "#734939FF", "#A6432DFF", "#9BA8AEFF")
names(myColors) <- color_levels
colScale <- scale_color_manual(name = "Superclass", values = myColors)


#making the stack plot!
#kegg ODDOM
stackPlot <- ggplot(K_sup, aes(fill=Superclass, y=Percent, x=Group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic") +
  xlab(NULL) + ylab("Percentage") + labs(fill = "Superclass") + 
  theme_minimal() +  
  ggtitle("KEGG ODDOM") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
  ) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = myColors) 
stackPlot


#### end ####


#### violin ####
# chemical characteristics of Kegg oddom, experimental oddom, and hbh oddom
# read in all formulas with metadata on group - kegg, hbh, or exp

Feat = read.csv('kegg_exp_oddom_formulas.csv',check.names = T,header=T)
Feat

# First, Sirius molecular formulas

#splits the molecular formula if alphabets are present in it such as 'C13' 'H22' 'S'
if (!require("CHNOSZ")) install.packages("CHNOSZ")
library(CHNOSZ)
library(data.table)

formulas <- makeup(Feat$Formula, count.zero = TRUE)
formulas <- data.table(compound = Feat$Formula)[
  ,names(formulas[[1]]) := transpose(formulas)
]

head(formulas)
formulas$group = Feat$Group
formulas

#change order to match ?

sf <- as.data.frame(formulas)
rownames(sf) <- NULL
colnames(sf) <- paste(colnames(sf))

head(sf)

df <- sf
head(df)
#C:N, O:C, H:C and average C oxidation columns using only Sirius formulas
df$CN = df$C/df$N
df$OC = df$O/df$C
df$HC = df$H/df$C
df$avCox <- -((1*df$H) - (3*df$N) - (2*df$O) + (5*df$P) -(2*df$S))/(df$C)
df$DBE <- (1 + 0.5*(2*df$C - df$H + df$N + df$P))
df$AI <- (1 + df$C - (0.5*df$O) - df$S - 0.5*(df$N + df$P + df$H))/(df$C - 0.5*(df$O) - df$N - df$S - df$P)
df$NOSC = 4 - (((4*df$C) + df$H - (2*df$O) - (3*df$N) + (5*df$P) - (2*df$S))/df$C)
df$DBE_C = df$DBE/df$C

head(df)

#making vk plot 
#Only selecting formulas that make sense molecularly (chemically)
df_vk <- df %>% filter(avCox > -4 & avCox < 4 & OC < 2 & HC < 2.5 & CN > 4)
write.csv(df_vk,"kegg_exp_oddom_elemental.csv")
head(df_vk)
range(df_vk$avCox)

#need to filter AI (any values less than 0 = 0)
df_vk <-read.csv("kegg_exp_oddom_elemental.csv",check.names = T,header=T)
df_vk$AI[df_vk$AI<0] = 0

data = df_vk
tail(data)

li = list("KEGG ODDOM","EXP ODDOM","EXP LDOM","EXP RDOM")
data$group = factor(data$group, levels=c("KEGG ODDOM","EXP ODDOM","EXP LDOM","EXP RDOM"))
data= data %>%
  arrange(group) 

#groups to compare

library(ggplot2)
library(ggpubr)

#pre subset, save all groups in data1
data1 = data

#subset data to exculde ldom and rdom
data <- subset(data, !(group %in% c("EXP LDOM", "EXP RDOM")))
tail(data)
my_comparisons <- list(c("KEGG ODDOM", "EXP ODDOM"))

compare_means(OC~group, data=data)

p1<-ggplot(data, aes(x=group, y=OC, fill =group ,color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98")) +
  scale_color_manual(values = c("#007200", "#FE3298")) +
  facet_wrap(~ "O/C", scales = "free_y")

p1

#Comparison of HC
compare_means(HC~group, data=data)

p2<-ggplot(data, aes(x=group, y=HC, fill =group ,color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98")) +
  scale_color_manual(values = c("#007200", "#FE3298")) +
  facet_wrap(~ "H/C", scales = "free_y")


#comparison of average C oxidation
compare_means(avCox~group, data=data)

p3<-ggplot(data, aes(x=group, y=NOSC, fill =group ,color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98")) +
  scale_color_manual(values = c("#007200", "#FE3298")) +
  facet_wrap(~ "NOSC", scales = "free_y")


#comparison of DBE
compare_means(DBE~group, data=data)

p4<-ggplot(data, aes(x=group, y=DBE, fill =group ,color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98")) +
  scale_color_manual(values = c("#007200", "#FE3298")) +
  facet_wrap(~ "DBE", scales = "free_y")


#comparison of dbe/c
compare_means(DBE_C~group, data=data)

p7<-ggplot(data, aes(x=group, y=DBE_C, fill =group, color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98")) +
  scale_color_manual(values = c("#007200", "#FE3298")) +
  facet_wrap(~ "DBE/C", scales = "free_y")

p7

#comparison of AImod
my_comparisons <- list(c("KEGG ODDOM", "EXP ODDOM"),c("KEGG ODDOM", "EXP LDOM"),c("KEGG ODDOM", "EXP RDOM"))

compare_means(AI~group, data=data1)

p5<-ggplot(data1, aes(x=group, y=AI, fill =group ,color= group)) +
  geom_violin(trim = FALSE, alpha = 0.7,linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill = NA,linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01,vjust = 1.5) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98", "#ADEDED", "#C2A5CF")) +
  scale_color_manual(values = c("#007200", "#FE3298", "#2166AC", "#8073AC")) +
  facet_wrap(~ "AImod", scales = "free_y")
p5


#can add this to shift astriks into view + scale_y_continuous(limit=c(0,1250))
#combine them all better

#just have to add legend by changing legend.position = 'right' for one and edit in illustrator
plot_grid(p5+ theme(legend.position = "none"),p4+ theme(legend.position = "none"),p7+ theme(legend.position = "none"),p2+ theme(legend.position = "none"),p3+ theme(legend.position = "none"),p1+ theme(legend.position = "none"), ncol = 3, align = "v")

#7x5 landscape

#### vank krev ####
VK3 <- ggplot(data, aes(x = OC, y = HC, fill=group, color = group)) +
  geom_point(shape = 21, size = 3.5, stroke = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),legend.position = c(0.85, 0.12)) +
  scale_fill_manual(values = c("#00bb00", "#FA9A98"),
                    labels = c("KEGG ODDOM (n=528)", "EXP ODDOM (n=214)")) +
  scale_color_manual(values = c("#007200", "#FE3298"),
                     labels = c("KEGG ODDOM (n=528)", "EXP ODDOM (n=214)")) +
  labs(x = "O/C", y = "H/C", size='Average Abundance') +
  scale_size(range=c(3,10))+scale_x_continuous(limits = c(0, 0.9), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) + stat_ellipse() +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
VK3


#7x5 landscape

#### end ####

#### PCA of rel ab DOM ####
#took relative abundance averagers between duplicates of each replicate
#see methods in DOM_featuretable_all_relative_analysis_SORT.xlsx

setwd("~/Desktop/HBH_DOM/DOM_R2")

ndoc = read.csv("PCA_rel_ab_DOM_3_25.csv")

#subset to have one identifier - sample_replicate
ndoc1 = subset(ndoc,select=-c(Sample_name, Day, Treatment, Replicate))
head(ndoc1)

#make numeric
ndoc2 <- mutate_all(ndoc1, function(x) as.numeric(as.character(x)))
head(ndoc2)


#bray-curtis distance 
ndist = vegdist(ndoc2, method="bray",na.rm = TRUE)
ndist

#PCOA
npcoa = cmdscale(ndist, k=2, eig=TRUE, add=TRUE)
npcoa

npositions = npcoa$points
npositions
colnames(npositions) = c("pcoa1", "pcoa2")

#put in sample identifiers and treatment identifiers
npositions = data.frame(npositions)
npositions$sample = ndoc$Sample_name
npositions$Treatment = ndoc$Treatment
npositions$Day = ndoc$Day
npositions

#put in identifiers as left most column
npositions1 = subset(npositions, select=c(sample, Treatment, Day ,pcoa1,pcoa2))
npositions1

npositions1$Day= factor(npositions1$Day, levels =c("0","2","14","42","70","96"))
#figure out percent explained by each axis
percent_explained = 100*npcoa$eig / sum(npcoa$eig)
percent_explained[1:2]

#PCA comparing treatments
PCA_DOM = npositions1 %>%
  as_tibble() %>%
  ggplot(aes(x=pcoa1, y=pcoa2, colour=Day, shape=Treatment, group = Day)) +
  geom_point(size=3) + 
  stat_ellipse(linetype=2, size=0.5) +
  theme_bw() +
  scale_color_manual(values=c("#1B9E77","skyblue2", "#D95F02", "#7570B3", "#E7298A", "#66A61E"))+
  theme_bw() +labs(x="PCo 1 (44.08%)", y="PCo 2 (21.67%)") +
  scale_shape_manual(values=c(16,17)) 
PCA_DOM

#### end ####

#### comparison of MLB, and DOM diversity supp ####

# lability index
#need to do presence vs absence for each compound at each timepoint within each treatment... what is the best threshold?
# 50% or greater to be present across duplicates and replicates at each timepoint?
library(dplyr)

#create a dataframe that filters each timepoint for presence or absence, and then calculates lability index
# based on this: https://www.sciencedirect.com/science/article/pii/S0146638023001134
# first, filter for DOM within range, 0 to 1.2 O/C ratio and above 1.5 H/C ratio
df = read.csv("ave_dup_DOM_all_r_mlb.csv")
df
df <- df %>% filter( OC <= 1.2 & OC >= 0)

#how many compounds?
tc = nrow(df)
tc
#how many compounds after going above 1.5?
dfh = df %>% filter(HC>=1.5)
hc = nrow(dfh)

mlb = 100*(hc/tc)
mlb

#subset specific to each timepoint and calculate %MLB 
#bloom
mlb_bloom1 = 100*(nrow(df %>% filter(bloom_1 > 0 & HC>=1.5)))/(nrow(subset(df, bloom_1 > 0)))
mlb_bloom2 = 100*(nrow(df %>% filter(bloom_2 > 0 & HC>=1.5)))/(nrow(subset(df, bloom_2 > 0)))
mlb_bloom3 = 100*(nrow(df %>% filter(bloom_3 > 0 & HC>=1.5)))/(nrow(subset(df, bloom_3 > 0)))

#hypoxic 
mlb_H0_1 = 100*(nrow(df %>% filter(X0_H_1 > 0 & HC>=1.5)))/(nrow(subset(df, X0_H_1 > 0)))
mlb_H0_2 = 100*(nrow(df %>% filter(X0_H_2 > 0 & HC>=1.5)))/(nrow(subset(df, X0_H_2 > 0)))
mlb_H0_3 = 100*(nrow(df %>% filter(X0_H_3 > 0 & HC>=1.5)))/(nrow(subset(df, X0_H_3 > 0)))

mlb_H12_1 = 100*(nrow(df %>% filter(X12_H_1 > 0 & HC>=1.5)))/(nrow(subset(df, X12_H_1 > 0)))
mlb_H12_2 = 100*(nrow(df %>% filter(X12_H_2 > 0 & HC>=1.5)))/(nrow(subset(df, X12_H_2 > 0)))
mlb_H12_3 = 100*(nrow(df %>% filter(X12_H_3 > 0 & HC>=1.5)))/(nrow(subset(df, X12_H_3 > 0)))

mlb_H42_1 = 100*(nrow(df %>% filter(X42_H_1 > 0 & HC>=1.5)))/(nrow(subset(df, X42_H_1 > 0)))
mlb_H42_2 = 100*(nrow(df %>% filter(X42_H_2 > 0 & HC>=1.5)))/(nrow(subset(df, X42_H_2 > 0)))
mlb_H42_3 = 100*(nrow(df %>% filter(X42_H_3 > 0 & HC>=1.5)))/(nrow(subset(df, X42_H_3 > 0)))

mlb_H70_1 = 100*(nrow(df %>% filter(X70_H_1 > 0 & HC>=1.5)))/(nrow(subset(df, X70_H_1 > 0)))
mlb_H70_2 = 100*(nrow(df %>% filter(X70_H_2 > 0 & HC>=1.5)))/(nrow(subset(df, X70_H_2 > 0)))
mlb_H70_3 = 100*(nrow(df %>% filter(X70_H_3 > 0 & HC>=1.5)))/(nrow(subset(df, X70_H_3 > 0)))

mlb_H96_1 = 100*(nrow(df %>% filter(X96_H_1 > 0 & HC>=1.5)))/(nrow(subset(df, X96_H_1 > 0)))
mlb_H96_2 = 100*(nrow(df %>% filter(X96_H_2 > 0 & HC>=1.5)))/(nrow(subset(df, X96_H_2 > 0)))
mlb_H96_3 = 100*(nrow(df %>% filter(X96_H_3 > 0 & HC>=1.5)))/(nrow(subset(df, X96_H_3 > 0)))

#saturated
mlb_S0_1 = 100*(nrow(df %>% filter(X0_S_1 > 0 & HC>=1.5)))/(nrow(subset(df, X0_S_1 > 0)))
mlb_S0_2 = 100*(nrow(df %>% filter(X0_S_2 > 0 & HC>=1.5)))/(nrow(subset(df, X0_S_2 > 0)))
mlb_S0_3 = 100*(nrow(df %>% filter(X0_S_3 > 0 & HC>=1.5)))/(nrow(subset(df, X0_S_3 > 0)))

mlb_S12_1 = 100*(nrow(df %>% filter(X12_S_1 > 0 & HC>=1.5)))/(nrow(subset(df, X12_S_1 > 0)))
mlb_S12_2 = 100*(nrow(df %>% filter(X12_S_2 > 0 & HC>=1.5)))/(nrow(subset(df, X12_S_2 > 0)))
mlb_S12_3 = 100*(nrow(df %>% filter(X12_S_3 > 0 & HC>=1.5)))/(nrow(subset(df, X12_S_3 > 0)))

mlb_S42_1 = 100*(nrow(df %>% filter(X42_S_1 > 0 & HC>=1.5)))/(nrow(subset(df, X42_S_1 > 0)))
mlb_S42_2 = 100*(nrow(df %>% filter(X42_S_2 > 0 & HC>=1.5)))/(nrow(subset(df, X42_S_2 > 0)))
mlb_S42_3 = 100*(nrow(df %>% filter(X42_S_3 > 0 & HC>=1.5)))/(nrow(subset(df, X42_S_3 > 0)))

mlb_S70_1 = 100*(nrow(df %>% filter(X70_S_1 > 0 & HC>=1.5)))/(nrow(subset(df, X70_S_1 > 0)))
mlb_S70_2 = 100*(nrow(df %>% filter(X70_S_2 > 0 & HC>=1.5)))/(nrow(subset(df, X70_S_2 > 0)))
mlb_S70_3 = 100*(nrow(df %>% filter(X70_S_3 > 0 & HC>=1.5)))/(nrow(subset(df, X70_S_3 > 0)))

mlb_S96_1 = 100*(nrow(df %>% filter(X96_S_1 > 0 & HC>=1.5)))/(nrow(subset(df, X96_S_1 > 0)))
mlb_S96_2 = 100*(nrow(df %>% filter(X96_S_2 > 0 & HC>=1.5)))/(nrow(subset(df, X96_S_2 > 0)))
mlb_S96_3 = 100*(nrow(df %>% filter(X96_S_3 > 0 & HC>=1.5)))/(nrow(subset(df, X96_S_3 > 0)))

#create dataframe with MLB
Time = c(0,0,0,12,12,12,42,42,42,70,70,70,96,96,96,0,0,0,12,12,12,42,42,42,70,70,70,96,96,96)
treatment = c('Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Hypoxic','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated','Saturated')
value = c(mlb_H0_1,mlb_H0_2,mlb_H0_3,
          mlb_H12_1,mlb_H12_2,mlb_H12_3,
          mlb_H42_1, mlb_H42_2, mlb_H42_3,
          mlb_H70_1,mlb_H70_2,mlb_H70_3,
          mlb_H96_1,mlb_H96_2,mlb_H96_3,
          mlb_S0_1,mlb_S0_2,mlb_S0_3,
          mlb_S12_1, mlb_S12_2, mlb_S12_3,
          mlb_S42_1,mlb_S42_2,mlb_S42_3,
          mlb_S70_1,mlb_S70_2,mlb_S70_3,
          mlb_S96_1,mlb_S96_2,mlb_S96_3)
sample = c('0_H_1','0_H_2','0_H_3','12_H_1','12_H_2','12_H_3','42_H_1','42_H_2','42_H_3','70_H_1','70_H_2','70_H_3','96_H_1','96_H_2','96_H_3','0_S_1','0_S_2','0_S_3','12_S_1','12_S_2','12_S_3','42_S_1','42_S_2','42_S_3','70_S_1','70_S_2','70_S_3','96_S_1','96_S_2','96_S_3')
type =  c('MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB','MLB')

mlb = data.frame(sample,value,type,treatment,Time)
mlb


# chemodiversity index - shannon diversity and simpson - but should use simpson I think
# created a file with rel abundance only shown if present in 3 or more samples per timepoint
# see present_absent_DOM_methods.xlsx for details
install.packages("vegan")
library(vegan)

df = read.csv("ave_dup_DOM_all_r_div.csv")

dom_div = read.csv("DOM_diversity_table.csv", header=T,row.names = 1)
head(dom_div)

# Calculate diversity indices
shannon_diversity         <- diversity(dom_div, index = "shannon")
simpson_diversity         <- diversity(dom_div, index = "simpson")
inverse_simpson_diversity <- diversity(dom_div, index = "inv")

# View table
print(diversity_table)

type= c('SH','SH','SH','SH','SH','SH','SH','SH','SH','SH')
sh = data.frame(sample,shannon_diversity,type,treatment,Time)
colnames(sh)[2] <- "value"
sh

type= c('SD','SD','SD','SD','SD','SD','SD','SD','SD','SD')
sd = data.frame(sample,simpson_diversity,type,treatment,Time)
colnames(sd)[2] <- "value"
sd

type= c('iSD','iSD','iSD','iSD','iSD','iSD','iSD','iSD','iSD','iSD')
isd = data.frame(sample,inverse_simpson_diversity,type,treatment,Time)
colnames(isd)[2] <- "value"
isd

cdiv= rbind(sh,sd,isd)
cdiv

# compile with lability index

all = rbind(mlb,cdiv)
all

write.csv(all,"DOM_lability_chemodiv.csv")

# make plot with four measurements of lability and diversity
type_labels <- c(
  "iSD" = "Inverse Simpson Index",
  "SD" = "Simpson Index",
  "SH" = "Shannon Diversity",
  "MLB" = "Molecular Lability Boundary")

all$type <- factor(all$type, levels = c("MLB", "SH", "iSD","SD"))


ggplot(all, aes(x = Time, y = value, color = treatment)) +
  geom_point() +
  geom_smooth(se = TRUE) +  # Add se=TRUE if you want confidence intervals
  facet_wrap(~type, nrow = 2, scales = "free", labeller = labeller(type = type_labels)) +
  scale_color_manual(
    values = c(
      "Saturated" = "skyblue3",
      "Hypoxic" = "firebrick3"
    )
  ) +
  theme_minimal()
#### end ####

#### FT-ICR-MS meta-analysis ####
# get methods from Qi

library(ggplot2)
library(ggpubr)

setwd("~/Desktop/HBH_DOM/DOM_R2")

data=read.table("HBH_MF4.txt",header=T,row.names=1)
ggplot(data, aes(x = O_C, y = H_C, fill = type1, color = type1)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.5, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c( "#ADEDED","#C2A5CF", "#FA9A98")) +
  scale_color_manual(values = c( "#2166AC","#8073AC", "#FE3298"))


data <- read.table("HBH_MF5.txt", header = TRUE, row.names = 1)
my_comparisons <- combn(unique(data$type1), 2, simplify = FALSE)
ggplot(data, aes(x = type1, y = value, fill = type1, color = type1)) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", size = 5, tip.length = 0.01) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c("#ADEDED", "#FA9A98", "#C2A5CF")) +
  scale_color_manual(values = c("#2166AC","#FE3298", "#8073AC")) +
  facet_wrap(~ type2, scales = "free_y")




# old violins: 
p1<-ggviolin(data, x="DOM", y="OC", fill ="DOM" ,ylab="Oxygen/Carbon", add = "boxplot", palette = c("#66FFFF", "#FE9898","#CAF3CD"), line.color="gray", line.size=0.1, 
#facet.by = "type1" + 
short.panel.labs = FALSE) +
  #facet_grid(vars(type1)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14), axis.text.y=element_text(size=14),axis.title.y = element_text(size = 14),legend.position = 'none')+stat_compare_means(comparisons=my_comparisons,aes(label=..p.signif..))+xlab(NULL)
p1
#### end ####

#### VOC bar plots ####
# can find methods for normalization and relative abundance calculations in "voc_meta_new_version_for_blank_removal.xlsx"
# also annotations methods in there

# load VOC data relative abundance...

setwd("~/Desktop/VOC/")

voc = read.csv("VOC_sig_3_25.csv")

head(voc)

data = subset(voc, select=c(M.1, Day, treatment, ave_percent_ab, std_dev, compound))
data = na.omit(data)
data

#bar plot with one box per compound (variety description)
# error bars as std dev
data$M.1 <- factor(data$M.1, levels = unique(data$M.1))  
facet_labels <- setNames(paste0(data$compound, " (day ", data$Day, ")"), data$M.1)


day = ggplot(data, aes(fill=treatment, y=ave_percent_ab, x=treatment)) + 
  geom_bar(position="dodge", stat="identity",width=0.7) +
  scale_fill_manual(values=c("firebrick3", "skyblue3")) +
  facet_wrap(~M.1, scales="free", labeller = labeller(M.1 = facet_labels)) +  
  theme_minimal() +
  theme(strip.text = element_text(size=12,face="bold")) +
  geom_errorbar(aes(x=treatment, ymin=ave_percent_ab-std_dev, ymax=ave_percent_ab+std_dev), width=0.2, colour="black", alpha=0.9, size=.4)
day


#### PCA of rel ab VOCs ####
library(ggplot2)

library(vegan)

library(cluster)

library(tidyr)

library(dplyr)
library(viridis)
library(hrbrthemes)
setwd("~/Desktop/VOC/")

nvoc = read.csv("PCA_rel_ab_VOC_3_25.csv")

head(nvoc)


#subset to have one identifier - sample_replicate
nvoc1 = subset(nvoc,select=-c(Sample_name, Day, Treatment, Replicate))
head(nvoc1)

#make numeric
nvoc2 <- mutate_all(nvoc1, function(x) as.numeric(as.character(x)))
head(nvoc2)


#bray-curtis distance 
ndist = vegdist(nvoc2, method="bray",na.rm = TRUE)
ndist

#PCOA
npcoa = cmdscale(ndist, k=2, eig=TRUE, add=TRUE)
npcoa

npositions = npcoa$points
npositions
colnames(npositions) = c("pcoa1", "pcoa2")

#put in sample identifiers and treatment identifiers
npositions = data.frame(npositions)
npositions$sample = nvoc$Sample_name
npositions$Treatment = nvoc$Treatment
npositions$Day = nvoc$Day
npositions

#put in identifiers as left most column
npositions1 = subset(npositions, select=c(sample, Treatment, Day ,pcoa1,pcoa2))
npositions1

npositions1$Day= factor(npositions1$Day, levels =c("0","2","14","42","96"))
#figure out percent explained by each axis
percent_explained = 100*npcoa$eig / sum(npcoa$eig)
percent_explained[1:2]

#PCA comparing treatments
PCA_VOC = npositions1 %>%
  as_tibble() %>%
  ggplot(aes(x=pcoa1, y=pcoa2, colour=Day, shape=Treatment, group = Day)) +
  geom_point(size=3) + 
  stat_ellipse(linetype=2, size=0.5) +
  theme_bw() +
  scale_color_manual(values=c("#1B9E77", "skyblue2","#D95F02", "#7570B3", "#66A61E"))+
  theme_bw() +labs(x="PCo 1 (44.52%)", y="PCo 2 (18.17%)") +
  scale_shape_manual(values=c(16,17)) 
PCA_VOC

#### end ####

#### PCA of rel ab microbes ####
# needs to have relative abundances of each ASV 
# see rel ab calculation for each ASV in meta_16s(AutoRecovered).xlsx

setwd("~/Desktop/HBH")

ndoc = read.csv("PCA_rel_ab_MC_3_25_fam.csv")
#subset to have one identifier 
ndoc1 = subset(ndoc,select=-c(Sample_name, Day, Treatment, Replicate))


#make numeric
ndoc2 <- mutate_all(ndoc1, function(x) as.numeric(as.character(x)))



#bray-curtis distance 
ndist = vegdist(ndoc2, method="bray",na.rm = TRUE)
ndist

#PCOA
npcoa = cmdscale(ndist, k=2, eig=TRUE, add=TRUE)
npcoa

npositions = npcoa$points
npositions
colnames(npositions) = c("pcoa1", "pcoa2")

#put in sample identifiers and treatment identifiers
npositions = data.frame(npositions)
npositions$sample = ndoc$Sample_name
npositions$Treatment = ndoc$Treatment
npositions$Day = ndoc$Day
npositions

#put in identifiers as left most column
npositions1 = subset(npositions, select=c(sample, Treatment, Day ,pcoa1,pcoa2))
npositions1

npositions1$Day= factor(npositions1$Day, levels =c("0","2","14","42","70","96"))
#figure out percent explained by each axis
percent_explained = 100*npcoa$eig / sum(npcoa$eig)
percent_explained[1:2]

#PCA comparing treatments
PCA_MC = npositions1 %>%
  as_tibble() %>%
  ggplot(aes(x=pcoa1, y=pcoa2, colour=Day, shape=Treatment, group = Day)) +
  geom_point(size=4) + 
  stat_ellipse(linetype=2, size=0.5) +
  theme_bw() +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"))+
  theme_bw() +labs(x="PCo 1 (37.13%)", y="PCo 2 (19.43%)") +
  scale_shape_manual(values=c(16,17)) 

PCA_MC 

# should add direction arrows for specific taxa groups (need PCA loadings for each taxa name)
# to do


#plot all PCAs next to eachother
PCA_DOM
PCA_VOC
PCA_MC

library(ggpubr)

# Arrange the three PCA plots in a single figure
combined_PCA <- ggarrange(PCA_DOM + theme(legend.position = "none"), PCA_VOC + theme(legend.position = "none"), PCA_MC,
                          ncol = 3, nrow = 1,  # Arrange in one row
                          labels = c("A", "B", "C"))  # Add labels if needed

# Display the combined plot
print(combined_PCA)

#### end ####




