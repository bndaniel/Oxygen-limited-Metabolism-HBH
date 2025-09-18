# R code for processing and analyzing metabolites for "Persistence of Oxygen Dependent Dissolved Organic Matter (ODDOM) induced by Oxygen-limited Metabolism in Hypoxic Marine Environments"
# Scripts adopted from https://github.com/Functional-Metabolomics-Lab/CCE_Data-Analysis/blob/main/Notebooks_CCE_Data/StackedBarPlot_Metabolites.ipynb

#### quick look at data ####
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("highcharter")) install.packages("highcharter")

#loading libraries
library(tidyverse)
library(highcharter)

setwd("Desktop/HBH_DOM")

print(list.files("."))

#load in data and metadata
FeatureTable <- read.csv('DOM_matrix.csv',check.names = T,header=T,row.names = 1)
Metadata <- read.csv('metadata.csv',check.names=F,header=T,row.names=1)

#whats going on in there?
head(FeatureTable)
dim(FeatureTable)

#making the metadata a factor
Metadata$d_t <- as.factor(Metadata$d_t) #converting to factor level
date_treatment <- levels(Metadata$d_t)
date_treatment

Metadata$date <- as.factor(Metadata$date) #converting to factor level
dates <- levels(Metadata$date)
dates

Metadata$treatment <- as.factor(Metadata$treatment) #converting to factor level
treatment <- levels(Metadata$treatment)
treatment

#making superclass a factor
FeatureTable$CANOPUS_ClassyFire.superclass <- as.factor(FeatureTable$CANOPUS_ClassyFire.superclass) #converting to factor level
levels(FeatureTable$CANOPUS_ClassyFire.superclass)

#subset by treatment
hft <- subset(FeatureTable, select = c(CANOPUS_ClassyFire.superclass, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
sft <- subset(FeatureTable, select = c(CANOPUS_ClassyFire.superclass, X0_S1a, X0_S1b, X0_S2a, X0_S2b, X0_S3a, X0_S3b, X12_S1a, X12_S2a, X12_S2b, X12_S3a, X12_S3b, X42_S1a, X42_S1b, X42_S2a, X42_S2b, X42_S3a, X42_S3b, X70_S1a, X70_S1b, X70_S2a, X70_S2b, X70_S3a, X70_S3b, X96_S1a, X96_S1b, X96_S2a, X96_S2b, X96_S3a, X96_S3b))

head(hft)
head(sft)
#### end ####

#### H_treatment_DOM_analysis ####
#create an average column based on the mean count for each compound
hft[,2:ncol(hft)] <- lapply(hft[,2:ncol(hft)],as.numeric)
hft$Avg0_H <- rowMeans(hft[,2:7])
hft$Avg12_H <- rowMeans(hft[,8:13])
hft$Avg42_H <- rowMeans(hft[,14:19])
hft$Avg70_H <- rowMeans(hft[,20:25])
hft$Avg96_H <- rowMeans(hft[,26:31])
head(hft)

#creating average column for each superclass specific to each sample
#Avg0_H
Aggregates_Level_Avg0_H<- aggregate(hft$Avg0_H, #Variable to be grouped
                             by=list(hft$CANOPUS_ClassyFire.superclass), #grouping element
                             FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg0_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg0_H$Group_sum)){
  x <- round((Aggregates_Level_Avg0_H$Group_sum[j]/sum(Aggregates_Level_Avg0_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg0_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg0_H$group <- 'Avg0_H'

#Avg12_H
Aggregates_Level_Avg12_H<- aggregate(hft$Avg12_H, #Variable to be grouped
                                    by=list(hft$CANOPUS_ClassyFire.superclass), #grouping element
                                    FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg12_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg12_H$Group_sum)){
  x <- round((Aggregates_Level_Avg12_H$Group_sum[j]/sum(Aggregates_Level_Avg12_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg12_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg12_H$group <- 'Avg12_H'

head(Aggregates_Level_Avg12_H)

#Avg42_H
Aggregates_Level_Avg42_H<- aggregate(hft$Avg42_H, #Variable to be grouped
                                     by=list(hft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg42_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg42_H$Group_sum)){
  x <- round((Aggregates_Level_Avg42_H$Group_sum[j]/sum(Aggregates_Level_Avg42_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg42_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg42_H$group <- 'Avg42_H'

head(Aggregates_Level_Avg42_H)

#Avg70_H
Aggregates_Level_Avg70_H<- aggregate(hft$Avg70_H, #Variable to be grouped
                                     by=list(hft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg70_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg70_H$Group_sum)){
  x <- round((Aggregates_Level_Avg70_H$Group_sum[j]/sum(Aggregates_Level_Avg70_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg70_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg70_H$group <- 'Avg70_H'

head(Aggregates_Level_Avg70_H)

#Avg96_H
Aggregates_Level_Avg96_H<- aggregate(hft$Avg96_H, #Variable to be grouped
                                     by=list(hft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg96_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg96_H$Group_sum)){
  x <- round((Aggregates_Level_Avg96_H$Group_sum[j]/sum(Aggregates_Level_Avg96_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg96_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg96_H$group <- 'Avg96_H'

head(Aggregates_Level_Avg96_H)
Overall_percent_H = data.frame()
Overall_percent_H = rbind(Overall_percent_H, Aggregates_Level_Avg0_H, Aggregates_Level_Avg12_H, Aggregates_Level_Avg42_H, Aggregates_Level_Avg70_H, Aggregates_Level_Avg96_H)
Overall_percent_H

write.csv(Overall_percent_H,"percents_H.csv")

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F","#FFFF99",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","green","#FB9A99","#CAB2D6",
                  "#E31A1C","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Overall_percent_H, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot

#### end ####

#### S_treatment_DOM_analysis ####
#majority of below is just for S treatment
#create an average column based on the mean count for each compound
sft[,2:ncol(sft)] <- lapply(sft[,2:ncol(sft)],as.numeric)
sft$Avg0_S <- rowMeans(sft[,2:7])
sft$Avg12_S <- rowMeans(sft[,8:13])
sft$Avg42_S <- rowMeans(sft[,14:19])
sft$Avg70_S <- rowMeans(sft[,20:25])
sft$Avg96_S <- rowMeans(sft[,26:31])
head(sft)

#creating average column for each superclass specific to each sample
#Avg0_S
Aggregates_Level_Avg0_S<- aggregate(sft$Avg0_S, #Variable to be grouped
                                    by=list(sft$CANOPUS_ClassyFire.superclass), #grouping element
                                    FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg0_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg0_S$Group_sum)){
  x <- round((Aggregates_Level_Avg0_S$Group_sum[j]/sum(Aggregates_Level_Avg0_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg0_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg0_S$group <- 'Avg0_S'

#Avg12_S
Aggregates_Level_Avg12_S<- aggregate(sft$Avg12_S, #Variable to be grouped
                                     by=list(sft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg12_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg12_S$Group_sum)){
  x <- round((Aggregates_Level_Avg12_S$Group_sum[j]/sum(Aggregates_Level_Avg12_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg12_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg12_S$group <- 'Avg12_S'

head(Aggregates_Level_Avg12_S)

#Avg42_S
Aggregates_Level_Avg42_S<- aggregate(sft$Avg42_S, #Variable to be grouped
                                     by=list(sft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg42_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg42_S$Group_sum)){
  x <- round((Aggregates_Level_Avg42_S$Group_sum[j]/sum(Aggregates_Level_Avg42_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg42_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg42_S$group <- 'Avg42_S'

head(Aggregates_Level_Avg42_S)

#Avg70_S
Aggregates_Level_Avg70_S<- aggregate(sft$Avg70_S, #Variable to be grouped
                                     by=list(sft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg70_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg70_S$Group_sum)){
  x <- round((Aggregates_Level_Avg70_S$Group_sum[j]/sum(Aggregates_Level_Avg70_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg70_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg70_S$group <- 'Avg70_S'

head(Aggregates_Level_Avg70_S)

#Avg96_S
Aggregates_Level_Avg96_S<- aggregate(sft$Avg96_S, #Variable to be grouped
                                     by=list(sft$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg96_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg96_S$Group_sum)){
  x <- round((Aggregates_Level_Avg96_S$Group_sum[j]/sum(Aggregates_Level_Avg96_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg96_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg96_S$group <- 'Avg96_S'

head(Aggregates_Level_Avg96_S)

Overall_percent_S = data.frame()
Overall_percent_S = rbind(Overall_percent_S, Aggregates_Level_Avg0_S, Aggregates_Level_Avg12_S, Aggregates_Level_Avg42_S, Aggregates_Level_Avg70_S, Aggregates_Level_Avg96_S)
Overall_percent_S

write.csv(Overall_percent_S,"percents_S.csv", row.namaes=F)

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F","#FFFF99",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","green","#FB9A99","#CAB2D6",
                  "#E31A1C","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Overall_percent_S, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Saturated") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot

#### end ####

#### only differentially abundant compounds and their superclasses ####
head(FeatureTable)

#need to seperate based on ttest < 0.05 
ftsub = FeatureTable
ftsub$ttest = lapply(ftsub$ttest, as.numeric)
ftsub = subset(ftsub, ttest<0.05)
head(ftsub)

#only inculde more abundant in hypoxic in h dataset and more abundant in saturated in s dataset
ftsub_h = subset(ftsub, ftsub$Higher.in.H.or.S. == "Hypoxic")
ftsub_s = subset(ftsub, ftsub$Higher.in.H.or.S. == "Saturated")

#subset by treatment
hft_sig <- subset(ftsub_h, select = c(CANOPUS_ClassyFire.superclass, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
sft_sig <- subset(ftsub_s, select = c(CANOPUS_ClassyFire.superclass, X0_S1a, X0_S1b, X0_S2a, X0_S2b, X0_S3a, X0_S3b, X12_S1a, X12_S1b, X12_S2a, X12_S3a, X12_S3b, X42_S1a, X42_S1b, X42_S2a, X42_S2b, X42_S3a, X42_S3b, X70_S1a, X70_S1b, X70_S2a, X70_S2b, X70_S3a, X70_S3b, X96_S1a, X96_S1b, X96_S2a, X96_S2b, X96_S3a, X96_S3b))

head(hft_sig)
head(sft_sig)

#majority of below is just for H treatment
#create an average column based on the mean count for each compound
hft_sig[,2:ncol(hft_sig)] <- lapply(hft_sig[,2:ncol(hft_sig)],as.numeric)
hft_sig$Avg0_H <- rowMeans(hft_sig[,2:7])
hft_sig$Avg12_H <- rowMeans(hft_sig[,8:13])
hft_sig$Avg42_H <- rowMeans(hft_sig[,14:19])
hft_sig$Avg70_H <- rowMeans(hft_sig[,20:25])
hft_sig$Avg96_H <- rowMeans(hft_sig[,26:31])
head(hft_sig)

#creating average column for each superclass specific to each sample
#Avg0_H
Aggregates_Level_Avg0_H<- aggregate(hft_sig$Avg0_H, #Variable to be grouped
                                    by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                    FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg0_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg0_H$Group_sum)){
  x <- round((Aggregates_Level_Avg0_H$Group_sum[j]/sum(Aggregates_Level_Avg0_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg0_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg0_H$group <- 'Avg0_H'

#Avg12_H
Aggregates_Level_Avg12_H<- aggregate(hft_sig$Avg12_H, #Variable to be grouped
                                     by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg12_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg12_H$Group_sum)){
  x <- round((Aggregates_Level_Avg12_H$Group_sum[j]/sum(Aggregates_Level_Avg12_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg12_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg12_H$group <- 'Avg12_H'

head(Aggregates_Level_Avg12_H)

#Avg42_H
Aggregates_Level_Avg42_H<- aggregate(hft_sig$Avg42_H, #Variable to be grouped
                                     by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg42_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg42_H$Group_sum)){
  x <- round((Aggregates_Level_Avg42_H$Group_sum[j]/sum(Aggregates_Level_Avg42_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg42_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg42_H$group <- 'Avg42_H'

head(Aggregates_Level_Avg42_H)

#Avg70_H
Aggregates_Level_Avg70_H<- aggregate(hft_sig$Avg70_H, #Variable to be grouped
                                     by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg70_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg70_H$Group_sum)){
  x <- round((Aggregates_Level_Avg70_H$Group_sum[j]/sum(Aggregates_Level_Avg70_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg70_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg70_H$group <- 'Avg70_H'

head(Aggregates_Level_Avg70_H)

#Avg96_H
Aggregates_Level_Avg96_H<- aggregate(hft_sig$Avg96_H, #Variable to be grouped
                                     by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg96_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg96_H$Group_sum)){
  x <- round((Aggregates_Level_Avg96_H$Group_sum[j]/sum(Aggregates_Level_Avg96_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg96_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg96_H$group <- 'Avg96_H'

head(Aggregates_Level_Avg96_H)
Overall_percent_H = data.frame()
Overall_percent_H = rbind(Overall_percent_H, Aggregates_Level_Avg0_H, Aggregates_Level_Avg12_H, Aggregates_Level_Avg42_H, Aggregates_Level_Avg70_H, Aggregates_Level_Avg96_H)
Overall_percent_H

write.csv(Overall_percent_H,"percents_H_sig.csv")

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","#FB9A99","#CAB2D6","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Overall_percent_H, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic_sig") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot




#majority of below is just for S treatment
#create an average column based on the mean count for each compound
sft_sig[,2:ncol(sft_sig)] <- lapply(sft_sig[,2:ncol(sft_sig)],as.numeric)
sft_sig$Avg0_S <- rowMeans(sft_sig[,2:7])
sft_sig$Avg12_S <- rowMeans(sft_sig[,8:13])
sft_sig$Avg42_S <- rowMeans(sft_sig[,14:19])
sft_sig$Avg70_S <- rowMeans(sft_sig[,20:25])
sft_sig$Avg96_S <- rowMeans(sft_sig[,26:31])
head(sft_sig)

#creating average column for each superclass specific to each sample
#Avg0_S
Aggregates_Level_Avg0_S<- aggregate(sft_sig$Avg0_S, #Variable to be grouped
                                    by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                    FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg0_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg0_S$Group_sum)){
  x <- round((Aggregates_Level_Avg0_S$Group_sum[j]/sum(Aggregates_Level_Avg0_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg0_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg0_S$group <- 'Avg0_S'

#Avg12_S
Aggregates_Level_Avg12_S<- aggregate(sft_sig$Avg12_S, #Variable to be grouped
                                     by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg12_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg12_S$Group_sum)){
  x <- round((Aggregates_Level_Avg12_S$Group_sum[j]/sum(Aggregates_Level_Avg12_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg12_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg12_S$group <- 'Avg12_S'

head(Aggregates_Level_Avg12_S)

#Avg42_S
Aggregates_Level_Avg42_S<- aggregate(sft_sig$Avg42_S, #Variable to be grouped
                                     by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg42_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg42_S$Group_sum)){
  x <- round((Aggregates_Level_Avg42_S$Group_sum[j]/sum(Aggregates_Level_Avg42_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg42_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg42_S$group <- 'Avg42_S'

head(Aggregates_Level_Avg42_S)

#Avg70_S
Aggregates_Level_Avg70_S<- aggregate(sft_sig$Avg70_S, #Variable to be grouped
                                     by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg70_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg70_S$Group_sum)){
  x <- round((Aggregates_Level_Avg70_S$Group_sum[j]/sum(Aggregates_Level_Avg70_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg70_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg70_S$group <- 'Avg70_S'

head(Aggregates_Level_Avg70_S)

#Avg96_S
Aggregates_Level_Avg96_S<- aggregate(sft_sig$Avg96_S, #Variable to be grouped
                                     by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg96_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg96_S$Group_sum)){
  x <- round((Aggregates_Level_Avg96_S$Group_sum[j]/sum(Aggregates_Level_Avg96_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg96_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg96_S$group <- 'Avg96_S'

head(Aggregates_Level_Avg96_S)

Overall_percent_S = data.frame()
Overall_percent_S = rbind(Overall_percent_S, Aggregates_Level_Avg0_S, Aggregates_Level_Avg12_S, Aggregates_Level_Avg42_S, Aggregates_Level_Avg70_S, Aggregates_Level_Avg96_S)
Overall_percent_S

write.csv(Overall_percent_S,"percents_S_sig.csv")

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","#FB9A99","#CAB2D6",
                  "#E31A1C","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Overall_percent_S, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Saturated_sig") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot
#### end ####

##### aggregate based on only treatment and not days ####
hft_sig$Avg_H <- rowMeans(hft_sig[,2:31])
head(hft_sig)

Aggregates_Level_Avg_H<- aggregate(hft_sig$Avg_H, #Variable to be grouped
                                     by=list(hft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg_H$Group_sum)){
  x <- round((Aggregates_Level_Avg_H$Group_sum[j]/sum(Aggregates_Level_Avg_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg_H$group <- 'Avg_H'

head(Aggregates_Level_Avg_H)

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","#FB9A99","#CAB2D6","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Aggregates_Level_Avg_H, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="average_hypoxic_sig") +
  xlab("Hypoxic") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot


#for saturated
sft_sig$Avg_S <- rowMeans(sft_sig[,2:31])
head(sft_sig)

Aggregates_Level_Avg_S<- aggregate(sft_sig$Avg_S, #Variable to be grouped
                                   by=list(sft_sig$CANOPUS_ClassyFire.superclass), #grouping element
                                   FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg_S) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg_S$Group_sum)){
  x <- round((Aggregates_Level_Avg_S$Group_sum[j]/sum(Aggregates_Level_Avg_S$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg_S$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg_S$group <- 'Avg_S'

head(Aggregates_Level_Avg_S)

#making the stack plot!
stack_colors <- c("#6A3D9A","#FDBF6F",'#660066',
                  "#B2DF8A","gray","darkgreen","#85c0ed",
                  "#1f78b4","#FB9A99","#CAB2D6",
                  "#E31A1C","#aa11aa","#B15928","#FF7F00")

stackPlot <- ggplot(Aggregates_Level_Avg_S, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="average_saturated_sig") +
  xlab("Saturated") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = stack_colors)
stackPlot

#### end ####

#### Zooming in on significantly different benzenoids to higher class levels... in hypoxic treatment ####

#making specific class a factor
FeatureTable$CANOPUS_ClassyFire.most.specific.class <- as.factor(FeatureTable$CANOPUS_ClassyFire.most.specific.class) #converting to factor level
levels(FeatureTable$CANOPUS_ClassyFire.most.specific.class)

#need to seperate based on ttest < 0.05 
ftsub = FeatureTable
ftsub$ttest = lapply(ftsub$ttest, as.numeric)
ftsub = subset(ftsub, ttest<0.05)
head(ftsub)

#only inculde more abundant in hypoxic in h dataset and more abundant in saturated in s dataset
ftsub_h = subset(ftsub, ftsub$Higher.in.H.or.S. == "Hypoxic")

#subset by treatment and being benzenoids
hft_ben <- subset(ftsub_h, select = c(CANOPUS_ClassyFire.superclass, CANOPUS_ClassyFire.most.specific.class, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
hft_ben2 = subset(hft_ben, hft_ben$CANOPUS_ClassyFire.superclass == "Benzenoids")
hft_ben2 = subset(hft_ben2, select = c(CANOPUS_ClassyFire.most.specific.class, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
head(hft_ben2)


#create an average column based on the mean count for each compound
hft_ben2[,2:ncol(hft_ben2)] <- lapply(hft_ben2[,2:ncol(hft_ben2)],as.numeric)
hft_ben2$Avg0_H <- rowMeans(hft_ben2[,2:7])
hft_ben2$Avg12_H <- rowMeans(hft_ben2[,8:13])
hft_ben2$Avg42_H <- rowMeans(hft_ben2[,14:19])
hft_ben2$Avg70_H <- rowMeans(hft_ben2[,20:25])
hft_ben2$Avg96_H <- rowMeans(hft_ben2[,26:31])
head(hft_ben2)

#creating average column for each superclass specific to each sample
#Avg0_H
Aggregates_Level_Avg0_H<- aggregate(hft_ben2$Avg0_H, #Variable to be grouped
                                    by=list(hft_ben2$CANOPUS_ClassyFire.most.specific.class), #grouping element
                                    FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg0_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg0_H$Group_sum)){
  x <- round((Aggregates_Level_Avg0_H$Group_sum[j]/sum(Aggregates_Level_Avg0_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg0_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg0_H$group <- 'Avg0_H'

#Avg12_H
Aggregates_Level_Avg12_H<- aggregate(hft_ben2$Avg12_H, #Variable to be grouped
                                     by=list(hft_ben2$CANOPUS_ClassyFire.most.specific.class), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg12_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg12_H$Group_sum)){
  x <- round((Aggregates_Level_Avg12_H$Group_sum[j]/sum(Aggregates_Level_Avg12_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg12_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg12_H$group <- 'Avg12_H'

head(Aggregates_Level_Avg12_H)

#Avg42_H
Aggregates_Level_Avg42_H<- aggregate(hft_ben2$Avg42_H, #Variable to be grouped
                                     by=list(hft_ben2$CANOPUS_ClassyFire.most.specific.class), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg42_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg42_H$Group_sum)){
  x <- round((Aggregates_Level_Avg42_H$Group_sum[j]/sum(Aggregates_Level_Avg42_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg42_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg42_H$group <- 'Avg42_H'

head(Aggregates_Level_Avg42_H)

#Avg70_H
Aggregates_Level_Avg70_H<- aggregate(hft_ben2$Avg70_H, #Variable to be grouped
                                     by=list(hft_ben2$CANOPUS_ClassyFire.most.specific.class), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg70_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg70_H$Group_sum)){
  x <- round((Aggregates_Level_Avg70_H$Group_sum[j]/sum(Aggregates_Level_Avg70_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg70_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg70_H$group <- 'Avg70_H'

head(Aggregates_Level_Avg70_H)

#Avg96_H
Aggregates_Level_Avg96_H<- aggregate(hft_ben2$Avg96_H, #Variable to be grouped
                                     by=list(hft_ben2$CANOPUS_ClassyFire.most.specific.class), #grouping element
                                     FUN=sum) #provides the sum of all Avg
colnames(Aggregates_Level_Avg96_H) <- c("Group_name","Group_sum")
Percent <- c() 
for(j in 1:length(Aggregates_Level_Avg96_H$Group_sum)){
  x <- round((Aggregates_Level_Avg96_H$Group_sum[j]/sum(Aggregates_Level_Avg96_H$Group_sum))*100,3) 
  Percent <- c(Percent,x)
}
Aggregates_Level_Avg96_H$Percent <- Percent  #including the percent contribution to Aggregates Level dataframe
Aggregates_Level_Avg96_H$group <- 'Avg96_H'

head(Aggregates_Level_Avg96_H)
Overall_percent_H = data.frame()
Overall_percent_H = rbind(Overall_percent_H, Aggregates_Level_Avg0_H, Aggregates_Level_Avg12_H, Aggregates_Level_Avg42_H, Aggregates_Level_Avg70_H, Aggregates_Level_Avg96_H)
Overall_percent_H

write.csv(Overall_percent_H,"benzene_classes.csv")


#making the stack plot!

stackPlot <- ggplot(Overall_percent_H, aes(fill=Group_name, y=Percent, x=group)) + 
  geom_bar(position="stack", stat="identity")+ 
  ggtitle(label="Hypoxic_sig") +
  xlab("Day") + ylab("Percentage") + labs(fill = "Superclass") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  #theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5)) 
stackPlot
#### end ####

#### create a list of all the phenol ether compounds ####

hft_phenol_ethers = subset(hft_ben2, hft_ben2$CANOPUS_ClassyFire.most.specific.class == "Phenol ethers")
head(hft_phenol_ethers)

write.csv(hft_phenol_ethers,"phenol_ethers.csv")

#subset by treatment and being benzenoids
head(ftsub_h)
#hft_ben <- subset(ftsub_h, select = c(row.m.z, SIRIUS_molecularFormula, CANOPUS_ClassyFire.superclass, CANOPUS_ClassyFire.most.specific.class, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
hft_phenol = subset(ftsub_h, hft_ben$CANOPUS_ClassyFire.most.specific.class == "Phenol ethers")
hft_phenol2 = subset(hft_phenol, select = c(row.m.z, SIRIUS_molecularFormula, CANOPUS_molecularFormula, CANOPUS_ClassyFire.superclass, CANOPUS_ClassyFire.most.specific.class, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))
head(hft_phenol2)

hft_phenol2 = subset(hft_phenol, select = c(row.m.z, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b))

#create an average column based on the mean count for each compound
hft_phenol2[,2:ncol(hft_phenol2)] <- lapply(hft_phenol2[,2:ncol(hft_phenol2)],as.numeric)
hft_phenol2$Avg0_H <- rowMeans(hft_phenol2[,2:7])
hft_phenol2$Avg12_H <- rowMeans(hft_phenol2[,8:13])
hft_phenol2$Avg42_H <- rowMeans(hft_phenol2[,14:19])
hft_phenol2$Avg70_H <- rowMeans(hft_phenol2[,20:25])
hft_phenol2$Avg96_H <- rowMeans(hft_phenol2[,26:31])
head(hft_phenol2)

hft_phenol2$SIRIUS_molecularFormula = hft_phenol$SIRIUS_molecularFormula
hft_phenol2$CANOPUS_molecularFormula = hft_phenol$CANOPUS_molecularFormula
hft_phenol2$CANOPUS_ClassyFire.most.specific.class = hft_phenol$CANOPUS_ClassyFire.most.specific.class

head(hft_phenol2)

hft_phenol2 = subset(hft_phenol2, select = c(row.m.z, SIRIUS_molecularFormula, CANOPUS_molecularFormula, CANOPUS_ClassyFire.most.specific.class, X0_H1a, X0_H1b, X0_H2a, X0_H2b, X0_H3a, X0_H3b, X12_H1a, X12_H2a, X12_H2b, X12_H3a, X12_H3b, X42_H1a, X42_H1b, X42_H2a, X42_H2b, X42_H3a, X42_H3b, X70_H1a, X70_H1b, X70_H2a, X70_H2b, X70_H3a, X70_H3b, X96_H1a, X96_H1b, X96_H2a, X96_H2b, X96_H3a, X96_H3b, Avg0_H, Avg12_H, Avg42_H, Avg70_H, Avg96_H))
head(hft_phenol2)

write.csv(hft_phenol2,"phenol_ethers.csv")


# create bar graph for each compound at each timepoint comparing abundance
# did it manually in excel

#### end ####

#### creating van-krevelen plot - separating molecular formula ####
if (!require("ggthemes")) install.packages("ggthemes")
if (!require("svglite")) install.packages("svglite")

library(dplyr)
library(tidyr)
library(stringr)
library(ggthemes)
library(ggplot2)
library(svglite)

head(FeatureTable)

# First, Sirius molecular formulas

#splits the molecular formula if alphabets are present in it such as 'C13' 'H22' 'S'
if (!require("CHNOSZ")) install.packages("CHNOSZ")
library(CHNOSZ)
library(data.table)

formulas <- makeup(FeatureTable$SIRIUS_molecularFormula, count.zero = TRUE)
formulas <- data.table(compound = FeatureTable$SIRIUS_molecularFormula)[
  ,names(formulas[[1]]) := transpose(formulas)
]

head(formulas)

#change order to match ?

sf <- as.data.frame(formulas)
rownames(sf) <- NULL
colnames(sf) <- paste(colnames(sf), "Sirius", sep = '_')

head(sf)

#do the same for canopus

formulas <- makeup(FeatureTable$CANOPUS_molecularFormula, count.zero = TRUE)
formulas <- data.table(compound = FeatureTable$CANOPUS_molecularFormula)[
  ,names(formulas[[1]]) := transpose(formulas)
]

head(formulas)

#change order to match ?

cf <- as.data.frame(formulas) #canopus formula
rownames(cf) <- NULL
colnames(cf) <- paste(colnames(cf), "Canopus", sep = '_')

head(cf)

#combine them
fin <- cbind(FeatureTable, sf, cf) #combining both sirius and canopus formulae into a final dataframe
head(fin)

df <- fin
head(df)
#C:N, O:C, H:C and average C oxidation columns using only Sirius formulas
df$CN = df$C_Sirius/df$N_Sirius
df$OC = df$O_Sirius/df$C_Sirius
df$HC = df$H_Sirius/df$C_Sirius
df$avCox <- -((1*df$H_Sirius) - (3*df$N_Sirius) - (2*df$O_Sirius) + (5*df$P_Sirius) -(2*df$S_Sirius))/(df$C_Sirius)
head(df)

#making vk plot 
#Only selecting formulas that make sense molecularly (chemically)
df_vk <- df %>% filter(avCox > -4 & avCox < 4 & OC < 2 & HC < 2.5 & CN > 4)
write.csv(df_vk,"Sirius_Canopus_elemental_info.csv")
head(df_vk)
range(df_vk$avCox)

#check out dat big data frame
dim(df_vk)
colnames(df_vk)


#### end ####

#### making VK plots: ####
# Relative abundance of Carbon
VK <- ggplot(df_vk, aes(x = OC, y = HC, col=C_Sirius)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "Relative Abundance of Carbon",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(df_vk))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 
#svglite("VK_Plot_Relative Abundance_C.svg")
VK
dev.off()

# based on higher in h or s
VK2 <- ggplot(df_vk, aes(x = OC, y = HC, col=`CANOPUS_ClassyFire.superclass`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(df_vk))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5))

#svglite("Van_Krevelen_Plot_Sirius_formula_withoutNA_with_SuperClass_Lables.svg")
VK2
dev.off()


# based on higher in h or s but only looking at differentially abundant
head(df_vk)
df_vk2 = subset(df_vk, ttest<0.05)
dim(df_vk2)

VK3 <- ggplot(df_vk2, aes(x = OC, y = HC, col=`Higher.in.H.or.S.`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(df_vk2))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) +
  scale_color_manual('Treatment',, values = stack_colors) +
  stat_ellipse()

VK3
dev.off()

# zoom in on it a bit
VK4 <- ggplot(df_vk2, aes(x = OC, y = HC, col=`Higher.in.H.or.S.`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(df_vk2))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) +
  scale_color_manual('Treatment',, values = stack_colors) + 
  stat_ellipse()

VK4

# van-krevelen of only compunds that are higher abundant in hypoxic treatment overall (based on averages) and increasing at a constant rate (correlation test)
corr_H <- read.csv('corr_H.csv',check.names = T,header=T,row.names = 1)
head(corr_H)

VK5 <- ggplot(corr_H, aes(x = OC, y = HC, col=`P.value`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(corr_H))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) +
  stat_ellipse()

VK5


# van-krevelen of only compunds that are higher abundant in saturated treatment overall (based on averages) and increasing at a constant rate (correlation test)
corr_S <- read.csv('corr_S.csv',check.names = T,header=T,row.names = 1)
head(corr_S)

VK6 <- ggplot(corr_S, aes(x = OC, y = HC, col=`P.value`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(corr_S))) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) +
  stat_ellipse()

VK6

#### end ####

#### plot AI vs Cox ####
colnames(corr_H)
VK7 <- ggplot(corr_S, aes(x = avCox, y = AI, col=`P.value`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(corr_S))) +
  labs(x = "avCox", y = "AI") +
  scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-3, 1, by = 0.3)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0.0, 40, by = 5)) +
  stat_ellipse()

VK7


VK8 <- ggplot(corr_H, aes(x = avCox, y = AI, col=`P.value`)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('Van Krevelen Plot of Sirius Formulas:',nrow(corr_H))) +
  labs(x = "avCox", y = "AI") +
  scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-3, 1, by = 0.3)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0.0, 40, by = 5)) +
  stat_ellipse()

VK8

#### end ####

#### Van-krevel of each treatment over time ####
# use this method for boxplot of AI and avCox comparison
head(FeatureTable)
dim(FeatureTable)
colnames(FeatureTable)

#create an average column based on the mean count for each compound
vk_time = FeatureTable
head(vk_time)
vk_time$Avg0_H <- rowMeans(vk_time[,8:13])
vk_time$Avg12_H <- rowMeans(vk_time[,14:19])
vk_time$Avg42_H <- rowMeans(vk_time[,20:25])
vk_time$Avg70_H <- rowMeans(vk_time[,26:31])
vk_time$Avg96_H <- rowMeans(vk_time[,32:37])
vk_time$Avg0_S <- rowMeans(vk_time[38:43])
vk_time$Avg12_S <- rowMeans(vk_time[44:49])
vk_time$Avg42_S <- rowMeans(vk_time[50:55])
vk_time$Avg70_S <- rowMeans(vk_time[56:61])
vk_time$Avg96_S <- rowMeans(vk_time[62:67])
head(vk_time)

# convert abundance values into relative abundance 
vk_time$percab0_H <- (vk_time$Avg0_H)/sum(vk_time$Avg0_H)
vk_time$percab12_H <- (vk_time$Avg12_H)/sum(vk_time$Avg12_H)
vk_time$percab42_H <- (vk_time$Avg42_H)/sum(vk_time$Avg42_H)
vk_time$percab70_H <- (vk_time$Avg70_H)/sum(vk_time$Avg70_H)
vk_time$percab96_H <- (vk_time$Avg96_H)/sum(vk_time$Avg96_H)
vk_time$percab0_S <- (vk_time$Avg0_S)/sum(vk_time$Avg0_S)
vk_time$percab12_S <- (vk_time$Avg12_S)/sum(vk_time$Avg12_S)
vk_time$percab42_S <- (vk_time$Avg42_S)/sum(vk_time$Avg42_S)
vk_time$percab70_S <- (vk_time$Avg70_S)/sum(vk_time$Avg70_S)
vk_time$percab96_S <- (vk_time$Avg96_S)/sum(vk_time$Avg96_S)
head(vk_time)


#Rank order method (not looking at top 50% - all compounds but highlighting with rank order)
vk_time$rank0_H = rank(vk_time$percab0_H)
vk_time$rank12_H = rank(vk_time$percab12_H)
vk_time$rank42_H = rank(vk_time$percab42_H)
vk_time$rank70_H = rank(vk_time$percab70_H)
vk_time$rank96_H = rank(vk_time$percab96_H)
vk_time$rank0_S = rank(vk_time$percab0_S)
vk_time$rank12_S = rank(vk_time$percab12_S)
vk_time$rank42_S = rank(vk_time$percab42_S)
vk_time$rank70_S = rank(vk_time$percab70_S)
vk_time$rank96_S = rank(vk_time$percab96_S)
head(vk_time)

#vkplots
#splits the molecular formula if alphabets are present in it such as 'C13' 'H22' 'S'
if (!require("CHNOSZ")) install.packages("CHNOSZ")
library(CHNOSZ)
library(data.table)

formulas <- makeup(vk_time$SIRIUS_molecularFormula, count.zero = TRUE)
formulas <- data.table(compound = vk_time$SIRIUS_molecularFormula)[
  ,names(formulas[[1]]) := transpose(formulas)
]

head(formulas)

#change order to match ?

sf <- as.data.frame(formulas)
rownames(sf) <- NULL
colnames(sf) <- paste(colnames(sf), "Sirius", sep = '_')

head(sf)
#combine them
fin <- cbind(vk_time, sf) #combining both sirius and canopus formulae into a final dataframe
head(fin)


df <- fin
head(df)
#C:N, O:C, H:C and average C oxidation columns using only Sirius formulas
df$CN = df$C_Sirius/df$N_Sirius
df$OC = df$O_Sirius/df$C_Sirius
df$HC = df$H_Sirius/df$C_Sirius
df$avCox <- -((1*df$H_Sirius) - (3*df$N_Sirius) - (2*df$O_Sirius) + (5*df$P_Sirius) -(2*df$S_Sirius))/(df$C_Sirius)
head(df)

#making vk plot 
#Only selecting formulas that make sense molecularly (chemically)
df_vk <- df %>% filter(avCox > -4 & avCox < 4 & OC < 2 & HC < 2.5 & CN > 4)
write.csv(df_vk,"Sirius_Canopus_elemental_info.csv")
head(df_vk)
range(df_vk$avCox)

#check out dat big data frame
dim(df_vk)
colnames(df_vk)

# plotting with rank as main interest
VK <- ggplot(df_vk, aes(x = OC, y = HC, col=rank0_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "rank",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('rank0_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK

VK2 <- ggplot(df_vk, aes(x = OC, y = HC, col=rank12_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "rank",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('rank12_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK2

VK3 <- ggplot(df_vk, aes(x = OC, y = HC, col=rank42_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "rank",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('rank42_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK3

VK4 <- ggplot(df_vk, aes(x = OC, y = HC, col=rank70_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "rank",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('rank70_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK4

VK5 <- ggplot(df_vk, aes(x = OC, y = HC, col=rank96_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "rank",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('rank96_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK5
#### end ####

#### only the compounds that make up the top 50% ####
head(df_vk)
# make new dataframe with only the compounds in each column that makes up the top 50%
top_0_H = df_vk[order(df_vk$percab0_H, decreasing = TRUE),]
top_12_H = df_vk[order(df_vk$percab12_H, decreasing = TRUE),]
top_42_H = df_vk[order(df_vk$percab42_H, decreasing = TRUE),]
top_70_H = df_vk[order(df_vk$percab70_H, decreasing = TRUE),]
top_96_H = df_vk[order(df_vk$percab96_H, decreasing = TRUE),]
top_0_S = df_vk[order(df_vk$percab0_S, decreasing = TRUE),]
top_12_S = df_vk[order(df_vk$percab12_S, decreasing = TRUE),]
top_42_S = df_vk[order(df_vk$percab42_S, decreasing = TRUE),]
top_70_S = df_vk[order(df_vk$percab70_S, decreasing = TRUE),]
top_96_S = df_vk[order(df_vk$percab96_S, decreasing = TRUE),]
head(top_0_H)


# Function definition for finding the column where cumulation of perc abundance equals 50
DOM50 <- function(data) {
  DOM50 <- cumsum(as.numeric(data)) <= 0.5
  yup <- sum(DOM50 == TRUE)
  return(yup)
}

#getting the column number for each dataframe and adding it to the right
DOM50(top_0_H$percab0_H) #347
DOM50(top_12_H$percab12_H) #154
DOM50(top_42_H$percab42_H) #128
DOM50(top_70_H$percab70_H) #177
DOM50(top_96_H$percab96_H) #135
DOM50(top_0_S$percab0_S) #427
DOM50(top_12_S$percab12_S) #76
DOM50(top_42_S$percab42_S) #78
DOM50(top_70_S$percab70_S) #84
DOM50(top_96_S$percab96_S) #129

#you can check to make sure!
sum(top_96_H$percab96_H[0:135]) #should be less than .5

#subset each dataframe to only include the compounds that make up the top 50%
top_0_H = top_0_H[0:347,]
top_12_H = top_12_H[0:154,]
top_42_H = top_42_H[0:128,]
top_70_H = top_70_H[0:177,]
top_96_H = top_96_H[0:135,]
top_0_S = top_0_S[0:427,]
top_12_S = top_12_S[0:76,]
top_42_S = top_42_S[0:78,]
top_70_S = top_70_S[0:84,]
top_96_S = top_96_S[0:129,]

#make a vk for each
head(top_0_H)
VK1 <- ggplot(top_0_H, aes(x = OC, y = HC, col=percab0_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "percent abundance",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('top_0_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK1

VK2 <- ggplot(top_12_H, aes(x = OC, y = HC, col=percab12_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "percent abundance",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('top_12_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK2

VK3 <- ggplot(top_42_H, aes(x = OC, y = HC, col=percab42_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "percent abundance",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('top_42_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK3

VK4 <- ggplot(top_70_H, aes(x = OC, y = HC, col=percab70_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "percent abundance",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('top_70_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK4

VK5 <- ggplot(top_96_H, aes(x = OC, y = HC, col=percab96_H)) +
  geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
  scale_color_gradient(name = "percent abundance",low = ("#B2DF8A"), high = ("black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#eeeeee"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"),
        plot.title = element_text(size = 16, face = "bold",hjust=0.5),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle(paste('top_96_H')) +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5)) 

VK5

#make a dataframe of each with only the compound and percent abundance

ya = merge(top_0_H, top_12_H, by= 'row.ID')


DOM50 <- function() {
  DOM50 <- cumsum(as.numeric(top_0_H$percab0_H)) <= 0.5
  yup = nrow(as.data.frame(subset(DOM50, DOM50 == "TRUE")))
  return(yup)
}
DOM50()



# tryna make it an operation with a for loop
top50=list(c('top_0_H$percab0_H')) #need this to be a list of the each from above
DOM50 <- function() {
  DOM50 <- cumsum(as.numeric(top_0_H$percab0_H)) <= 0.5
  yup = nrow(as.data.frame(subset(DOM50, DOM50 == "TRUE")))
  return(yup)
}
DOM50(top_0_H$percab0_H)


# Function definition
DOM50 <- function(data) {
  DOM50 <- cumsum(as.numeric(data)) <= 0.5
  yup <- sum(DOM50 == TRUE)
  return(yup)
}

# List of variables
variable_list <- c("top_0_H$percab0_H", "top_12_H$percab12_H", "top_42_H$percab42_H", "top_70_H$percab70_H", "top_96_H$percab96_H")

# Create an empty vector to store results
results <- numeric(length(variable_list))

# Loop through variables and apply the function
for (i in seq_along(variable_list)) {
  results[i] <- DOM50(data = vk_time, variable = variable_list[i])
}

# Print results
print(results)


#### venn diagram of present vs absent taxa in each treatment (from day 12 to end) ####
#details of presence vs absence methods are in "meta16s(autorecovered).xlsx"
install.packages("VennDiagram")
library(VennDiagram)

x <- list(
  "A" = 1:111,
  "B" = 13:130)

venn.diagram(x, filename = "venn-2-dimensions.png", scaled=FALSE)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


display_venn(x, scale=FALSE)

display_venn(
  x,
  category.names = c("Hypoxic" , "Saturated"),
  fill = c("#E69F00", "#56B4E9"), scale=FALSE
)

#### end ####
