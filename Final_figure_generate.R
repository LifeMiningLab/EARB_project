######################################
#####         LIBRARIES          #####
######################################

library(tidyverse)
#library(viridis)
library(dplyr)
library("ggsci")
library("ggplot2")
#library("gridExtra")
library(gtools)
library(data.table)
library(ggpubr)
#library(lattice)
library(ggtext)
library(paletteer)
library(BioCircos)
library(circlize)
library(pheatmap)
#library(reshape2)
library(vegan)

################################################################################

######################################
#####         FIGURE 1           #####
######################################

#### (fig1-a) sample mag quality check ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata") # Please change the directory #
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_MAGs_quality.Rdata") # load name: total_MAGs_stats

# categorize MAGs into 'high' or 'low' MAG
total_MAGs_stats <- total_MAGs_stats %>% mutate("condition" = ifelse(completeness > 70 & contamination < 5, 'HQ', 'LQ'))

table(total_MAGs_stats$condition) # low 5268 / high 2585

# Generating figure 1b
ggscatter(total_MAGs_stats, x ="contamination", y = "completeness",
          color = 'condition',palette = 'npg', xlim = c(0,10), main = "MAG stats", size = 0.3, xlab = 'Contamination (T, ratio)', ylab = 'Completeness (Q, %)', legend.title='Quality (Q > 70, T < 5)') +
  geom_hline(yintercept = 70, linetype = "dashed"  , color = "black") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(colour = guide_legend(title.position = "top"))

rm(list=ls())

#### (fig1-c) number of bins and number of AMR ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata") # Please change the directory #
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

# calculate number of MAGs and AMR genes for each sample
tmp <- n_bin_n_amr %>% group_by(Run) %>% summarise(Bin=n())
tmp2 <- n_bin_n_amr %>% select(Run, AMR) %>% group_by(Run) %>% mutate(AMR=sum(AMR)) %>% unique()

n_bin_n_amr <- tmp %>% left_join(tmp2, by='Run') %>% left_join(metadata, by='Run')

remove(tmp,tmp2)

# calculate AMR per MAG ratio
n_bin_n_amr <- n_bin_n_amr %>%
  mutate('AMR_per_MAG' = `AMR`/`Bin`)

# rename Bin to MAG for clarify
n_bin_n_amr <- n_bin_n_amr %>% rename(MAG=Bin)

# Genereating figure 1c
ggboxplot(data=n_bin_n_amr, x='condition', y=c('MAG','AMR'), combine = T, ylab = 'Number of MAGs and AMR genes',  fill = "condition", palette = "npg", x.text.angle = 45, legend = 'none', outlier.size=0.6, size = 0.1) +
  ylim(c(0,100))+
  rremove('xlab') +
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6)+
  theme(strip.text.x = element_text(size = 6)) +
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "D0", size = 2, label.y = 95)


rm(list=ls())

#### (fig1-d) AMR dynamics of individuals ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

n_bin_n_amr <- n_bin_n_amr %>% left_join(metadata, by='Run') 

n_bin_n_amr$sample <- gsub('ERAS','H',n_bin_n_amr$sample) # change ERAS to H for simplifying

# calculate mean AMR gene possession per MAG in each sample
amr_dynamics <- n_bin_n_amr %>% group_by(Run) %>% mutate('Average_AMR' = mean(AMR)) %>% select(Run, sample, condition, Average_AMR) %>% unique()
amr_dynamics$sample <- factor(amr_dynamics$sample, level=c(unique(mixedsort(amr_dynamics$sample))))

# Generating figure 1d
ggline(amr_dynamics, "condition","Average_AMR", color = "sample", palette = get_palette(palette = "npg", 12), ylab="Total number of AMR genes per total MAG", legend.title='Host') +
  rremove('xlab') +
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.margin=margin(0,0,0,0)) +  # Set top, right, bottom, and left margins
  guides(colour = guide_legend(title.position = "top"))

rm(list=ls())

#### (fig1-e) susceptable and tolerance individuals line plots with arrow ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

n_bin_n_amr <- n_bin_n_amr %>% left_join(metadata, by='Run') 

n_bin_n_amr$sample <- gsub('ERAS','H',n_bin_n_amr$sample) # change ERAS to H for simplifying

amr_dynamics <- n_bin_n_amr %>% group_by(Run) %>% mutate('Average_AMR' = mean(AMR)) %>% select(Run, sample, condition, Average_AMR) # no unique this time
amr_dynamics$sample <- factor(amr_dynamics$sample, level=c(unique(mixedsort(amr_dynamics$sample))))

# calculate number of MAGs and average_AMR genes per MAG
temp <- amr_dynamics %>% unique()
amr_dynamics <- amr_dynamics %>% group_by(Run,sample,condition) %>% summarise(n=n()) %>% left_join(temp, by=c('Run','sample','condition'))

remove(temp)

# select host who response quickly or slowly to antibiotic treatment
amr_dynamics <- amr_dynamics %>% filter(sample %in% c('H2', 'H4', 'H5', 'H12'))

amr_dynamics <- droplevels(amr_dynamics %>% filter(sample %in% c('H2', 'H4', 'H5', 'H12')))
amr_dynamics$condition <- factor(amr_dynamics$condition, levels = c("D0","D4","D8","D42","D180"))

# categorize host into Tolerant or Susceptible
amr_dynamics <- amr_dynamics %>% mutate(susceptability = ifelse(sample %in% c('H2','H4'), 'Tolerant','Susceptible'))

# information dataframe for plotting arrow plot
host_arrows <- amr_dynamics %>%
  arrange(sample, condition) %>%
  group_by(sample) %>%
  mutate(xend = lead(n),
         yend = lead(Average_AMR)) %>%
  ungroup() %>%
  filter(!is.na(xend) & !is.na(yend))

# Generating figure 1e #
ggplot() +
  geom_segment(data = host_arrows, 
               aes(x = n, y = Average_AMR, xend = xend, yend = yend, color = sample), 
               arrow = arrow(type = "closed", length = unit(0.2, "cm")),
               linewidth = 0.5) +
  geom_point(data = amr_dynamics, 
             aes(x = n, y = Average_AMR, shape = condition, fill = condition), 
             pch = 21) +
  theme_bw() +
  facet_wrap(~ susceptability) +
  scale_color_manual(values = get_palette(palette = "npg", '9')) +
  scale_shape_manual(values = 1:length(unique(amr_dynamics$condition))) +
  labs(x = "Number of MAGs",
       y = "Average number of AMR genes per MAG",
       shape = "Time point",  # setting legend title for shape
       color = "Host") +      # setting legend title for color
  theme(
    text = element_text(size = 6),
    legend.key.size = unit(0.3, 'cm'),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold")
  ) +
  guides(
    color = guide_legend(title.position = "top"),
    shape = guide_legend(title.position = "top")
  )

rm(list = ls())

#### (fig1-f) number of AMR in MAG ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

n_bin_n_amr <- n_bin_n_amr %>% left_join(metadata, by='Run') 

# collect EARB bin, which possess more than 17 AMR genes
EARB_bin <- n_bin_n_amr %>%
  filter(AMR >= 17) 

# Generating figure 1f
ggboxplot(n_bin_n_amr, x='condition', y='AMR', fill = 'condition', legend = 'none', palette = 'npg',outlier.size=0.6, size = 0.1, ylab = 'Number of AMR genes in MAGs') +
  ylim(c(0,70)) +
  geom_hline(yintercept = 17, linetype = 2) +
  rremove('xlab') +
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6)+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "D0",label.y = 65, size=3) + 
  geom_point(data = EARB_bin, size=0.6, color='red') 

rm(list = ls())


################################################################################

######################################
#####         FIGURE 2           ##### 
######################################

#### ultimate information of AMR rgi result ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

EARB_list <- n_bin_n_amr %>% filter(AMR>16) %>% select(Run, Bin)
none_EARB_list <- n_bin_n_amr %>% filter(AMR<16 & AMR>0) %>% select(Run, Bin)

amr_gene_rgi_result %>% select(Run,Bin) %>% inner_join(EARB_list, by=c('Run','Bin')) # 803 EARB
amr_gene_rgi_result %>% select(Run,Bin) %>% inner_join(none_EARB_list, by=c('Run','Bin')) # 832 none_EARB

rm(list=ls())

#### (fig2-a) EARB and SARB relative frequency analysis ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

amr_gene_rgi_result <- amr_gene_rgi_result %>% select(Run, Bin, `Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class
total_AMR_gene <- amr_gene_rgi_result %>% select(Run, Bin) %>% group_by(Run, Bin) %>% summarise(total_AMR_gene=n())

# Create an empty data frame with the columns
amr_gene_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(amr_gene_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(amr_gene_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = amr_gene_rgi_result$Run[i], Bin = amr_gene_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    amr_gene_drugclass <- bind_rows(amr_gene_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(amr_gene_drugclass) # number of row: 4333

amr_gene_drugclass <- amr_gene_drugclass %>% left_join(total_AMR_gene, by=c('Run', 'Bin'))
rm(total_AMR_gene)

# calculate number of drugclass in each MAGs
total_drugclass <- amr_gene_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup() # number of row: 1508

EARB_list <- n_bin_n_amr %>% filter(AMR > 16) %>% select(Run, Bin)
SARB_list <- n_bin_n_amr %>% filter(AMR < 16 & AMR > 0) %>% select(Run, Bin)

## EARB ##
EARB_drugclass <- total_drugclass %>% inner_join(EARB_list, by=c('Run', 'Bin'))
sum(EARB_drugclass$total_drugclass) # 2889 AMR drugclass

#number of carriers#
EARB_drugclass_carrier <- EARB_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n())

EARB_drugclass_carrier %>% print(n=23) # total number of drug classes: 23

#rel_frequency calculate#
EARB_drugclass_sum <- EARB_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount' = sum(total_drugclass))

EARB_drugclass_sum %>% arrange(amount) %>% print(n=23)

EARB_drugclass_sum$rel_frequency <- EARB_drugclass_sum$amount / sum(EARB_drugclass_sum$amount) * 100
EARB_drugclass_sum$condition <- c("EARB")

EARB_drugclass <- EARB_drugclass_sum %>% left_join(EARB_drugclass_carrier, by = 'Drug_Class') %>% 
  mutate('average' = amount / carrier) %>% 
  mutate('value' = carrier * rel_frequency / 20)

EARB_drugclass %>% arrange(desc(value)) %>% print(n=50)

## SARB bacteria ##

SARB_drugclass <- total_drugclass %>% inner_join(SARB_list, by=c('Run', 'Bin'))
sum(SARB_drugclass$total_drugclass) # 1444 none_EARB drugclass

#number of carriers#
SARB_drugclass_carrier <- SARB_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n())

SARB_drugclass_carrier %>% print(n=26) # total number of drug classes: 26

#rel_frequency calculate#
SARB_drugclass_sum <- SARB_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount' = sum(total_drugclass))

SARB_drugclass_sum %>% arrange(amount) %>% print(n=26)

SARB_drugclass_sum$rel_frequency <- SARB_drugclass_sum$amount / sum(SARB_drugclass_sum$amount) * 100
SARB_drugclass_sum$condition <- c("SARB")

SARB_drugclass <- SARB_drugclass_sum %>% left_join(SARB_drugclass_carrier, by='Drug_Class') %>% 
  mutate('average' = amount / carrier) %>% 
  mutate('value' = carrier * rel_frequency / 479)

SARB_drugclass %>% arrange(desc(value)) %>% print(n=50)

## merge two results ##
drugclass_compare <- EARB_drugclass %>% select(c(Drug_Class, value)) %>% 
  full_join(SARB_drugclass, by='Drug_Class', suffix = c('.EARB', '.SARB')) %>% 
  select(-c(amount, carrier, condition, average))
drugclass_compare[is.na(drugclass_compare)] <- 0

print(drugclass_compare, n=31) # total drug class: 31

drugclass_compare <- drugclass_compare %>% mutate('gap' = value.EARB - value.SARB) %>% arrange(desc(gap)) %>% print(n=50)
drugclass_compare <- drugclass_compare %>% select(-rel_frequency)

EARB_drugclass_proportion <- data.frame('drugclass' = drugclass_compare %>% select(Drug_Class), 
                                        'value' = drugclass_compare$value.EARB, 
                                        'condition' = 'EARB') %>% arrange(desc(value))

SARB_drugclass_proportion <- data.frame('drugclass' = drugclass_compare %>% select(Drug_Class), 
                                        'value' = drugclass_compare$value.SARB, 
                                        'condition' = 'SARB')

drugclass_proportion <- rbind(EARB_drugclass_proportion, SARB_drugclass_proportion)
rm(EARB_drugclass_proportion, SARB_drugclass_proportion)

custom_colors <- ifelse(drugclass_proportion$Drug_Class[1:31] %in% c("carbapenem","aminoglycoside antibiotic","glycopeptide antibiotic"), 
                        "red", "black")

# Generating fig 2a
ggbarplot(drugclass_proportion, x = "Drug_Class", y = "value", fill = "condition", x.text.angle = 90, facet.by = 'condition',
          position = position_dodge2(preserve = "single"), legend=c('none'), # set dodge position
          palette = c("#00AFBB", "#E7B800", "green"), ylab = 'AMR prevalence value') + # set color palette 
  rremove('xlab') +
  font("xy.text", size = 6) +
  font("legend.title", size = 6) +
  font("legend.text", size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("ylab", size = 6) +
  theme(axis.text.x = element_markdown(colour = ifelse(drugclass_proportion$Drug_Class[1:31] %in% c("carbapenem","aminoglycoside antibiotic","glycopeptide antibiotic"), 
                                                       "red", "black")))

rm(list = ls())

#### (fig2-b,c) EARB Circos plot (M) ####

# EARB and amr gene #

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo

EARB_list <- n_bin_n_amr %>% filter(AMR>16) %>% select(Run, Bin) # number of EARB: 20
#SARB_list <- n_bin_n_amr %>% filter(AMR<16 & AMR>0) %>% select(Run, Bin) # number of SARB: 479

EARB_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin, Best_Hit_ARO, `Drug Class`, `Resistance Mechanism`) %>% inner_join(EARB_list, by=c('Run','Bin')) # 803 in multi-AMR

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

EARB_rgi_result <- EARB_rgi_result %>% left_join(antibiotics_abun_taxo, by=c('Run','Bin'))



## EARB vs AMR gene circos plot ##

genus_list <- (EARB_rgi_result %>% select(genus) %>% unique())$genus # 4 
amr_gene_list <- (EARB_rgi_result %>% select(Best_Hit_ARO) %>% unique())$Best_Hit_ARO # 91

final_result <- data.frame('amr_name'=amr_gene_list)

for (i in 1:length(genus_list)) {
  genus_name <- genus_list[i]
  
  tmp <- data.frame('amr_name'=character(), 'number'=numeric())
  for (k in 1:length(amr_gene_list)) {
    amr_name <- amr_gene_list[k]
    
    tmp_number <- EARB_rgi_result %>% filter(genus == genus_name & Best_Hit_ARO == amr_name) %>% nrow()
    
    tmp2 <- data.frame('amr_name' = amr_name, 'number' = tmp_number)
    
    tmp <- rbind(tmp,tmp2)
  }
  tmp <- tmp %>% rename(!!genus_name := number)
  final_result <- cbind(final_result, tmp %>% select(!!genus_name))
}

rownames(final_result) <- final_result$amr_name
final_result <- final_result[,-1]
final_result_matrix <- as.matrix(final_result)

npg_colors <- c("#009E73", "#D55E00", "#0072B2", "#E69F00")

named_vector <- character(95)

named_vector[1:4] <- npg_colors
nrow(final_result_matrix)
named_vector[5:95] <- 'grey'

names(named_vector)[1:4] <- colnames(final_result_matrix)
names(named_vector)[5:95] <- rownames(final_result_matrix)

circos.clear()


# Find rows where the count of zeros is <= 1
rows_with_less_or_equal_one_zero <- apply(final_result_matrix, 1, function(row) sum(row == 0) <= 1)

# Create a color matrix with default colors (e.g., "white" for no highlight)
color_matrix <- matrix("blue", nrow = nrow(final_result_matrix), ncol = ncol(final_result_matrix))

# Fill in your highlight color for rows with 0 count <= 1
color_matrix[rows_with_less_or_equal_one_zero, ] <- "red"

# Generate the chord diagram
chordDiagram(final_result_matrix, col = color_matrix , grid.col = named_vector, annotationTrack = c('grid'))

circos.clear()

## EARB vs Drug class ##

EARB_drugclass_rgi_result <- EARB_rgi_result %>% select(genus,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class

# Create an empty data frame with the columns
EARB_drugclass <- tibble(genus = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(EARB_drugclass_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(EARB_drugclass_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(genus = EARB_drugclass_rgi_result$genus[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    EARB_drugclass <- bind_rows(EARB_drugclass, new_row)
  }
}

genus_list <- (EARB_rgi_result %>% select(genus) %>% unique())$genus # 4 
drugclass_list <- (EARB_drugclass %>% select(Drug_Class) %>% unique())$Drug_Class # 23

final_result <- data.frame('Drug_Class'=drugclass_list)

for (i in 1:length(genus_list)) {
  genus_name <- genus_list[i]
  
  tmp <- data.frame('Drug_Class'=character(), 'number'=numeric())
  for (k in 1:length(drugclass_list)) {
    drugclass_name <- drugclass_list[k]
    
    tmp_number <- EARB_drugclass %>% filter(genus == genus_name & Drug_Class == drugclass_name) %>% nrow()
    
    tmp2 <- data.frame('Drug_Class' = amr_name, 'number' = tmp_number)
    
    tmp <- rbind(tmp,tmp2)
  }
  tmp <- tmp %>% rename(!!genus_name := number)
  final_result <- cbind(final_result, tmp %>% select(!!genus_name))
}

rownames(final_result) <- final_result$Drug_Class
final_result <- final_result[,-1]
final_result_matrix <- as.matrix(final_result)

npg_colors <- c("#009E73", "#D55E00", "#0072B2", "#E69F00")

named_vector <- character(27)

named_vector[1:4] <- npg_colors
nrow(final_result_matrix) # 23
named_vector[5:27] <- 'grey'

names(named_vector)[1:4] <- colnames(final_result_matrix)
names(named_vector)[5:27] <- rownames(final_result_matrix)

circos.clear()

# Find rows where the count of zeros is <= 1
rows_with_less_or_equal_one_zero <- apply(final_result_matrix, 1, function(row) sum(row == 0) <= 1)

# Create a color matrix with default colors (e.g., "white" for no highlight)
color_matrix <- matrix("blue", nrow = nrow(final_result_matrix), ncol = ncol(final_result_matrix))

# Fill in your highlight color for rows with 0 count <= 1
color_matrix[rows_with_less_or_equal_one_zero, ] <- "red"

# Generate the chord diagram
chordDiagram(final_result_matrix, col = color_matrix , grid.col = named_vector, annotationTrack = c('grid'))

plot(1, 1, type="n", xlab="", ylab="", xlim=c(0, 2), ylim=c(0, 2), xaxt='n', yaxt='n')
legend("center", legend=c(expression(italic("Escherichia")), 
                          expression(italic("Klebsiella")), 
                          expression(italic("Cronobacter")), 
                          expression(italic("Enterobacter"))), 
       fill=npg_colors, cex = 0.8)
legend("bottom", legend = c(expression(genera >= 3), expression(genera < 3)), col = c("red", "blue"), lty = 1, cex = 0.8, inset=c(0, 0.03))
title('EARB drugclass / EARB AMR gene sharing', cex.main=0.6)


circos.clear()

#### (fig2-d,e) SARB Circos plot (M) ####

# SARB and amr gene #

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo

SARB_list <- n_bin_n_amr %>% filter(AMR<16 & AMR>0) %>% select(Run, Bin) # number of SARB: 479

SARB_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin, Best_Hit_ARO, `Drug Class`, `Resistance Mechanism`) %>% inner_join(SARB_list, by=c('Run','Bin')) # 832 in SARB AMR

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

SARB_rgi_result <- SARB_rgi_result %>% left_join(antibiotics_abun_taxo, by=c('Run','Bin'))



## SARB vs AMR gene circos plot ##

phylum_list <- (SARB_rgi_result %>% select(phylum) %>% unique())$phylum # 8
amr_gene_list <- (SARB_rgi_result %>% select(Best_Hit_ARO) %>% unique())$Best_Hit_ARO # 77

final_result <- data.frame('amr_name'=amr_gene_list)

for (i in 1:length(phylum_list)) {
  phylum_name <- phylum_list[i]
  
  tmp <- data.frame('amr_name'=character(), 'number'=numeric())
  for (k in 1:length(amr_gene_list)) {
    amr_name <- amr_gene_list[k]
    
    tmp_number <- SARB_rgi_result %>% filter(phylum == phylum_name & Best_Hit_ARO == amr_name) %>% nrow()
    
    tmp2 <- data.frame('amr_name' = amr_name, 'number' = tmp_number)
    
    tmp <- rbind(tmp,tmp2)
  }
  tmp <- tmp %>% rename(!!phylum_name := number)
  final_result <- cbind(final_result, tmp %>% select(!!phylum_name))
}

rownames(final_result) <- final_result$amr_name
final_result <- final_result[,-1]
final_result_matrix <- as.matrix(final_result)

npg_colors <- c("#009E73", "#D55E00", "#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

named_vector <- character(85)

named_vector[1:8] <- npg_colors
nrow(final_result_matrix)
named_vector[9:85] <- 'grey'

names(named_vector)[1:8] <- colnames(final_result_matrix)
names(named_vector)[9:85] <- rownames(final_result_matrix)

circos.clear()


# Find rows where the count of zeros is <= 2
rows_with_less_or_equal_two_zero <- apply(final_result_matrix, 1, function(row) sum(row == 0) <= 2)

# Create a color matrix with default colors (e.g., "white" for no highlight)
color_matrix <- matrix("blue", nrow = nrow(final_result_matrix), ncol = ncol(final_result_matrix))

# Fill in your highlight color for rows with 0 count <= 1
color_matrix[rows_with_less_or_equal_two_zero, ] <- "red"

# Generate the chord diagram
chordDiagram(final_result_matrix, col = color_matrix , grid.col = named_vector, annotationTrack = c('grid'))

circos.clear()

## SARB vs Drug class ##

SARB_drugclass_rgi_result <- SARB_rgi_result %>% select(phylum,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class

# Create an empty data frame with the columns
SARB_drugclass <- tibble(phylum = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(SARB_drugclass_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(SARB_drugclass_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(phylum = SARB_drugclass_rgi_result$phylum[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    SARB_drugclass <- bind_rows(SARB_drugclass, new_row)
  }
}

phylum_list <- (SARB_rgi_result %>% select(phylum) %>% unique())$phylum # 8
drugclass_list <- (SARB_drugclass %>% select(Drug_Class) %>% unique())$Drug_Class # 26

final_result <- data.frame('Drug_Class'=drugclass_list)

for (i in 1:length(phylum_list)) {
  phylum_name <- phylum_list[i]
  
  tmp <- data.frame('Drug_Class'=character(), 'number'=numeric())
  for (k in 1:length(drugclass_list)) {
    drugclass_name <- drugclass_list[k]
    
    tmp_number <- SARB_drugclass %>% filter(phylum == phylum_name & Drug_Class == drugclass_name) %>% nrow()
    
    tmp2 <- data.frame('Drug_Class' = amr_name, 'number' = tmp_number)
    
    tmp <- rbind(tmp,tmp2)
  }
  tmp <- tmp %>% rename(!!phylum_name := number)
  final_result <- cbind(final_result, tmp %>% select(!!phylum_name))
}

rownames(final_result) <- final_result$Drug_Class
final_result <- final_result[,-1]
final_result_matrix <- as.matrix(final_result)

npg_colors <- c("#009E73", "#D55E00", "#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

named_vector <- character(34)

named_vector[1:8] <- npg_colors
nrow(final_result_matrix)
named_vector[9:34] <- 'grey'

names(named_vector)[1:8] <- colnames(final_result_matrix)
names(named_vector)[9:34] <- rownames(final_result_matrix)

circos.clear()

# Find rows where the count of zeros is <= 2
rows_with_less_or_equal_two_zero <- apply(final_result_matrix, 1, function(row) sum(row == 0) <= 2)

# Create a color matrix with default colors (e.g., "white" for no highlight)
color_matrix <- matrix("blue", nrow = nrow(final_result_matrix), ncol = ncol(final_result_matrix))

# Fill in your highlight color for rows with 0 count <= 1
color_matrix[rows_with_less_or_equal_two_zero, ] <- "red"

# Generate the chord diagram
chordDiagram(final_result_matrix, col = color_matrix , grid.col = named_vector, annotationTrack = c('grid'))

circos.clear()
plot(1, 1, type="n", xlab="", ylab="", xlim=c(0, 2), ylim=c(0, 2), xaxt='n', yaxt='n')
legend("center", legend=gsub(pattern = 'p__',replacement = '',phylum_list), fill=npg_colors, cex = 0.8)
legend("bottom", legend = c(expression(phyla >= 6), expression(phyla < 6)), col = c("red", "blue"), lty = 1, cex = 0.8, inset=c(0, 0.03))

title('SARB drugclass / SARB AMR gene sharing', cex.main=0.6)

circos.clear()

rm(list=ls())

################################################################################

######################################
#####         FIGURE 3           #####
######################################

#### (fig3-b) scoring by portion of kegg ####

setwd("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/fig4_CPA_keggresults/")

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata") # Please change the directory #


EARB_data <- n_bin_n_amr %>% filter(AMR >=17) %>% select(Run) %>% unique()

sample_list=NULL

for (i in 1:nrow(EARB_data)) {
  sample_list[i] <- paste("./",EARB_data$Run[i], sep = "")
}

paste("Total running sample number is", length(sample_list)) # 16


merged_result <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("name","score"))

for (sample in sample_list) {
  
  file_list <- list.files(path = sample, pattern = "keggcount")
  i = 1
  
  for (keggcount in file_list) {
    if (i == 1) {
      keggcount_table <- fread(file = file.path(sample,keggcount))
      keggcount_table <- keggcount_table[,c(2,1)]
      colnames(keggcount_table) <- c("KO",keggcount)
      i =i+1
      next
    }
    tmp <- fread(file = file.path(sample,keggcount))
    colnames(tmp) <- c(keggcount,"KO")
    
    keggcount_table <- keggcount_table %>% full_join(tmp, by = "KO")
    i=i+1
  }
  keggcount_table <- as.data.frame(keggcount_table)
  
  
  rownames(keggcount_table) <- keggcount_table$KO
  keggcount_table <- keggcount_table[,-1]
  keggcount_table[is.na(keggcount_table)] <- 0
  keggcount_table <- keggcount_table/rowSums(keggcount_table)
  
  tmp <- data.frame("name"=colnames(keggcount_table) %>% gsub(pattern = ".keggcount", replacement = ""),"score"=colSums(keggcount_table))
  
  merged_result <- rbind(merged_result,tmp)
  
}

EARB_list <- n_bin_n_amr %>% filter(AMR >=17) %>% select(-AMR)
EARB_list <- EARB_list %>% mutate(multi_amr = "yes")

final_result <- merged_result %>% separate(name, into = c("Run","Bin"), sep = '_')

final_result <- final_result %>% full_join(EARB_list,by=c("Run","Bin"))

final_result$multi_amr <- replace_na(final_result$multi_amr, "no")

EARB_score <- final_result %>% filter(multi_amr == "yes")


final_result <- final_result %>% left_join(metadata, by='Run')
final_result <- final_result %>% group_by(Run,condition) %>% arrange(condition)
final_result$condition <- factor(final_result$condition, levels = c("D4","D8","D42"))
final_result <- final_result %>% group_by(Run,condition) %>% arrange(condition)

ggboxplot(final_result, x="Run",y="score", fill = "Run", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Community power score', palette = c(rep("#E64B35FF",2), rep("#4DBBD5FF", 11), rep("#00A087FF" ,3))) +
  rotate_x_text(angle=45) +
  ggtitle("Community power analysis") +
  rremove('x.text') +
  rremove('xlab') +
  font("ylab", size = 6)+
  font("y.text", size=6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5)) +
  geom_point(data = EARB_score, size=0.6, color='red')

# community power score based on the category #

n_bin_n_amr <- n_bin_n_amr %>% mutate(category=ifelse(AMR == 0, 'NONE',ifelse(AMR < 17, 'SARB','EARB')))

final_result <- final_result %>% left_join(n_bin_n_amr, by=c('Run','Bin'))

final_result$category <- factor(final_result$category, levels = c('NONE','SARB','EARB'))

ggboxplot(final_result, x="category",y="score", fill = "category", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Community power score', palette = c(rep("#91D1C2FF",1), rep("#F39B7FFF", 1), rep("#DC0000FF",1))) +
  rotate_x_text(angle=45)+
  rremove('xlab') +
  font("ylab", size = 6)+
  font("y.text", size=6) +
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm'))+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "NONE",size=3) 

rm(list=ls())

#### (fig3-c) Choose other bacteria for correlation analysis ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/annotated.mOTU.rel_abund.rarefied.Rdata") # load name: annotated_motu_dataframe

cor_two_microbe <- function(data,query,threshold,pvalue,method){
  result <- data.frame(NULL)
  for (i in 1:454){
    
    tmp <- data %>% dplyr::select(c(query, i))
    target_name <- colnames(data)[i]
    if (target_name == query) {
      result <- rbind.data.frame(result,c(colnames(data)[i], 1))
      next
    }
    colnames(tmp) <- c(query, 'target')
    
    if (tmp %>% filter(query != 0 & target != 0) %>% nrow() >= threshold) {
      if (cor.test(data[[query]], data[,i], method=method)$p.value < pvalue) {
        result <- rbind.data.frame(result,c(colnames(data)[i], cor.test(data[[query]], data[,i], method=method, )$estimate))
      }
    } else {
      result <- rbind.data.frame(result,c(colnames(data)[i], 0))
      colnames(result) <- c('targets',query)
    }
  }
  return(result)
} 

compare_specie <- c("Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bilophila wadsworthia",  "Parabacteroides distasonis", "Escherichia coli",
                    "Klebsiella oxytoca","Klebsiella variicola/pneumoniae")

for (file in 1:length(compare_specie)) {
  name <- compare_specie[file]
  
  assign(paste("result",file,sep=""), cor_two_microbe(annotated_motu_dataframe,name,5,1,'spearman'))
  
}

result_list <- list(result1, result2, result3, result4, result5, result6, result7)

result <- Reduce(function(x, y) inner_join(x, y, by = "targets"), result_list)

rownames(result) <- result$targets
result <- result[,-1]

result <- mutate_all(result, function(x) as.numeric(as.character(x)))

result <- result[rowSums(result)!=0,]

filtered_result <- result %>% 
  filter(apply(., 1, function(x) !all(x > -0.3 & x < 0.3)))
rownames(filtered_result) <- gsub('Klebsiella variicola/pneumoniae','Klebsiella pneumoniae', rownames(filtered_result))
colnames(filtered_result) <- gsub('Klebsiella variicola/pneumoniae','Klebsiella pneumoniae', colnames(filtered_result))

filtered_result <- subset(filtered_result, !grepl("motu_linkage_group", rownames(filtered_result)))

newnames <- lapply(
  colnames(filtered_result),
  function(x) bquote(italic(.(x)))
)

pheatmap(filtered_result, show_rownames = F,fontsize = 6, angle_col = 90, treeheight_row = 15, treeheight_col = 10, border_color = NA, labels_col=as.expression(newnames)) 

# Do not remove environments this time.

#### (fig3-d) Positive or Negative related OTU counts ####

# Need previous data, filtered_result #

compare_specie <- gsub('Klebsiella variicola/pneumoniae', 'Klebsiella pneumoniae', compare_specie)

for (file in 1:length(compare_specie)) {
  name <- compare_specie[file]
  print(paste("NUMBER OF POSITIVE CORRELATED SPECIES WITH", name))
  print(filtered_result %>% select(name) %>% filter(.[[name]] > 0.3) %>% nrow())
  
  print(paste("NUMBER OF NEGATIVE CORRELATED SPECIES WITH", name))
  print(filtered_result %>% select(name) %>% filter(.[[name]] < -0.3) %>% nrow())
}

positive <- c(10, 10, 12, 21, 14, 21, 22)
negative <- c(3, 8, 7, 6, 21, 28, 13)

# Combine vectors into matrix using cbind()
mat <- cbind(positive, negative)

# Assign row names to the matrix
rownames(mat) <- c("Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bilophila wadsworthia", 
                   "Parabacteroides distasonis", "Escherichia coli", 
                   "Klebsiella oxytoca", "Klebsiella pneumoniae")

# Convert matrix to data.frame
df <- as.data.frame(mat)
colnames(df) <- c("positive ( > 0.3 )", "negative ( < -0.3 )")
df$specie <- rownames(df)

df$Condition <- c("Non-EARB","Non-EARB","Non-EARB","Non-EARB","EARB","EARB","EARB")

df$Condition <- factor(df$`Condition`, levels = c("Non-EARB", "EARB"))

ggbarplot(data=df, x='specie', y=c("positive ( > 0.3 )", "negative ( < -0.3 )"), combine = T, ylab = 'Positive or negative related OTU count',  fill = "Condition", palette = "npg", legend = 'bottom', size = 0.1) +
  rremove('xlab') +
  rremove('x.text')+
  font("ylab", size = 6)+
  font("y.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6)+
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(strip.text.x = element_text(size = 6)) 

rm(list=ls())

#### (fig3-e,f) correlation between EARB and number of OTU ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/annotated.mOTU.rel_abund.rarefied.Rdata") # load name: annotated_motu_dataframe

annotated_motu_dataframe <- annotated_motu_dataframe %>% select(-c(Run, condition))
annotated_motu_dataframe <- as.data.frame(t(annotated_motu_dataframe)) 
annotated_motu_dataframe$Species <- rownames(annotated_motu_dataframe)

filtered_annotated_motu_dataframe <- subset(annotated_motu_dataframe, !grepl("motu_linkage_group", annotated_motu_dataframe$Species))

filtered_annotated_motu_dataframe <- filtered_annotated_motu_dataframe %>% select(-Species)
rownames(filtered_annotated_motu_dataframe) <- gsub('Klebsiella variicola/pneumoniae','Klebsiella pneumoniae', rownames(filtered_annotated_motu_dataframe))

EARB_species <- c("Escherichia coli", "Klebsiella oxytoca","Klebsiella pneumoniae")
None_EARB <- c("Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bilophila wadsworthia",  "Parabacteroides distasonis")

result <- data.frame(Sample = character(0), `EARB_sum` = numeric(), `None_EARB sum` = numeric(), `OTU_number` = integer())

for (i in 1:ncol(filtered_annotated_motu_dataframe)) {
  name <- colnames(filtered_annotated_motu_dataframe)[i]
  tmp <- filtered_annotated_motu_dataframe %>% mutate(Species=rownames(filtered_annotated_motu_dataframe)) %>% select(Species, !!name)
  
  EARB_sum <- sum(tmp %>% filter(Species %in% EARB_species) %>% select(!!name))
  None_EARB_sum <- sum(tmp %>% filter(Species %in% None_EARB) %>% select(!!name))
  OTU_num <- tmp %>% filter(!!sym(name)>0) %>% nrow()
  
  tmp_result <- data.frame('Sample' = name, 'EARB_sum' = EARB_sum, 'None_EARB_sum' = None_EARB_sum, 'OTU_number' = OTU_num)
  result <- rbind(result, tmp_result)
}

result <- separate(result, Sample, c('Host', 'Condition'))
EARB_result <- result %>% filter(EARB_sum != 0)
EARB_result$EARB_sum <- log2(EARB_result$EARB_sum)

# Generating scatter plot #
ggscatter(EARB_result, x = 'EARB_sum', y = 'OTU_number',
          add = "loess", conf.int = TRUE, cor.coef = FALSE,
          cor.method = 'spearman', size = 1) +
  labs(x = expression(paste("EARB relative abundnace sum (", log[2], ")")),
       y = 'Number of OTUs') +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, 'cm'),
    plot.title = element_text(size = 6, face = "bold"),
    legend.position = "top"
  )

# Box plot
EARB_result <- EARB_result %>% mutate('group'=ifelse(EARB_sum < -10, 'rare', ifelse(EARB_sum >= -5, 'rich', 'normal')))
EARB_result$group <- factor(EARB_result$group, levels = c('rare','normal','rich'))

ggboxplot(EARB_result, x="group",y="OTU_number", fill = "group", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Specie richness', xlab = 'EARB abundance', palette = c(rep("#91D1C2FF",1), rep("#F39B7FFF", 1), rep("#DC0000FF",1))) +
  rotate_x_text(angle=45) +
  ggtitle("EARB abundnace - specie richness") +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "rare",size=3, label.y=68) # label.y full = 180

## None_EARB vs mOTU richness ##
None_EARB_result <- result %>% filter(None_EARB_sum != 0)
None_EARB_result$None_EARB_sum <- log2(None_EARB_result$None_EARB_sum)

# Generating scatter plot #
ggscatter(None_EARB_result, x = 'None_EARB_sum', y = 'OTU_number',
          add = "loess", conf.int = TRUE, cor.coef = FALSE, 
          cor.method = 'spearman', size = 1) +
  labs(x = expression(paste("None EARB relative abundnace sum (", log[2], ")")), 
       y = 'Number of OTUs') +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, 'cm'),
    plot.title = element_text(size = 6, face = "bold"),
    legend.position = "top"
  )


# Box plot #
None_EARB_result <- None_EARB_result %>% mutate('group'=ifelse(None_EARB_sum < -10, 'rare', ifelse(None_EARB_sum >= -5, 'rich', 'normal')))
None_EARB_result$group <- factor(None_EARB_result$group, levels = c('rare','normal','rich'))

ggboxplot(None_EARB_result, x="group",y="OTU_number", fill = "group", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Specie richness', xlab = 'None EARB abundance', palette = c(rep("#91D1C2FF",1), rep("#F39B7FFF", 1), rep("#DC0000FF",1))) +
  rotate_x_text(angle=45) +
  ggtitle("None EARB abundance - specie richness") +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "rare",size=2.5, label.y=71)

rm(list=ls())

################################################################################

######################################
#####         FIGURE 4           #####
######################################

#### (fig4-a) significant 한 pathway 찾아보기 ####

setwd("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/fig5_pathway_module_analysis/")

sample_list <- list.dirs()
sample_list <- sample_list[-1]

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/pathway_list.Rdata") #load name: pathway_list
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/AMR_pathway_SJ.Rdata") #load name: pathway_SJ


EARB_bin <- n_bin_n_amr %>% filter(AMR>=17)

pattern <- EARB_bin %>% select(Run) %>% unique() %>% t() %>% as.vector()
pattern <- paste(pattern, collapse="|")

sample_list <- sample_list[grepl(pattern, sample_list)] # compare only samples which contain EARB_bin

EARB_colname <- paste(EARB_bin$Run,EARB_bin$Bin,sep = '_')

pathway <- fread(file = "./ERR1995223/ERR1995223_merged.tsv") %>% select("ko","pathway")
ko_pathway_table_wo_EARB <- pathway
ko_pathway_table_EARB <- pathway

for (i in sample_list) {
  file_list <- list.files(path = i, pattern = "merged")
  file_list <- file_list[!str_detect(file_list,'module')]
  
  for (k in file_list) {
    merged_file <- fread(file = file.path(i,k))
    
    names.use <- names(merged_file)[(names(merged_file) %in% EARB_colname)]
    EARB <- merged_file %>% select(c(names.use))
    ko_pathway_table_EARB <- cbind(ko_pathway_table_EARB, EARB)
    
    None_EARB <- merged_file %>% select(-c("ko","pathway",names.use))
    ko_pathway_table_wo_EARB <- cbind(ko_pathway_table_wo_EARB, None_EARB)
  }
}

rm(EARB,i,k,file_list,names.use, pattern, sample_list,merged_file, pathway, None_EARB)

merged_pathway <- ko_pathway_table_EARB %>% left_join(ko_pathway_table_wo_EARB, by=c("ko","pathway")) %>% as.data.frame()
rownames(merged_pathway) <- merged_pathway$pathway
merged_pathway <- merged_pathway[,-c(1,2)]

metadata <- data.frame(bin=c(colnames(merged_pathway)), condition=c(rep("EARB",each=20),rep("None_EARB",each=217)))

sig_result = data.frame(Pathway = character(), p.value = numeric(), log2fold = numeric())

for (line in 1:nrow(merged_pathway)) {
  p.value <- wilcox.test(as.numeric(merged_pathway[line,1:20]),as.numeric(merged_pathway[line,21:237]),
                            paired = FALSE,
                            conf.level = 0.95,
                            alternative="two.sided")$p.value %>% as.numeric()
  
  log2fold <- log2(mean(as.numeric(merged_pathway[line,1:20]))/mean(as.numeric(merged_pathway[line,21:237]))) %>% as.numeric()
  
  sig_result <- rbind(sig_result, c(rownames(merged_pathway[line,]),p.value, log2fold))

}


merged_pathway %>%
  select(c(1:20)) %>%
  rowSums()/20 -> merged_pathway$EARB_sum
merged_pathway %>% select(c(21:237)) %>% rowSums()/217 -> merged_pathway$None_EARB_sum

colnames(sig_result) <- c("Pathway","p.value","log2fold")


sig_result$EARB_sum <- merged_pathway$EARB_sum
sig_result$None_EARB_sum <- merged_pathway$None_EARB_sum
sig_result$pathway_length <- lengths(pathway_list)


sig_result$log2fold <- as.numeric(sig_result$log2fold)
sig_result$p.value <- as.numeric(sig_result$p.value)

sig_result <- sig_result %>% filter((EARB_sum/pathway_length > 0.30 | None_EARB_sum/pathway_length > 0.30) & (log2fold > 0.5 | log2fold < -0.5)) %>% arrange(desc(log2fold))
sig_result <- sig_result  %>% filter(Pathway!='Peroxisome')

sig_result <- sig_result %>% left_join(pathway_SJ, by='Pathway')

# Calculate number of bars for each category
category_order <- sig_result %>%
  group_by(Category) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(Category)

# Factor the Category variable based on this order
sig_result$Category <- factor(sig_result$Category, levels = category_order)


ggbarplot(sig_result, x = "Pathway", y = "log2fold",
          fill = "Category",         
          color = "black",           
          palette = "npg",           
          sort.val = "asc",          
          x.text.angle = 90,         
          size = 0.3,
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  labs(x = 'Pathway', 
       y = expression(log["2"]~fold))+
  theme(strip.text.y = element_blank())+
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 5)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") 

rm(list=ls())

#### (fig4-b) significant 한 module 찾아보기 ####

setwd("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/fig5_pathway_module_analysis/")

sample_list <- list.dirs()
sample_list <- sample_list[-1]

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("~/Documents/analysis/EARB_project/final_data_for_figure/module_list.Rdata") #load name: module_list
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/AMR_module_SJ.Rdata") #load name: module_SJ

EARB_bin <- n_bin_n_amr %>% filter(AMR>=17)

pattern <- EARB_bin %>% select(Run) %>% unique() %>% t() %>% as.vector()
pattern <- paste(pattern, collapse="|")

sample_list <- sample_list[grepl(pattern, sample_list)] # compare only samples which contain EARB_bin

EARB_colname <- paste(EARB_bin$Run,EARB_bin$Bin,sep = '_')

module <- fread(file = "./ERR1995223/ERR1995223_module_merged.tsv") %>% rename(Module=Pathway) %>% select("Module")
module_table_wo_EARB <- module
module_table_EARB <- module

for (i in sample_list) {
  file_list <- list.files(path = i, pattern = "module_merged")

  for (k in file_list) {
    merged_file <- fread(file = file.path(i,k))
    
    names.use <- names(merged_file)[(names(merged_file) %in% EARB_colname)]
    EARB <- merged_file %>% select(c(names.use))
    module_table_EARB <- cbind(module_table_EARB, EARB)
    
    None_EARB <- merged_file %>% select(-c("Pathway",names.use))
    module_table_wo_EARB <- cbind(module_table_wo_EARB, None_EARB)
  }
}


rm(EARB,i,k,file_list,names.use, pattern, sample_list,merged_file, module, None_EARB)

merged_module <- cbind(module_table_EARB, module_table_wo_EARB[,-1]) %>% as.data.frame()
metadata <- data.frame(bin=c(colnames(merged_module)[-1]), condition=c(rep("EARB",each=20),rep("None_EARB",each=217)))

sig_result = data.frame(Module = character(), p.value = numeric(), log2fold = numeric())

for (line in 1:nrow(merged_module)) {
  p.value <- wilcox.test(as.numeric(merged_module[line,2:21]),as.numeric(merged_module[line,22:238]),
                            paired = FALSE,
                            conf.level = 0.95,
                            alternative="two.sided")$p.value %>% as.numeric()
  
  log2fold <- log2(mean(as.numeric(merged_module[line,2:21]))/mean(as.numeric(merged_module[line,22:238]))) %>% as.numeric()
  
  sig_result <- rbind(sig_result, c(as.character(merged_module[line,"Module"]),p.value, log2fold))
}

colnames(sig_result) <- c("Module","p.value","log2fold")

merged_module %>%
  select(c(2:21)) %>%
  rowSums()/20 -> merged_module$EARB_sum
merged_module %>% select(c(22:238)) %>% rowSums()/217 -> merged_module$None_EARB_sum

sig_result$EARB_sum <- merged_module$EARB_sum
sig_result$None_EARB_sum <- merged_module$None_EARB_sum
sig_result$module_length <- lengths(module_list)

sig_result$log2fold <- as.numeric(sig_result$log2fold)
sig_result$p.value <- as.numeric(sig_result$p.value)

sig_result <- sig_result %>% filter(EARB_sum/module_length > 0.8 | None_EARB_sum/module_length > 0.8) %>% filter(log2fold > 0.5) %>% arrange(desc(log2fold))


# Calculate number of bars for each category
category_order <- module_SJ %>%
  group_by(Category) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(Category)

# Factor the Category variable based on this order
module_SJ$Category <- factor(module_SJ$Category, levels = category_order)
sig_result <- sig_result %>% left_join(module_SJ, by='Module')

# Calculate number of bars for each category
category_order <- sig_result %>%
  group_by(Category) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(Category)

# Factor the Category variable based on this order
sig_result$Category <- factor(sig_result$Category, levels = category_order)


ggbarplot(sig_result, x = "Module", y = "log2fold",
          fill = "Category",           
          color = "black",            
          palette = paletteer_d("ggsci::category20_d3"),      
          sort.val = "asc",          
          x.text.angle = 90,       
          size = 0.3,
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  labs(x = 'Module', 
       y = expression(log["2"]~fold))+
  theme(strip.text.y = element_blank())+
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 5)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") 


rm(list=ls())




################################################################################

######################################
#####         FIGURE 5           ##### 
######################################

#### (fig5-a) ANI_matrix ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/ANI_matrix.Rdata") #load name: ANI_matrix

colnames(ANI_matrix) <- gsub('_D',' D',colnames(ANI_matrix))
rownames(ANI_matrix) <- gsub('_D',' D',rownames(ANI_matrix))

pheatmap(ANI_matrix,display_numbers = F, fontsize = 6)

### ONLY ECOLI HEATMAP ###

E.coli_list <- ANI_matrix %>% select(`H7 D8.1`) %>% filter(`H7 D8.1`>95) %>% rownames()

ANI_matrix <- ANI_matrix %>% select(E.coli_list) %>% mutate(tmp=rownames(ANI_matrix)) %>% filter(tmp %in% E.coli_list) %>% select(-tmp)

ANI_matrix <- ANI_matrix %>% select(c("H7 D8.1",  "H2 D42",  "H1 D8.1",  "H9 D8",  "H10 D8", "H5 D4",  "H8 D8",  "H12 D4", "H12 D8", "H5 D8",  "H4 D42",  "H6 D8" ))
ANI_matrix <- ANI_matrix[c("H7 D8.1",  "H2 D42",  "H1 D8.1",  "H9 D8",  "H10 D8", "H5 D4",  "H8 D8",  "H12 D4", "H12 D8", "H5 D8",  "H4 D42",  "H6 D8"),]

colnames(ANI_matrix) <- gsub('D8.1','D8',colnames(ANI_matrix)) 
rownames(ANI_matrix) <- gsub('D8.1','D8',rownames(ANI_matrix))

pheatmap(ANI_matrix,display_numbers = F, cluster_cols = F, cluster_rows = F, fontsize = 6, main = "intra-species ANI matrix")

# Do not remove environments this time.

#### (fig5-b) Top 5 ANI score of 23, 24 / 47, 48 ####

ANI_compare <- as.data.frame(ANI_matrix) %>% select(c("H12 D4", "H12 D8","H5 D4","H5 D8"))

result <- as.data.frame(matrix(ncol=3,nrow=0))
for (i in 1:ncol(ANI_compare)) {
  name <- colnames(ANI_compare)[i]

  tmp_df <- ANI_compare %>% select(`name`) %>% top_n(n=5)
  
  tmp <- data.frame('ANI score' = tmp_df[[1]], 'E.coli' = name, 'Close E.coli' = rownames(tmp_df))
  result <- rbind(result,tmp)
}

# Create a named vector of colors
result$Close.E.coli <- as.factor(result$Close.E.coli)
colors <- ifelse(unique(result$Close.E.coli) %in% c("H12 D4", "H12 D8"), "red",
                 ifelse(unique(result$Close.E.coli) %in% c("H5 D4","H5 D8"), "blue", "black"))
names(colors) <- unique(result$Close.E.coli)

ggbarplot(result, x = 'Close.E.coli', y = 'ANI.score', fill = 'E.coli',
          x.text.angle = 45, xlab = 'Top 5 closest E. coli',
          position = position_dodge2(preserve = "single"), 
          legend = "none", 
          palette = get_palette(palette = 'npg',4), 
          ylab = 'ANI score', 
          ylim = c(98, 100)) +
  facet_wrap(~ E.coli, scales = 'free_x', ncol=2) +
  theme(strip.text.x = element_text(size = 6))+
  font("xy.text", size = 6) +
  font("legend.title", size = 6) +
  font("legend.text", size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6)

rm(list=ls())


#### (fig5-d) Top 20 genes-SNP ratio ####

# H5_D8 core #
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D8_core.Rdata")

H5_D8_core <- H5_D8_core %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .)))

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D8_core)) {
  position_depth <- sum(H5_D8_core[[position]])
  nucleotide_freq <- sort(H5_D8_core[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D8_core)[which(H5_D8_core[[position]] == nucleotide_freq[1])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D8_core)[which(H5_D8_core[[position]] == nucleotide_freq[1])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor', 
                              'nucl'=rownames(H5_D8_core)[which(H5_D8_core[[position]] == nucleotide_freq[2])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major 48: 74 / minor 48: 74/ one allele 48: 166

result_48_53  <- result

# H5_D4 core #

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D4_core.Rdata")

H5_D4_core <- H5_D4_core %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .))) # 804

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D4_core)) {
  position_depth <- sum(H5_D4_core[[position]])
  nucleotide_freq <- sort(H5_D4_core[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D4_core)[which(H5_D4_core[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D4_core)[which(H5_D4_core[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor',
                              'nucl'=rownames(H5_D4_core)[which(H5_D4_core[[position]] == nucleotide_freq[2])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major 47: 691 / minor 47: 691/ one allele 48: 113

result_47_53  <- result

# Merge H5 D4 and H5 D8 #

merged_result <- rbind(result_47_53, result_48_53)
head(merged_result)

merged_result$true_condition <- paste(merged_result$status, merged_result$condition)

gghistogram(merged_result, x="SNV_ratio", fill = 'condition', facet.by= 'status', ylab= 'Top 20 genes - Allele frequency', xlab = 'SNV_ratio', pallet = get_palette('npg',6),
            bins = 50) +
  rotate_x_text(angle=45) +
  rremove('xlab') + # If you want to remove the x-label, though it's usually useful
  font("ylab", size = 8)+
  font("xy.text", size = 8)+
  font("legend.title", size = 8)+
  font("legend.text",size = 8) +
  font("title", size = 8, face = "bold") +
  theme(legend.key.size = unit(0.5, 'cm')) 

rm(list=ls())

#### (fig5-e) AMR gene - SNP ratio ####

# H5 D8 AMR # 

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D8_AMR.Rdata") # 1841

H5_D8_AMR <- H5_D8_AMR %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .))) 

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D8_AMR)) {
  position_depth <- sum(H5_D8_AMR[[position]])
  nucleotide_freq <- sort(H5_D8_AMR[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D8_AMR)[which(H5_D8_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D8_AMR)[which(H5_D8_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor', 
                              'nucl'=rownames(H5_D8_AMR)[which(H5_D8_AMR[[position]] == nucleotide_freq[2])], 'status'='H5 D8')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major: 104 / minor: 104/ one allele: 1737

result_48_48amr  <- result

# H5 D4 AMR # 
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D4_AMR.Rdata") # 1841

H5_D4_AMR <- H5_D4_AMR %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .))) 

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D4_AMR)) {
  position_depth <- sum(H5_D4_AMR[[position]])
  nucleotide_freq <- sort(H5_D4_AMR[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor', 
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[2])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major: 1841 / minor: 1841

result_47_48amr  <- result

# H5 D0 AMR # 
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D0_AMR.Rdata") # 1841

H5_D0_AMR <- H5_D0_AMR %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .))) 

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D0_AMR)) {
  position_depth <- sum(H5_D0_AMR[[position]])
  nucleotide_freq <- sort(H5_D0_AMR[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D0_AMR)[which(H5_D0_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D0')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D0_AMR)[which(H5_D0_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D0')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor', 
                              'nucl'=rownames(H5_D0_AMR)[which(H5_D0_AMR[[position]] == nucleotide_freq[2])], 'status'='H5 D0')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major: 75 / minor: 76 / one allele: 333

### median depth ###
one_allele_list <- result %>% filter(condition=='major') %>% pull(position) 

# position_depths 계산
position_depths <- sapply(one_allele_list, function(pos) sum(H5_D0_AMR[[as.character(pos)]]))

# 중앙값 계산
median_depth <- median(position_depths)
print(paste("Median position depth:", median_depth))
# 기본 통계 확인
summary(position_depths)

### ###

# reason of one more minor #
result %>% filter(position=='k113_10916_length_229362_cov_50.0807_113_1197')
H5_D0_AMR %>% select('k113_10916_length_229362_cov_50.0807_113_1197') 

result_44_48amr  <- result

rm(temp_result, result, nucleotide_freq, position, position_depth)

####

result <- result_47_48amr %>% inner_join(result_48_48amr, by=c('position','nucl'), suffix=c('47_48amr','48_48amr'))

temp <- result %>% filter(condition47_48amr=='major') %>% select(c('SNV_ratio48_48amr', 'condition48_48amr'))
temp$condition <- c("Strain 1 (in D8)")
colnames(temp) <- c("SNV_ratio","day4_condition","condition")

temp2 <- result %>% filter(condition47_48amr=='minor') %>% select(c('SNV_ratio48_48amr', 'condition48_48amr'))
temp2$condition <- c("Strain 2 (in D8)")
colnames(temp2) <- c("SNV_ratio","day4_condition","condition")


result <- result_47_48amr %>% select(c(SNV_ratio,condition))
temp <- temp %>% select(c(SNV_ratio,condition))
temp2 <- temp2 %>% select(c(SNV_ratio,condition))

result <- rbind(result, temp, temp2)

result$condition <- gsub('major', 'Strain 1 (in D4)', result$condition)
result$condition <- gsub('minor', 'Strain 2 (in D4)', result$condition)

result48 <- result

#### day 0

result <- result_44_48amr %>% inner_join(result_47_48amr, by=c('position','nucl'), suffix=c('44_48amr','47_48amr'))

result %>% filter(str_detect(condition44_48amr, "major|one allele")) %>% count(condition47_48amr) # 400 major&one allele in day0 becomes minor in day4, 8 vice versa
result %>% filter(str_detect(result$condition44_48amr, 'minor')) %>% count(condition47_48amr) # 61 minor in day0 becomes major in day4, 8 vice versa 

temp <- result %>% filter(str_detect(result$condition44_48amr, 'major')) %>% select(c('SNV_ratio47_48amr', 'condition47_48amr'))
temp$condition <- c("major_in_day0")
colnames(temp) <- c("SNV_ratio","day4_condition","condition")

temp2 <- result %>% filter(str_detect(result$condition44_48amr, 'minor')) %>% select(c('SNV_ratio47_48amr', 'condition47_48amr'))
temp2$condition <- c("minor_in_day0")
colnames(temp2) <- c("SNV_ratio","day4_condition","condition")


result <- result_44_48amr %>% select(c(SNV_ratio,condition))
temp <- temp %>% select(c(SNV_ratio,condition))
temp2 <- temp2 %>% select(c(SNV_ratio,condition))

result <- rbind(result, temp, temp2)

result$condition <- factor(result$condition, levels = c("one allele","major","minor","major_in_day0","minor_in_day0"),
                           labels = c("one allele","Strain 2 (in D0)", "Strain 1 (in D0)","Strain 2 (in D4)","Strain 1(in D4)"))

result04 <- result


result48 # day 4 and 8 result
result04 # day 0 and 4 result
result04 <- result04 %>% filter(condition!='one allele') # exclude D4 one allele position for clarity

result48 <- result48 %>% mutate(day = ifelse(str_detect(condition, 'D4'), 'D4', 'D8'))
result04 <- result04 %>% mutate(day = ifelse(str_detect(condition, 'D4'), 'D4', 'D0'))
result04 <- result04 %>% filter(day!='D4')

result <- rbind(result48, result04)
result$day <- factor(result$day, levels = c('D0','D4','D8'))
result <- result %>% mutate(strain = ifelse(str_detect(condition,'Strain 1'),'Strain 1','Strain 2'))

ggboxplot(result, x="day", y="SNV_ratio", color = "condition", palette =c(rep(c("#E64B35FF","#4DBBD5FF"), each=3)), ylab = 'AMR gene - Allele frequency',
          add = "jitter", size=0.5, add.params = list(size=0.2), legend=c('right'), facet.by = 'strain') +
  rremove('xlab')+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm'))


rm(list=ls())


################################################################################

######################################
#####         FIGURE 6           ##### 
######################################

#### liver metadata ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/liver_metadata.Rdata") #load name: liver_metadata

#### infant metadata ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/infant_metadata.Rdata") #load name: infant_metadata

#### (fig6-a) Drug class distribution for validation ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr

EARB_list <- n_bin_n_amr %>% filter(AMR>16) %>% select(Run, Bin)

amr_gene_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) %>% inner_join(EARB_list, by=c('Run','Bin')) # select drug class

# Create an empty data frame with the columns
amr_gene_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(amr_gene_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(amr_gene_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = amr_gene_rgi_result$Run[i], Bin = amr_gene_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    amr_gene_drugclass <- bind_rows(amr_gene_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(amr_gene_drugclass) # number of row: 2889

# calculate number of drugclass in each MAGs
EARB_drugclass <- amr_gene_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup() # number of row: 1508

sum(EARB_drugclass$total_drugclass) # 2889 

#number of carriers#
EARB_drugclass_carrier <- EARB_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n()) 

EARB_drugclass_carrier %>% print(n=23) # total number of drug classes: 23

#rel_frequency calculate#
EARB_drugclass_sum <- EARB_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount'=sum(total_drugclass))

EARB_drugclass_sum$rel_frequency <- EARB_drugclass_sum$amount / sum(EARB_drugclass_sum$amount) * 100
EARB_drugclass_sum$condition <- c("EARB")


EARB_drugclass <- EARB_drugclass_sum %>% left_join(EARB_drugclass_carrier, by='Drug_Class') %>% mutate('average'=amount/carrier) %>% mutate('value'=carrier*rel_frequency/20)
EARB_drugclass$Drug_Class <- reorder(EARB_drugclass$Drug_Class,EARB_drugclass$rel_frequency)

#### liver information of AMR rgi result 

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/liver_rgi_result.Rdata") #load name: liver_rgi_result

liver_rgi_result <- liver_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class

# Create an empty data frame with the columns
liver_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(liver_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(liver_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = liver_rgi_result$Run[i], Bin = liver_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    liver_drugclass <- bind_rows(liver_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(liver_drugclass) # number of row: 5192

# calculate number of drugclass in each MAGs
liver_drugclass <- liver_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup() # number of row: 592

sum(liver_drugclass$total_drugclass) # 5192

#number of carriers#
liver_drugclass_carrier <- liver_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n()) 

liver_drugclass_carrier %>% print(n=24) # total number of drug classes: 24

#rel_frequency calculate#
liver_drugclass_sum <- liver_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount'=sum(total_drugclass))

liver_drugclass_sum$rel_frequency <- liver_drugclass_sum$amount / sum(liver_drugclass_sum$amount) * 100
liver_drugclass_sum$condition <- c("liver_EARB")


liver_drugclass <- liver_drugclass_sum %>% left_join(liver_drugclass_carrier, by='Drug_Class') %>% mutate('average'=amount/carrier) %>% mutate('value'=carrier*rel_frequency/28)
liver_drugclass$Drug_Class <- reorder(liver_drugclass$Drug_Class, liver_drugclass$rel_frequency)

####

## merge two results ##
liver_healthy_compare <- liver_drugclass %>% select(c(Drug_Class, value)) %>% 
  full_join(EARB_drugclass, by='Drug_Class', suffix = c('.liver','.healthy')) %>% select(-c(amount, carrier, condition, average, rel_frequency))

liver_healthy_compare[is.na(liver_healthy_compare)] <- 0

print(liver_healthy_compare,n=35)

liver_healthy_compare <- liver_healthy_compare %>% mutate('gap'=value.liver - value.healthy) %>% arrange(desc(gap)) 

tmp <- data.frame('drugclass' = liver_healthy_compare %>% select(Drug_Class), 'value' = liver_healthy_compare$value.liver, 'condition' = 'Liver') %>% arrange(desc(value))
tmp2 <- data.frame('drugclass' = liver_healthy_compare %>% select(Drug_Class), 'value' = liver_healthy_compare$value.healthy, 'condition' = 'Healthy')

liver_healthy_proportion <- rbind(tmp, tmp2)
liver_healthy_proportion$Drug_Class <- as.character(liver_healthy_proportion$Drug_Class)

liver_healthy_colors <- ifelse(liver_healthy_proportion$Drug_Class[1:24] %in% c("carbapenem","aminoglycoside antibiotic","glycopeptide antibiotic"), 
                              "red", "black")

ggbarplot(liver_healthy_proportion, x = "Drug_Class", y = "value", fill = "condition", x.text.angle = 90, facet.by = 'condition',
          position = position_dodge2(preserve = "single"), legend=c('none'), # set dodge position
          palette = c("#00AFBB", "#E7B800"), ylab = 'AMR prevalence value') + # set color palette 
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  theme(axis.text.x = element_markdown(colour = liver_healthy_colors))

#### infant information of AMR rgi result 

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/infant_rgi_result.Rdata") #load name: infant_rgi_result

infant_rgi_result <- infant_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class

# Create an empty data frame with the columns
infant_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(infant_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(infant_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = infant_rgi_result$Run[i], Bin = infant_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    infant_drugclass <- bind_rows(infant_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(infant_drugclass) # number of row: 5045

# calculate number of drugclass in each MAGs
infant_drugclass <- infant_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup() # number of row: 788

sum(infant_drugclass$total_drugclass) # 5045

#number of carriers#
infant_drugclass_carrier <- infant_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n()) 

infant_drugclass_carrier %>% print(n=25) # total number of drug classes: 25

#rel_frequency calculate#
infant_drugclass_sum <- infant_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount'=sum(total_drugclass))

infant_drugclass_sum$rel_frequency <- infant_drugclass_sum$amount / sum(infant_drugclass_sum$amount) * 100
infant_drugclass_sum$condition <- c("infant_EARB")

infant_drugclass <- infant_drugclass_sum %>% left_join(infant_drugclass_carrier, by='Drug_Class') %>% mutate('average'=amount/carrier) %>% mutate('value'=carrier*rel_frequency/39)
infant_drugclass$Drug_Class <- reorder(infant_drugclass$Drug_Class, infant_drugclass$rel_frequency)

## merge two results ##
infant_healthy_compare <- infant_drugclass %>% select(c(Drug_Class, value)) %>% 
  full_join(EARB_drugclass, by='Drug_Class', suffix = c('.infant','.healthy')) %>% select(-c(amount, carrier, condition, average, rel_frequency))

infant_healthy_compare[is.na(infant_healthy_compare)] <- 0

print(infant_healthy_compare,n=25)

infant_healthy_compare <- infant_healthy_compare %>% mutate('gap'=value.infant - value.healthy) %>% arrange(desc(gap))

tmp <- data.frame('drugclass' = infant_healthy_compare %>% select(Drug_Class), 'value' = infant_healthy_compare$value.infant, 'condition' = 'infant') %>% arrange(desc(value))
tmp2 <- data.frame('drugclass' = infant_healthy_compare %>% select(Drug_Class), 'value' = infant_healthy_compare$value.healthy, 'condition' = 'Healthy')

infant_healthy_proportion <- rbind(tmp, tmp2)
infant_healthy_proportion$Drug_Class <- as.character(infant_healthy_proportion$Drug_Class)

infant_healthy_colors <- ifelse(infant_healthy_proportion$Drug_Class[1:25] %in% c("carbapenem","aminoglycoside antibiotic","glycopeptide antibiotic"), 
                                "red", "black")

ggbarplot(infant_healthy_proportion, x = "Drug_Class", y = "value", fill = "condition", x.text.angle = 90, facet.by = 'condition',
          position = position_dodge2(preserve = "single"), legend=c('none'), # set dodge position
          palette = c("#00AFBB", "#E7B800"), ylab = 'AMR prevalence value') + # set color palette 
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  theme(axis.text.x = element_markdown(colour = infant_healthy_colors))

#### ecoli information of AMR rgi result 

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/ecoli_rgi_result.Rdata") #load name: ecoli_rgi_result

ecoli_rgi_result <- ecoli_rgi_result %>% select(Run,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class

# Create an empty data frame with the columns
ecoli_drugclass <- tibble(Run = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(ecoli_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(ecoli_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = ecoli_rgi_result$Run[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    ecoli_drugclass <- bind_rows(ecoli_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(ecoli_drugclass) # number of row: 109,304

# calculate number of drugclass in each MAGs
ecoli_drugclass <- ecoli_drugclass %>% group_by(Run, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup() # number of row: 12,210

sum(ecoli_drugclass$total_drugclass) # 109,304

#number of carriers#
ecoli_drugclass_carrier <- ecoli_drugclass %>% select(-total_drugclass) %>% unique() %>% group_by(Drug_Class) %>% summarise(carrier=n()) 

ecoli_drugclass_carrier %>% print(n=28) # total number of drug classes: 28

#rel_frequency calculate#
ecoli_drugclass_sum <- ecoli_drugclass %>% select(Drug_Class, total_drugclass) %>% group_by(Drug_Class) %>% summarise('amount'=sum(total_drugclass))

ecoli_drugclass_sum$rel_frequency <- ecoli_drugclass_sum$amount / sum(ecoli_drugclass_sum$amount) * 100
ecoli_drugclass_sum$condition <- c("ecoli_EARB")


ecoli_drugclass <- ecoli_drugclass_sum %>% left_join(ecoli_drugclass_carrier, by='Drug_Class') %>% mutate('average'=amount/carrier) %>% mutate('value'=carrier*rel_frequency/556)
ecoli_drugclass$Drug_Class <- reorder(ecoli_drugclass$Drug_Class, ecoli_drugclass$rel_frequency)

####

## merge two results ##
ecoli_healthy_compare <- ecoli_drugclass %>% select(c(Drug_Class, value)) %>% 
  full_join(EARB_drugclass, by='Drug_Class', suffix = c('.ecoli','.healthy')) %>% select(-c(amount, carrier, condition, average, rel_frequency))

ecoli_healthy_compare[is.na(ecoli_healthy_compare)] <- 0

print(ecoli_healthy_compare,n=28)

ecoli_healthy_compare <- ecoli_healthy_compare %>% mutate('gap'=value.ecoli - value.healthy) %>% arrange(desc(gap))

tmp <- data.frame('drugclass' = ecoli_healthy_compare %>% select(Drug_Class), 'value' = ecoli_healthy_compare$value.ecoli, 'condition' = 'ecoli') 
tmp2 <- data.frame('drugclass' = ecoli_healthy_compare %>% select(Drug_Class), 'value' = ecoli_healthy_compare$value.healthy, 'condition' = 'Healthy') %>% arrange(desc(value))

ecoli_healthy_proportion <- rbind(tmp2, tmp)
ecoli_healthy_proportion$Drug_Class <- as.character(ecoli_healthy_proportion$Drug_Class)

ecoli_healthy_colors <- ifelse(ecoli_healthy_proportion$Drug_Class[1:28] %in% c("carbapenem","aminoglycoside antibiotic","glycopeptide antibiotic"), 
                               "red", "black")

ggbarplot(ecoli_healthy_proportion, x = "Drug_Class", y = "value", fill = "condition", x.text.angle = 90, facet.by = 'condition',
          position = position_dodge2(preserve = "single"), legend=c('none'), # set dodge position
          palette = c("#00AFBB", "#E7B800"), ylab = 'AMR prevalence value') + # set color palette 
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  theme(axis.text.x = element_markdown(colour = ecoli_healthy_colors))


#### total drug class distribution 

total_drugclass_distribution <- rbind(ecoli_healthy_proportion, liver_healthy_proportion %>% filter(condition!='Healthy'), infant_healthy_proportion %>% filter(condition!='Healthy'))

total_drugclass_distribution$condition <- factor(total_drugclass_distribution$condition, levels=c('Healthy','ecoli','infant','Liver'), labels=c('Healthy','RUTI E.coli','Preterm infant','Liver cirrhosis'))

ggbarplot(total_drugclass_distribution, x = "Drug_Class", y = "value", fill = "condition", x.text.angle = 90, facet.by = 'condition', nrow=1,
          position = position_dodge2(preserve = "single"), legend=c('none'), xlab='Drug class', # set dodge position
          palette = c("#00AFBB", "#E7B800","#33CC33","#FF6347"), ylab = 'AMR prevalence value') + # set color palette 
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  theme(axis.text.x = element_markdown(colour = ecoli_healthy_colors))

rm(list=ls())

#### (fig6-b) CPA for Preterm infant ####

setwd("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/fig6_infant_CPA_result/")

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/infant_metadata.Rdata") #load name: infant_metadata

infant_EARB_data <- infant_metadata %>% filter(number_of_bin > 5 & AMR>16) %>% select(Run) %>% unique() # 24

sample_list=NULL

for (i in 1:nrow(infant_EARB_data)) {
  sample_list[i] <- paste("./", infant_EARB_data$Run[i], sep = "")
}

paste("Total running sample number is", length(sample_list)) # 24

merged_result <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("name","score"))

for (sample in sample_list) {
  
  file_list <- list.files(path = sample, pattern = "keggcount")
  i = 1
  
  for (keggcount in file_list) {
    if (i == 1) {
      keggcount_table <- fread(file = file.path(sample,keggcount))
      keggcount_table <- keggcount_table[,c(2,1)]
      colnames(keggcount_table) <- c("KO",keggcount)
      i =i+1
      next
    }
    tmp <- fread(file = file.path(sample,keggcount))
    colnames(tmp) <- c(keggcount,"KO")
    
    keggcount_table <- keggcount_table %>% full_join(tmp, by = "KO")
    i=i+1
  }
  keggcount_table <- as.data.frame(keggcount_table)
  
  
  rownames(keggcount_table) <- keggcount_table$KO
  keggcount_table <- keggcount_table[,-1]
  keggcount_table[is.na(keggcount_table)] <- 0
  keggcount_table <- keggcount_table/rowSums(keggcount_table)
  
  tmp <- data.frame("name"=colnames(keggcount_table) %>% gsub(pattern = ".keggcount", replacement = ""),"score"=colSums(keggcount_table))
  
  merged_result <- rbind(merged_result,tmp)
  
}

infant_EARB_list <- infant_metadata %>% filter(number_of_bin > 5 & AMR > 16) %>% select(Run, Bin)
infant_EARB_list <- infant_EARB_list %>% mutate(multi_amr = "yes")

final_result <- merged_result %>% separate(name, into = c("Run","Bin"), sep = '_')

final_result <- final_result %>% full_join(infant_EARB_list, by=c("Run","Bin"))

final_result$multi_amr <- replace_na(final_result$multi_amr, "no")

infant_EARB_score <- final_result %>% filter(multi_amr == "yes")

infant_metadata <- infant_metadata %>% select(Run, host_day_of_life, Bin)

day_of_life <- c(infant_metadata %>% select(host_day_of_life) %>% arrange(host_day_of_life) %>% unique() %>% c())$host_day_of_life

final_result <- final_result %>% left_join(infant_metadata,by=c('Run','Bin'))
final_result$host_day_of_life <- factor(final_result$host_day_of_life, levels = day_of_life)
final_result <- final_result %>% arrange(host_day_of_life)

ggboxplot(final_result, x="Run",y="score", fill = "Run", legend = 'none', size = 0.1, outlier.size = 0.6, ylab = 'Community power score') +
  rotate_x_text(angle=45) +
  ggtitle("Community power analysis (Preterm infant)") +
  rremove('xlab')+
  rremove('x.text')+
  font("ylab", size = 6)+
  # font("xy.text", size = 6)+
  font('y.text',size=6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  geom_point(data = infant_EARB_score, size=0.6, color='red')


# ggboxplot(final_result, x="multi_amr",y="score", fill = "multi_amr", ylab = 'Community power score', legend='none', ylim =c(0,1500), size = 0.1, outlier.size=0.6, xlab = 'Multi AMR') +
#   stat_compare_means(label= "p.signif", method = "t.test", ref.group = "no", comparisons = list(c("yes","no")), label.y = 1350, vjust = 0, size=3,) +
#   # stat_compare_means(label = "p.signif", method = "t.test") +
#   font("xlab", size = 6)+
#   font("ylab", size = 6)+
#   font("xy.text", size = 6)+
#   font("legend.title", size = 6)+
#   font("legend.text",size = 6) +
#   font("title", size = 6, face = "bold") +
#   theme(legend.key.size = unit(0.3, 'cm'))

rm(list=ls())

#### (fig6-b) CPA for liver ####

setwd("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/fig6_liver_CPA_result/")

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/liver_metadata.Rdata") #load name: liver_metadata

liver_EARB_data <- liver_metadata %>% filter(number_of_bin > 5 & AMR>16) %>% select(Run) %>% unique() # 28

sample_list=NULL

for (i in 1:nrow(liver_EARB_data)) {
  sample_list[i] <- paste("./", liver_EARB_data$Run[i], sep = "")
}

paste("Total running sample number is", length(sample_list)) # 28

merged_result <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("name","score"))

for (sample in sample_list) {
  
  file_list <- list.files(path = sample, pattern = "keggcount")
  i = 1
  
  for (keggcount in file_list) {
    if (i == 1) {
      keggcount_table <- fread(file = file.path(sample,keggcount))
      keggcount_table <- keggcount_table[,c(2,1)]
      colnames(keggcount_table) <- c("KO",keggcount)
      i =i+1
      next
    }
    tmp <- fread(file = file.path(sample,keggcount))
    colnames(tmp) <- c(keggcount,"KO")
    
    keggcount_table <- keggcount_table %>% full_join(tmp, by = "KO")
    i=i+1
  }
  keggcount_table <- as.data.frame(keggcount_table)
  
  
  rownames(keggcount_table) <- keggcount_table$KO
  keggcount_table <- keggcount_table[,-1]
  keggcount_table[is.na(keggcount_table)] <- 0
  keggcount_table <- keggcount_table/rowSums(keggcount_table)
  
  tmp <- data.frame("name"=colnames(keggcount_table) %>% gsub(pattern = ".keggcount", replacement = ""),"score"=colSums(keggcount_table))
  
  merged_result <- rbind(merged_result,tmp)
  
}

liver_EARB_list <- liver_metadata %>% filter(number_of_bin > 5 & AMR > 16) %>% select(Run, Bin)
liver_EARB_list <- liver_EARB_list %>% mutate(multi_amr = "yes")

final_result <- merged_result %>% separate(name, into = c("Run","Bin"), sep = '_')

final_result <- final_result %>% full_join(liver_EARB_list, by=c("Run","Bin"))

final_result$multi_amr <- replace_na(final_result$multi_amr, "no")

liver_EARB_score <- final_result %>% filter(multi_amr == "yes")

liver_metadata <- liver_metadata %>% select(Run, MELD, Bin)

MELD <- c(liver_metadata %>% select(MELD) %>% arrange(MELD) %>% unique() %>% c())$MELD

final_result <- final_result %>% left_join(liver_metadata,by=c('Run','Bin'))
final_result$MELD <- factor(final_result$MELD, levels = MELD)
final_result <- final_result %>% arrange(MELD)

ggboxplot(final_result, x="Run",y="score", fill = "Run", legend = 'none', size = 0.1, outlier.size = 0.6, ylab = 'Community power score') +
  rotate_x_text(angle=45) +
  ggtitle("Community power analysis (Liver cirrhosis)") +
  rremove('xlab') +
  rremove('x.text') +
  font("ylab", size = 6)+
  # font("xy.text", size = 6)+
  font('y.text', size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  geom_point(data = liver_EARB_score, size=0.6, color='red')

# ggboxplot(final_result, x="multi_amr",y="score", fill = "multi_amr", ylab = 'community power score', legend='none', ylim=c(0,1500), size = 0.1, outlier.size=0.6, xlab='Multi AMR') +
#   stat_compare_means(label= "p.signif", method = "t.test", ref.group = "no", comparisons = list(c("yes","no")), label.y = 1300, vjust = -0.5, size=3) +
#   # stat_compare_means(label = "p.signif", method = "t.test", label.y = 1100) +
#   font("xlab", size = 6)+
#   font("ylab", size = 6)+
#   font("xy.text", size = 6)+
#   font("legend.title", size = 6)+
#   font("legend.text",size = 6) +
#   font("title", size = 6, face = "bold") +
#   theme(legend.key.size = unit(0.3, 'cm'))

rm(list=ls())

#### (fig6-c) pie chart (Hard coding) ####

piechart <- data.frame(number_of_samples = c(303,96,60,532,7,321,53,211), Type = c(rep(c('EARB','None EARB'),4)), condition = rep(c('Infant (n = 399)','Liver (n = 592)','Liver saliva (n = 328)','Liver stool (n = 264)'),each=2))

tmp <- data.frame(number_of_samples = c(556), Type = c('EARB'), condition = c('RUTI E. coli (n = 556)'))

piechart <- rbind(tmp,piechart)
piechart$condition <- factor(piechart$condition, levels= c('RUTI E. coli (n = 556)','Infant (n = 399)','Liver (n = 592)','Liver saliva (n = 328)','Liver stool (n = 264)'))

piechart <- piechart %>%
  group_by(condition) %>%
  mutate(relative_abundance = number_of_samples / sum(number_of_samples))

my_labeller <- function(variable_value){
  gsub(" \\(", "\n(", variable_value)
}

# ggplot으로 파이 차트 생성
ggplot(piechart, aes(x="", y=relative_abundance, fill=Type)) + 
  scale_fill_npg() +
  geom_bar(width = 1, stat = "identity", color="white") +
  geom_text(aes(label=paste0(round(relative_abundance*100, 1), '%')), 
            position = position_stack(vjust = 0.5), size = 2.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ condition, labeller = as_labeller(my_labeller), nrow=1)+ 
  theme_void() +
  theme(legend.position = "bottom") +
  labs(fill = "Multi_AMR") +
  guides(fill = guide_legend(title = "Condition")) +
  theme(
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    plot.title = element_text(size = 6, face = "bold"),
    legend.key.size = unit(0.3, 'cm')
  )

rm(list=ls())

# 어도비에서 E.coli 바꿔!



######################################
#####       Supple Figure        ##### 
######################################


######################################
#####         Extend 1           #####
######################################

#### (ext fig1-b) taxo analysis by condition and abundance ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% select("Run","Bin","abundance","multi_amr", "phylum")

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

# calculate relative abundnace at pylum level in each sample
all_abundance <- antibiotics_abun_taxo %>% group_by(Run, phylum) %>% summarise(count=sum(abundance)) %>% mutate(relative_abundance=count/sum(count)*100) 

# summing relative abundance at condition level then recalculate relative abundance in condition level
all_abundance <- all_abundance %>% left_join(metadata, by='Run') %>% ungroup() %>% group_by(condition, phylum) %>% summarise(count=sum(relative_abundance)) 

all_abundance <- all_abundance %>% group_by(condition) %>% mutate(relative_abundance=count/sum(count)*100) %>% 
  mutate(top=ifelse(count %in% tail(sort(count),5), 'top5','others')) # only consider the top 5 abundant phyla

all_abundance <- all_abundance %>% mutate(phylum=ifelse(top=="others","others",phylum))

all_abundance$phylum <- gsub('p__','',all_abundance$phylum)

# Generating figure 2b
ggbarplot(all_abundance, x = "condition", y = "relative_abundance",
          fill = 'phylum', color = "phylum", legend = c('right'), ylab = 'Relative abundance',
          palette = get_palette('npg',11)) +
  rremove('xlab') +
  # font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5))

rm(list=ls())

#### (ext fig1-c) Enterobacteriaceae vs Veillonellaceae ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

antibiotics_abun_taxo <- antibiotics_abun_taxo <- antibiotics_abun_taxo %>% select("Run","Bin","abundance","multi_amr", "family")

# calculate relative abundnace at family level in each sample
all_abundance <- antibiotics_abun_taxo %>% group_by(Run, family) %>% summarise(count=sum(abundance)) %>% mutate(relative_abundance=count/sum(count)*100)

# summing relative abundance at condition level then recalculate relative abundance in condition level
all_abundance <- all_abundance %>% left_join(metadata, by='Run') %>% ungroup() %>% group_by(condition, family) %>% summarise(count=sum(relative_abundance))
all_abundance <- all_abundance %>% group_by(condition) %>% mutate(relative_abundance=count/sum(count)*100) %>% 
  mutate(top=ifelse(count %in% tail(sort(count),5), 'top5','others'))

# select Veillonellaceae and Enterobacteriaceae and make others to 'others'
all_abundance <- all_abundance %>% mutate(family=ifelse(family=="f__Veillonellaceae" | family=="f__Enterobacteriaceae",family,'others')) %>% mutate(family=ifelse(is.na(family), 'others',family))


all_abundance$family <- gsub('f__','',all_abundance$family)
all_abundance$family <- factor(all_abundance$family, levels = c('Enterobacteriaceae','Veillonellaceae','others'))

# Generating figure 2c
ggbarplot(all_abundance, x = "condition", y = "relative_abundance",
          fill = 'family', color = "family", legend = "top", ylab = 'Relative abundance',
          palette = get_palette('npg', 3)) +
  rremove('xlab') +
  font("ylab", size = 6) +
  font("xy.text", size = 6) +
  font("legend.title", size = 6) +
  font("legend.text", size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.text = element_text(face = 'italic'),  # Apply italic to all legend text
        legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5))

rm(list=ls())

#### (ext fig1-d) multi_amr top5 percentage ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

antibiotics_abunance <- antibiotics_abun_taxo %>% select("Run","Bin","abundance","multi_amr")

EARB_count <- antibiotics_abunance %>% group_by(Run) %>% filter(multi_amr=='EARB') %>% summarise(n=n())
sample_list <- antibiotics_abunance$Run %>% unique()

abundance_EARB_sample_top5 <- data.frame(NULL)

# compare relative abundnace between EARB and top5 microbes
for (name in sample_list) {
  current_data <- antibiotics_abunance %>% filter(Run == name)
  
  # Select top5 non-EARB bacteria
  top5_indices <- current_data %>%
    filter(multi_amr != "EARB") %>% 
    arrange(desc(abundance)) %>%
    slice(1:5) %>%
    pull(Bin)
  
  # calculate relative_abundance and categorize condition
  tmp <- current_data %>%
    mutate(relative_abundance = abundance / sum(abundance) * 100) %>%
    mutate(condition = ifelse(multi_amr == 'EARB', 'EARB',
                              ifelse(Bin %in% top5_indices, 'top5', 'others')))
  
  abundance_EARB_sample_top5 <- rbind(abundance_EARB_sample_top5, tmp)
}

# rename for clarity
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% left_join(metadata, by='Run')
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% rename(top5=condition.x, Day=condition.y)

abundance_EARB_sample_top5$top5 <- gsub('top5', 'TOP 5', abundance_EARB_sample_top5$top5)
abundance_EARB_sample_top5$top5 <- gsub('others', 'others', abundance_EARB_sample_top5$top5)

abundance_EARB_sample_top5$Day <- factor(abundance_EARB_sample_top5$Day, levels = c("D0", "D4", "D8", "D42", "D180"))
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% group_by(Run,Day) %>% arrange(desc(Day))

# select EARB existed sample only (If you skip this, you can get a graph with all 57 samples)
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% filter(Run %in% c(EARB_count$Run))

# Generating figure 2d
ggbarplot(abundance_EARB_sample_top5, x = "Run", y = "relative_abundance",
          fill = 'top5', color = "top5", legend = c('right'),
          palette = get_palette('npg',3), ylab = 'Relative abundance', legend.title="condition") +
  rremove('y.text') +
  rremove('ylab') +
  font("xlab", size = 6)+
  font('x.text', size =6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5)) +
  rotate()

rm(list=ls())




################################################################################


#### (Sub fig1-a) Good and bad percentage of three binning tools (Hard coding) ####

piechart <- data.frame(number_of_samples = c(1769,850,1799,682,1700,1053), quality = c(rep(c('LQ','HQ'),3)), condition = rep(c('Concoct (n=2619)','Maxbin2 (n=2481)','Metabat2 (n=2753)'),each=2))

piechart <- piechart %>%
  group_by(condition) %>%
  mutate(relative_abundance = number_of_samples / sum(number_of_samples))

ggplot(piechart, aes(x="", y=relative_abundance, fill=quality)) + 
  scale_fill_npg() +
  geom_bar(width = 1, stat = "identity", color="white") +
  geom_text(aes(label=paste0(round(relative_abundance*100,1), '%')),
            position=position_stack(vjust=0.5),size=2.5) +
  coord_polar("y", start=0) + 
  facet_wrap(~condition) + 
  theme_void() +
  theme(legend.position = "bottom") +
  labs(fill = "Multi_AMR") +
  guides(fill=guide_legend(title="MAG quality")) +
  theme(
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    plot.title = element_text(size = 6, face = "bold"),
    legend.key.size = unit(0.3, 'cm')
  )

rm(list=ls())


#### (Sub fig1-b) mag quality check by day ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_MAGs_quality.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

all_mags_stats <- total_MAGs_stats %>% left_join(metadata, by="Run")

all_mags_stats <- all_mags_stats %>% mutate("quality" = ifelse(completeness > 70 & contamination < 5, 'HQ', 'LQ'))

all_mags_stats <- all_mags_stats %>% group_by(quality,condition) %>% summarise(number=n())

all_mags_stats <- all_mags_stats %>% group_by(condition) %>% mutate(percentage=number/sum(number)*100) #percentage

all_mags_stats$condition <- factor(all_mags_stats$condition, levels = c("D0", "D4", "D8", "D42", "D180"))

ggbarplot(all_mags_stats, x = "condition", y = "percentage",
          fill = 'quality', color = "quality",
          palette = get_palette('npg',2), x.text.angle = 45, ylab='Percentage', legend.title='Quality (Q > 70, T < 5)') +
  scale_y_continuous(expand=c(0,0)) +
  rremove('xlab') +
  # font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5)) 

rm(list=ls())


#### (Sub fig1-c) Multi-AMR bacteria ratio heatmap ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

amr_per_bin <- n_bin_n_amr %>% left_join(metadata, by='Run')

amr_per_bin$condition <- factor(amr_per_bin$condition, levels = c("D0","D4","D8","D42","D180"))

multi_amr_ratio <- amr_per_bin %>% group_by(Run,sample, condition) %>% summarize(bin=n()) %>% left_join(amr_per_bin %>% filter(AMR>=17) %>% group_by(Run, condition) %>% summarize(EARB=n()), 
                                                                                                        by=c('Run','condition')) %>% ungroup()

multi_amr_ratio$EARB <-  replace_na(multi_amr_ratio$EARB, 0)
multi_amr_ratio$sample <- gsub('ERAS','H',multi_amr_ratio$sample)

multi_amr_ratio$multi_amr_ratio <- multi_amr_ratio$EARB/multi_amr_ratio$bin


multi_amr_ratio <- multi_amr_ratio %>%
  group_by(sample, condition) %>%
  summarise(multi_amr_ratio = mean(multi_amr_ratio, na.rm = TRUE)) %>%
  spread(key = condition, value = multi_amr_ratio)

# Convert tibble to data frame
wide_data <- as.data.frame(multi_amr_ratio)

# Optionally, set row names and remove the sample column
rownames(wide_data) <- wide_data$sample
wide_data$sample <- NULL  # remove the sample column
wide_data[is.na(wide_data)] <- 0
print(wide_data)
annotation_row <- data.frame('Susceptiblity'=c(rep('susceptible',2), rep('tolerant',2), rep('medium',8)))
rownames(annotation_row) <- c('H12','H5','H2','H4','H3','H7','H10','H6','H11','H8','H1','H9')

pheatmap(as.matrix(wide_data),cluster_cols = F, fontsize = 6, treeheight_row = 20,  color = colorRampPalette(c('white','red'))(50),cutree_rows = 3
         , annotation_row = annotation_row,annotation_names_row = F)

rm(list=ls())

#### (Sub fig1-d) PCoA plot of mOTU data ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/annotated.mOTU.rel_abund.rarefied.Rdata")

annotated_motu_dataframe$condition <- factor(annotated_motu_dataframe$condition, levels = c("D0", "D4", "D8", "D42", "D180"))
annotated_motu_table <- annotated_motu_dataframe %>% select(-c(Run,condition)) %>% t() %>% as.data.frame()

tmp_metadata <- data.frame('name'=rownames(annotated_motu_dataframe), 'Run'=annotated_motu_dataframe$Run, 'Condition'=annotated_motu_dataframe$condition)

dat.dist <- vegdist(t(annotated_motu_table),method = "bray")

dat.pcoa=cmdscale(dat.dist)

colnames(dat.pcoa) <- c("X","Y")

dat.binded <- cbind(dat.pcoa,tmp_metadata)

tmp <- dat.binded %>% filter(Condition %in% c('D0', 'D180'))
tmp$Run <- factor(tmp$Run, levels = paste0("ERAS", 1:12), labels = paste0("H", 1:12))

distances <- tmp %>%
  group_by(Run) %>%
  summarize(
    Distance = sqrt(sum((X[Condition == "D180"] - X[Condition == "D0"])^2 +
                          (Y[Condition == "D180"] - Y[Condition == "D0"])^2))
  ) %>% arrange(desc(Distance))

# 결과 출력
print(distances)

connections <- tmp %>% 
  filter(Run %in% c("H12", "H5") & (grepl("Dag0", name) | grepl("Dag180", name))) %>%
  group_by(Run) %>%
  summarize(start_x = first(X), start_y = first(Y),
            end_x = last(X), end_y = last(Y))

# Create the scatter plot with lines or arrows
ggplot(tmp, aes(x = X, y = Y, colour = Run)) +
  geom_point() +
  geom_segment(data = connections, aes(x = start_x, y = start_y, 
                                       xend = end_x, yend = end_y), 
               arrow = arrow(type = "closed", length = unit(0.2, "cm")),
               size = 0.5) +
  scale_colour_paletteer_d(palette = "ggsci::category20_d3") +
  theme_bw() +
  labs(x = "X",
       y = "Y",
       shape = "Time point",  # setting legend title for shape
       color = "Host") +      # setting legend title for color
  theme(
    text = element_text(size = 6),
    legend.key.size = unit(0.3, 'cm'),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold")
  ) +
  guides(
    color = guide_legend(title.position = "top"),
    shape = guide_legend(title.position = "top")
  )

rm(list=ls())

#### (Sub fig2-a) Empirical Cumulative Density Function (N) ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata")

n_bin_n_amr <- n_bin_n_amr %>% filter(AMR != 0)

ecdf_func <- ecdf(n_bin_n_amr$AMR)

AMR_amount <- n_bin_n_amr %>% select(AMR)

ggplot(AMR_amount, aes(AMR)) + stat_ecdf(geom = "point")+
  labs(title="Empirical Cumulative \n Density Function",
       y = "Density", x="Number of AMR")+
  geom_vline(aes(xintercept = 17), color = "blue", linetype = "dashed")+
  geom_hline(aes(yintercept = 0.96), color = "blue", linetype = "dashed")+
  theme_classic() +
  font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 5)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") 

rm(list=ls())

#### (Sub fig2-b) quality of EARB and Non-EARB ####

load('/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/good_MAGs_stats.Rdata')
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata")

good_MAGs_stats <- good_MAGs_stats %>% left_join(n_bin_n_amr, by=c('Run','Bin')) %>% 
  mutate('condition'=ifelse(AMR>16,'EARB',ifelse(AMR>0,'Non-EARB','None carrier'))) %>%
  filter(condition %in% c('EARB','Non-EARB')) 
good_MAGs_stats$condition <- factor(good_MAGs_stats$condition, levels = c('EARB','Non-EARB'))

ggboxplot(good_MAGs_stats, x="condition",y="completeness", fill = "condition", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Completeness') +
  rremove("xlab")+
  font("xy.text", size = 6)+
  font("ylab", size = 6)

ggboxplot(good_MAGs_stats, x="condition",y="contamination", fill = "condition", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Contamination') +
  rremove("xlab")+
  font("xy.text", size = 6)+
  font("ylab", size = 6)

rm(list=ls())


#### (Sub fig3-a,b) SARB, None EARB relative abundance stacked bar ####
# None multi-AMR ( 0 < amr < 17 ) #

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata")

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% select("Run","Bin","abundance","multi_amr", "phylum")

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

# calculate relative abundnace at pylum level in each sample
all_abundance <- antibiotics_abun_taxo %>% group_by(Run, Bin) %>% summarise(count=sum(abundance)) %>% mutate(relative_abundance=count/sum(count)*100)
all_abundance <- all_abundance %>% left_join(metadata, by="Run")

all_abundance <- all_abundance %>% left_join(n_bin_n_amr, by=c('Run','Bin'))
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% select(Run, Bin, phylum)

all_abundance <- all_abundance %>% left_join(antibiotics_abun_taxo,by=c("Run","Bin"))

all_abundance <- all_abundance %>% mutate('status'=ifelse(AMR==0, 'None AMR',ifelse(AMR<17, 'SARB','EARB')))

table(all_abundance$status)

all_abundance <- all_abundance %>% group_by(status, phylum) %>% summarise(count=sum(relative_abundance)) %>% mutate(relative_abundance=count/sum(count)*100)

all_abundance$phylum <- ifelse(is.na(all_abundance$phylum), 'unclassified', all_abundance$phylum)

SARB_abundance <- all_abundance %>% filter(status=='SARB')
SARB_abundance$phylum <- gsub('p__','',SARB_abundance$phylum)
SARB_abundance$phylum <- gsub('_I','',SARB_abundance$phylum)

ggbarplot(SARB_abundance, x = "status", y = "relative_abundance",
          fill = 'phylum', color = "phylum", legend = c('right'), palette = get_palette('npg', 10), 
          ylab = 'Relative abundance ( % )', legend.title="Phylum", ylim = c(0,100)) +
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("ylab", size = 6) +
  rremove("xlab")

NA_abundance <- all_abundance %>% filter(status=='None AMR')

NA_abundance$phylum <- gsub('p__','',NA_abundance$phylum)
NA_abundance$phylum <- gsub('_I','',NA_abundance$phylum)
NA_abundance$status <- 'None EARB'

ggbarplot(NA_abundance, x = "status", y = "relative_abundance",
          fill = 'phylum', color = "phylum", legend = c('right'), palette = get_palette('npg', 10),
          ylab = 'Relative abundance ( % )', legend.title="Phylum", ylim = c(0,100)) +
  font("xy.text", size = 6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("ylab", size = 6) +
  rremove("xlab")

rm(list=ls())

#### (Sub fig3-c) taxo analysis using kraken2 ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

taxo_kraken <- fread("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotic_phylum_kraken.tsv", header = F)
colnames(taxo_kraken) <- c('Run','phylum','abundance')

taxo_kraken$Run <- gsub('_mpa_pylum.mpa','', taxo_kraken$Run)
taxo_kraken$phylum <- gsub("d__Bacteria\\|",'', taxo_kraken$phylum)

taxo_kraken <- taxo_kraken %>% group_by(Run, phylum) %>% mutate(abundance=sum(abundance)) %>% unique()

# calculate relative abundnace at pylum level in each sample
all_abundance <- taxo_kraken %>% group_by(Run, phylum) %>% summarise(count=sum(abundance)) %>% mutate(relative_abundance=count/sum(count)*100) 

# summing relative abundance at condition level then recalculate relative abundance in condition level
all_abundance <- all_abundance %>% left_join(metadata, by='Run') %>% ungroup() %>% group_by(condition, phylum) %>% summarise(count=sum(relative_abundance)) 

all_abundance <- all_abundance %>% group_by(condition) %>% mutate(relative_abundance=count/sum(count)*100) %>% 
  mutate(top=ifelse(count %in% tail(sort(count),5), 'top5','others')) # only consider the top 5 abundant phyla

all_abundance <- all_abundance %>% mutate(phylum=ifelse(top=="others","others",phylum))

all_abundance$phylum <- gsub('p__','',all_abundance$phylum)

# Generating figure 2b
ggbarplot(all_abundance, x = "condition", y = "relative_abundance",
          fill = 'phylum', color = "phylum", legend = c('right'), ylab = 'Relative abundance',
          palette = get_palette('npg',11)) +
  rremove('xlab') +
  # font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5))

rm(list=ls())
#### (Sub fig3-d) top5 day 별로 ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") # load name: antibiotics_abun_taxo
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_metadata.Rdata")

antibiotics_abunance <- antibiotics_abun_taxo %>% select("Run","Bin","abundance","multi_amr")

sample_list <- antibiotics_abunance$Run %>% unique()

abundance_EARB_sample_top5 <- data.frame(NULL)

# compare relative abundnace between EARB and top5 microbes
for (name in sample_list) {
  current_data <- antibiotics_abunance %>% filter(Run == name)
  
  # Select top5 non-EARB bacteria
  top5_indices <- current_data %>%
    filter(multi_amr != "EARB") %>% 
    arrange(desc(abundance)) %>%
    slice(1:5) %>%
    pull(Bin)
  
  # calculate relative_abundance and categorize condition
  tmp <- current_data %>%
    mutate(relative_abundance = abundance / sum(abundance) * 100) %>%
    mutate(condition = ifelse(multi_amr == 'EARB', 'EARB',
                              ifelse(Bin %in% top5_indices, 'top5', 'others')))
  
  abundance_EARB_sample_top5 <- rbind(abundance_EARB_sample_top5, tmp)
}

# rename for clarity
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% left_join(metadata, by='Run')
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% rename(top5=condition.x, Day=condition.y)

abundance_EARB_sample_top5$top5 <- gsub('top5', 'TOP 5', abundance_EARB_sample_top5$top5)
abundance_EARB_sample_top5$top5 <- gsub('others', 'others', abundance_EARB_sample_top5$top5)

abundance_EARB_sample_top5$Day <- factor(abundance_EARB_sample_top5$Day, levels = c("D0", "D4", "D8", "D42", "D180"))
abundance_EARB_sample_top5 <- abundance_EARB_sample_top5 %>% group_by(Day) %>% arrange(Day)

test <- abundance_EARB_sample_top5 %>% group_by(Run, top5) %>% mutate(relative_abundnace_sum=sum(relative_abundance)) %>% select(Run, top5,relative_abundnace_sum, Day) %>% unique()

test <- test %>% ungroup() %>% group_by(Day, top5) %>% mutate(sum = sum(relative_abundnace_sum)) %>% select(Day, sum) %>% unique()

test <- test %>% ungroup() %>% group_by(Day) %>% mutate(total_day_sum=sum(sum)) %>% mutate(relative_abundance = sum / total_day_sum) %>% select(-total_day_sum)

# Generating figure 3d
ggbarplot(test, x = "Day", y = "relative_abundance",
          fill = 'top5', color = "top5", legend = c('right'),
          palette = get_palette('npg',3), ylab = 'Relative abundance', legend.title="condition") +
  # rremove('y.text') +
  font("y.text", size=6)+
  rremove('ylab') +
  font("xlab", size = 6)+
  font('x.text', size =6) +
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5)) +
  rotate()

rm(list=ls())

#### (Sub fig4-a,b) taxo related with specific drug class? ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") #load name: n_bin_n_amr

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

amr_gene_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class
total_AMR_gene <- amr_gene_rgi_result %>% select(Run, Bin) %>% group_by(Run,Bin) %>% summarise(total_AMR_gene=n())

# Create an empty data frame with the columns
amr_gene_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(amr_gene_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(amr_gene_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = amr_gene_rgi_result$Run[i], Bin = amr_gene_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    amr_gene_drugclass <- bind_rows(amr_gene_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(amr_gene_drugclass) # number of row: 4333

amr_gene_drugclass <- amr_gene_drugclass %>% left_join(total_AMR_gene, by=c('Run','Bin'))
rm(total_AMR_gene)


merged_drug_class <- amr_gene_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup()

SARB_w_amr <- n_bin_n_amr %>% filter(AMR>0 & AMR<17) # 479

merged_drug_class <- merged_drug_class %>% left_join(antibiotics_abun_taxo, by=c("Run","Bin")) # all include
merged_drug_class <- merged_drug_class %>% inner_join(SARB_w_amr %>% select(-AMR), by=c("Run","Bin")) # only SARB with amr gene

phylum_count <- merged_drug_class %>% select(c("Run","Bin", "phylum")) %>% unique()  %>% group_by(phylum) %>% summarise(n=n()) # phylum 등장 횟수

test <- merged_drug_class %>% select(c("total_drugclass","Drug_Class","phylum")) %>% left_join(phylum_count, by='phylum')

test <- test %>%
  select(c("total_drugclass","Drug_Class","phylum","n")) %>%
  group_by(phylum, Drug_Class) %>%
  summarise(total_amr = sum(total_drugclass, na.rm = TRUE)/n) %>%
  arrange(Drug_Class) %>% 
  unique()


# Data transformation to wide format
wide_data <- test %>%
  spread(key = Drug_Class, value = total_amr, fill = 0)

annotation_data <- wide_data$phylum

# Removing the phylum column
wide_data$phylum <- NULL

wide_data <- log(wide_data + 1)
rownames(wide_data) <- annotation_data

# Generating the heatmap
pheatmap(wide_data, fontsize = 6)

rm(list=ls())

### Now for EARB - specie level ###

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata") #load name: amr_gene_rgi_result; total amr genes = 1635
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata") #load name: n_bin_n_amr
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata") #load name: n_bin_n_amr

antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Firmicutes"), "p__Firmicutes",phylum)) # Firmicutes subtypes to Firmicutes
antibiotics_abun_taxo <- antibiotics_abun_taxo %>% mutate(phylum=ifelse(str_detect(phylum,"p__Desulfobacterota"), "p__Desulfobacterota",phylum)) # p__Desulfobacterota_I -> p__Desulfobacterota 

amr_gene_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class
total_AMR_gene <- amr_gene_rgi_result %>% select(Run, Bin) %>% group_by(Run,Bin) %>% summarise(total_AMR_gene=n())

# Create an empty data frame with the columns
amr_gene_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(amr_gene_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(amr_gene_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = amr_gene_rgi_result$Run[i], Bin = amr_gene_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    amr_gene_drugclass <- bind_rows(amr_gene_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(amr_gene_drugclass) # number of row: 4333

amr_gene_drugclass <- amr_gene_drugclass %>% left_join(total_AMR_gene, by=c('Run','Bin'))
rm(total_AMR_gene)


merged_drug_class <- amr_gene_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup()

EARB_w_amr <- n_bin_n_amr %>% filter(AMR>16) # 20

merged_drug_class <- merged_drug_class %>% left_join(antibiotics_abun_taxo, by=c("Run","Bin")) # all include
merged_drug_class <- merged_drug_class %>% inner_join(EARB_w_amr %>% select(-AMR), by=c("Run","Bin")) # only EARB with amr gene

specie_count <- merged_drug_class %>% select(c("Run","Bin", "specie")) %>% unique()  %>% group_by(specie) %>% summarise(n=n()) # specie 등장 횟수

test <- merged_drug_class %>% select(c("total_drugclass","Drug_Class","specie")) %>% left_join(specie_count, by='specie')

test <- test %>%
  select(c("total_drugclass","Drug_Class","specie","n")) %>%
  group_by(specie, Drug_Class) %>%
  summarise(total_amr = sum(total_drugclass, na.rm = TRUE)/n) %>%
  arrange(Drug_Class) %>% 
  unique()


# Data transformation to wide format
wide_data <- test %>%
  spread(key = Drug_Class, value = total_amr, fill = 0)

annotation_data <- wide_data$specie

# Removing the specie column
wide_data$specie <- NULL

wide_data <- log(wide_data + 1)
rownames(wide_data) <- annotation_data

# Generating the heatmap
pheatmap(wide_data, fontsize = 6)

rm(list=ls())

#### (Sub fig5-a) resistome 기준 clustering ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_abundance_taxo.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/antibiotics_rgi_result.Rdata")
load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/n_bin_n_amr.Rdata")

amr_gene_rgi_result <- amr_gene_rgi_result %>% select(Run,Bin,`Drug Class`) %>% mutate(`Drug Class` = str_replace_all(`Drug Class`, "; ", "\n")) # select drug class
total_AMR_gene <- amr_gene_rgi_result %>% select(Run, Bin) %>% group_by(Run,Bin) %>% summarise(total_AMR_gene=n())

# Create an empty data frame with the columns
amr_gene_drugclass <- tibble(Run = character(), Bin = character(), Drug_Class = character())

# Loop over each row in 'test' data frame
for (i in 1:nrow(amr_gene_rgi_result)) {
  
  # Split the 'Drug Class' column into individual drugs
  drug_classes <- unlist(strsplit(as.character(amr_gene_rgi_result$`Drug Class`[i]), "\n"))
  
  # Loop over each drug_class
  for (drug_class in drug_classes) {
    
    # Create a new row as a tibble
    new_row <- tibble(Run = amr_gene_rgi_result$Run[i], Bin = amr_gene_rgi_result$Bin[i], Drug_Class = drug_class)
    
    # Add the new row to the output data frame
    amr_gene_drugclass <- bind_rows(amr_gene_drugclass, new_row)
  }
}
rm(new_row)

# Print the output data frame
print(amr_gene_drugclass) # number of row: 4333

amr_gene_drugclass <- amr_gene_drugclass %>% group_by(Run, Bin, Drug_Class) %>% summarise(total_drugclass=n()) %>% ungroup()

test <- amr_gene_drugclass %>% left_join(n_bin_n_amr %>% filter(AMR>16) %>% select(-AMR) %>% mutate(EARB='yes'), by=c('Run','Bin'))

test[is.na(test)] <- 'no'

test <- test %>% left_join(antibiotics_abun_taxo, by=c('Run','Bin'))

final_result <- test %>% select(Run, Bin, Drug_Class, total_drugclass, EARB, specie)

specie_count <- final_result %>% select(c("Run","Bin", "specie")) %>% unique()  %>% group_by(specie) %>% summarise(n=n()) # specie 등장 횟수

test <- final_result %>% select(c("total_drugclass","Drug_Class","specie")) %>% left_join(specie_count, by='specie')

test <- test %>%
  select(c("total_drugclass","Drug_Class","specie","n")) %>%
  group_by(specie, Drug_Class) %>%
  summarise(total_amr = sum(total_drugclass, na.rm = TRUE)/n) %>%
  arrange(drug_class) %>% 
  unique()

# Data transformation to wide format
wide_data <- test %>%
  spread(key = Drug_Class, value = total_amr, fill = 0)

annotation_data <- wide_data$specie

# Removing the phylum column
wide_data$specie <- NULL

wide_data <- log(wide_data + 1)
rownames(wide_data) <- annotation_data

#### kmean clustering

library(fpc)
library(factoextra)

km.res <- kmeans(wide_data, 2, nstart = 5)

# Rename clusters
km.res$cluster[km.res$cluster == 1] <- "EARB"
km.res$cluster[km.res$cluster == 2] <- "SARB"

fviz_cluster(km.res, wide_data, ellipse.type = "norm", geom="point", main = 'k-means clustering')

km.res$cluster[km.res$cluster=='EARB']

rm(list=ls())

#### (Sub fig6-a) NMF ####
setwd('/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/')

metadataRData = "antibiotic_metadata.Rdata"
load(metadataRData)

geneCountFile = "Drugclass_sample_for_NMF.tsv"
geneCountTab = read.delim(geneCountFile)
geneCountMat = geneCountTab[,-1]
rownames(geneCountMat) = geneCountTab[,1]

require(NMF)
nmfSeed("nn")
ranks = 2

nOut_rank2 = nmf(geneCountMat, 2, seed=1) 

basis_rank2 = basis(nOut_rank2)
coef_rank2 = coef(nOut_rank2)

colnames(basis_rank2) = c("sig1", "sig2")
rownames(coef_rank2) = c("sig1", "sig2")

require(ggplot2)
require(reshape2)
basis_rank2_melt = melt(basis_rank2)
coef_rank2_melt = melt(coef_rank2)
coef_rank2_melt$Subject = tmp_metadata$sample[match(coef_rank2_melt$Var2, tmp_metadata$Run)]
coef_rank2_melt$day = tmp_metadata$condition[match(coef_rank2_melt$Var2, tmp_metadata$Run)]
coef_rank2_melt$day = factor(coef_rank2_melt$day, levels = c("Dag0", "Dag4", "Dag8", "Dag42", "Dag180"))

coef_rank2_melt$day <- gsub('Dag','Day', coef_rank2_melt$day)
coef_rank2_melt$Subject <- gsub('ERAS','H', coef_rank2_melt$Subject)
coef_rank2_melt$day <- gsub('Day','D', coef_rank2_melt$day)
coef_rank2_melt$day <- factor(coef_rank2_melt$day, levels = c('D0','D4','D8','D42','D180'))
coef_rank2_melt$Subject <- factor(coef_rank2_melt$Subject, 
                                  levels = paste("H", 1:12, sep=""))

ggbarplot(data=coef_rank2_melt, x='Subject', y='value', ylab = 'NMF value',  fill = "Var1", palette = "npg", legend = 'bottom', size = 0.1) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  facet_grid(~ day)+
  font("xlab", size = 6) +
  font("x.text", size = 6) +
  font("ylab", size = 6)+
  font("y.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6)+
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(strip.text.x = element_text(size = 6)) +
  labs(fill = 'coef rank')

rm(list=ls())

#### (Sub fig7-a) EARB heatmap ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/ANI_matrix.Rdata") #load name: ANI_matrix

colnames(ANI_matrix) <- gsub('_D',' D',colnames(ANI_matrix))
rownames(ANI_matrix) <- gsub('_D',' D',rownames(ANI_matrix))

pheatmap(ANI_matrix,display_numbers = F, fontsize = 6, main = "inter-species ANI matrix")

# Don't remove environment #

#### (Sub fig7-b) EARB E.coli ANI score compair with H5, H8 ####

E.coli_list <- ANI_matrix %>% select(`H7 D8.1`) %>% filter(`H7 D8.1`>95) %>% rownames()

ANI_matrix <- ANI_matrix %>% select(E.coli_list) %>% mutate(tmp=rownames(ANI_matrix)) %>% filter(tmp %in% E.coli_list) %>% select(-tmp)

ANI_matrix <- ANI_matrix %>% select(c("H7 D8.1",  "H2 D42",  "H1 D8.1",  "H9 D8",  "H10 D8", "H5 D4",  "H8 D8",  "H12 D4", "H12 D8", "H5 D8",  "H4 D42",  "H6 D8" ))
ANI_matrix <- ANI_matrix[c("H7 D8.1",  "H2 D42",  "H1 D8.1",  "H9 D8",  "H10 D8", "H5 D4",  "H8 D8",  "H12 D4", "H12 D8", "H5 D8",  "H4 D42",  "H6 D8"),]

colnames(ANI_matrix) <- gsub('D8.1','D8',colnames(ANI_matrix)) 
rownames(ANI_matrix) <- gsub('D8.1','D8',rownames(ANI_matrix))

ANI_compare <- as.data.frame(ANI_matrix) %>% select(c("H12 D4", "H12 D8","H5 D4","H5 D8"))

# 행 이름을 열로 변환
ANI_compare$Sample <- rownames(ANI_compare)

# melt 함수 적용
melted_data <- melt(ANI_compare, id.vars = "Sample", measure.vars = c("H12 D4", "H12 D8"))

melted_data <- melted_data %>%
  mutate(HOST_ID = as.numeric(sub("H([0-9]+).*", "\\1", Sample))) %>%
  arrange(HOST_ID) %>%
  select(-HOST_ID)

ggbarplot(melted_data, x = 'Sample', y = 'value', fill = 'variable',
          x.text.angle = 45, xlab = 'EARB E. coli strains',
          position = position_dodge2(preserve = "single"), 
          legend = "top", 
          palette = get_palette(palette = 'npg',4), 
          ylab = 'ANI score', 
          ylim = c(96, 100)) +
  theme(strip.text.x = element_text(size = 6))+
  font("xy.text", size = 6) +
  font("legend.title", size = 6) +
  font("legend.text", size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  labs(fill = bquote("Reference EARB " * italic("E. coli") * " strain"),
       x = bquote("EARB " * italic("E. coli") * " strains"))



# melt 함수 적용
melted_data <- melt(ANI_compare, id.vars = "Sample", measure.vars = c("H5 D4", "H5 D8"))

melted_data <- melted_data %>%
  mutate(HOST_ID = as.numeric(sub("H([0-9]+).*", "\\1", Sample))) %>%
  arrange(HOST_ID) %>%
  select(-HOST_ID)

ggbarplot(melted_data, x = 'Sample', y = 'value', fill = 'variable',
          x.text.angle = 45, xlab = 'EARB E. coli strains',
          position = position_dodge2(preserve = "single"), 
          legend = "top", 
          palette = get_palette(palette = 'npg',4), 
          ylab = 'ANI score', 
          ylim = c(96, 100)) +
  theme(strip.text.x = element_text(size = 6))+
  font("xy.text", size = 6) +
  font("legend.title", size = 6) +
  font("legend.text", size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  font("xylab", size = 6) +
  labs(fill = bquote("Reference EARB " * italic("E. coli") * " strain"),
       x = bquote("EARB " * italic("E. coli") * " strains"))

rm(list=ls())

#### (Sub fig7-c) 각 contig의 평균적인 depth가 궁금하군요. (N) ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D4_AMR_all_position.Rdata") #load name: H5_D4_AMR_AP

H5_D4_AMR_AP <- H5_D4_AMR_AP %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .)))

colname <- colnames(H5_D4_AMR_AP) # 1841
modified_colname <- gsub("_[^_]*$", "", colname) %>% unique() # column 이름중에서 position 정보를 지우고 contig의 이름들만 남겼다.
modified_colname <- modified_colname[-55]

result <- data.frame(matrix(nrow=0, ncol=3))

for (name in modified_colname) {
  test <- H5_D4_AMR_AP %>% select(starts_with(name)) %>% summarise(across(everything(), sum, na.rm = TRUE))
  column_sums <- colSums(test, na.rm = TRUE)
  sum_data <- data.frame(column_names = names(column_sums), sum_values = column_sums)
  sum_data$contig <- name
  
  print(name)
  print(summary(sum_data$sum_values))
  
  result <- rbind(result, sum_data)
}

contigs <- c('k113_10194_length_520744_cov_43.6593_402','k113_16363_length_314213_cov_52.9977_282','k113_18543_length_373955_cov_48.0311_3','k113_31323_length_358434_cov_43.9983_245','k113_7590_length_24589_cov_63.5384_4')

result$category <- ifelse(result$contig %in% contigs, "SNP not exist","SNP exist")

ggplot(result, aes(x = contig, y = sum_values, fill = category)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("SNP exist" = 'blue', 'SNP not exist'='red'), name = "Contigs") +
  theme(axis.text.x = element_blank(), 
        strip.text = element_blank()) +
  labs(x = "",
       y = "Sequencing depth", ) +
  font("ylab", size = 6)+
  font("y.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm'), legend.position = 'bottom', axis.title.x = element_blank())

rm(list=ls())


#### (Sub fig8-a) Proportion of Major, Minor, and None for Each Sample (it takes some time) ####

load("/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/H5_D4_AMR.Rdata") # 1841

H5_D4_AMR <- H5_D4_AMR %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .))) 

result <- data.frame(matrix(nrow=0, ncol=5))

for (position in colnames(H5_D4_AMR)) {
  position_depth <- sum(H5_D4_AMR[[position]])
  nucleotide_freq <- sort(H5_D4_AMR[[position]], decreasing = T)
  
  if (position_depth==0) {
    print('position_depth is zero')
    break
  }
  
  if (nucleotide_freq[1] == 0) {
    print('major position is zero')
    break
  } else if (nucleotide_freq[2] == 0) {
    print('minor position is zero')
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='one allele',
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    next
  } else {
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[1]/position_depth, 'condition'='major',
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[1])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
    
    temp_result <- data.frame('position'=position, 'SNV_ratio'=nucleotide_freq[2]/position_depth, 'condition'='minor', 
                              'nucl'=rownames(H5_D4_AMR)[which(H5_D4_AMR[[position]] == nucleotide_freq[2])], 'status'='H5 D4')
    result <- rbind(result, temp_result)
  }
  
}

table(result$condition) # major: 1841 / minor: 1841

result_47_48amr  <- result

major_nucl_pos <- result_47_48amr %>% filter(condition=='major') %>% select(position, nucl) %>% separate(position, into = c("contig", "pos"), sep = "_(?!.*_)")
minor_nucl_pos <- result_47_48amr %>% filter(condition=='minor') %>% select(position, nucl) %>% separate(position, into = c("contig", "pos"), sep = "_(?!.*_)")

major_nucl_pos$contig <- gsub('\\.','_', major_nucl_pos$contig)
minor_nucl_pos$contig <- gsub('\\.','_', minor_nucl_pos$contig)

file_list <- list.files(path = "/Users/jwbaek9506/Documents/analysis/EARB_project/final_data_for_figure/SUP 8", full.names = T)
final <- NULL

for (file in file_list) {
  
  test <- fread(file = file, fill=T)
  test <- as.data.frame(test)
  test <- test %>% select(-tem)
  rownames(test) <- test$nt
  test <- test[,-1]
  test <- test %>% select(-ncol(test))
  test <- test %>% mutate(across(everything(), ~ ifelse(. < sum(.) * 0.05, 0, .)))
  colnames(test) <- gsub(pattern = '\\.', replacement = '_',x = colnames(test))
  
  nucl_pos <- major_nucl_pos %>% full_join(minor_nucl_pos, by=c('contig','pos'), suffix=c("_major","_minor"))
  
  result <- data.frame(matrix(nrow=0, ncol=4))
  final_result <- data.frame(matrix(nrow=0, ncol=3))
  
  for (i in 1:ncol(test)) {
    name <- colnames(test)[i]
    max_row <- which.max(test[, name])
    max_row_name <- rownames(test)[max_row]
    split_string <- strsplit(name, "_")[[1]]
    part1 <- paste(split_string[1:length(split_string)-1], collapse = "_")
    part2 <- split_string[length(split_string)]
    tmp <- nucl_pos %>% filter(contig == part1 & pos == part2) %>% mutate('result'=ifelse(max_row_name == toupper(nucl_major), 'major', ifelse(max_row_name == toupper(nucl_minor), "minor","none")))
    result <- rbind(result,tmp)
  }
  final_result <- rbind(final_result, data.frame('major'=table(result$result)['major'],'minor'=table(result$result)['minor'],'none'=table(result$result)['none'], row.names = regmatches(file, regexpr("SRR\\d+", file))))
  final <- rbind(final, final_result)
}

temp <- final

final <- final %>% mutate(none = replace_na(none, 0))

# Assuming your dataframe is named 'final'
final$Sample <- rownames(final)

# Calculate the proportion of 'major' for each sample
final <- final %>%
  mutate(total = major + minor + none,
         major_prop = major / total) %>%
  arrange(desc(major_prop))  # Sort by major proportion in descending order

# Reshape the data from wide to long format
final_long <- pivot_longer(final, cols = c(major, minor, none), 
                           names_to = "Category", values_to = "Count")

# Create the sorted stacked bar chart
ggplot(final_long, aes(x = reorder(Sample, -major_prop), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of Major, Minor, and None for Each Sample",
       x = "Sample", y = "Proportion of SNP profile") +
  theme_minimal() +
  scale_fill_manual(values = c("major" = "#FF9999", "minor" = "#66CC66", "none" = "#6699CC")) +
  rremove('xlab') +
  rremove('x.text')+
  font("y.text", size=6) +
  font('ylab', size=6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm'))

final_long %>% filter(major_prop > 0.5) %>% select(Sample) %>% unique() # 360

rm(list=ls())

