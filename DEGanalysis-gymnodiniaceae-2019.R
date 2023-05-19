#Differential expression during 2018
rm(list = ls())
source("~/Dropbox/R_general_working_directory/Dictionary_of_Colors.R")

library(tidyverse); library(edgeR); library(reshape2); library(WGCNA); library(patchwork);library(jishonoiro);library(vegan);library(corrplot);library(DESeq2)
setwd("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2019_assembly/")
#save.image(file = "DEA_2018-edgeR-DESeq2.Rdata")
#load("DEA_2018-edgeR-DESeq2.Rdata")

one41= c("#FF5200", "#0D2B52","#A6FF47")
pal1 = c("#8c6510", "#0D2B52","#FF7399")
color_pal = c("black", "light grey", "steel blue")

# Load and process the raw count table ------------------------------------------
#filter to the dominant diatom taxa from 2018 with enough reads to meaningfully do so (Figure S#-diatoms). I left off Rhizosolenia because there just weren't enough reads.
Dino_raw_df = read.delim("raw-count-data-metaT-SMP2019.txt") %>%
  filter(str_detect(Taxonomy, regex("Gymnodiniaceae", ignore_case = T))) %>%  #filter for dinophyceae (akashiwo, cochlodinium) reads using string detect
  mutate(tax = case_when(str_detect(Taxonomy, "Gymnodiniaceae") ~"Gymnodiniales")) %>% 
  select(tax, KO, everything()) %>% 
  unite("ID", c(tax,KO), sep = "-") %>% 
  dplyr::select(-SMPier.Day1.C.04.13.2019_S11, -SMPier.Day12.C.04.24.2019_S13, -Taxonomy)

colsplit(Dino_raw_df$Taxonomy, pattern = ";", names = c("L1","L2","L3","L4","L5","L6","L7","L8")) %>% 
  mutate(L7 = ifelse(L7 == "", "XXX", L7),
         L6 = ifelse(L6 == "", "XXX", L6),
         L5 = ifelse(L5 == "", "XXX", L5))  %>% 
  dplyr::count(L5)

Dino_raw_df %>% nrow() # there are approx 8000 genes transcripts
colSums(Dino_raw_df[3:14]) %>% sum() # approx 1.2 Million transcripts across 8000 genes
colSums(Diatom_raw_df[3:17]) %>% sum() # approx 2.5 Million transcripts across 26,000 genes
nrow(Dino_raw_df)
#Aggregate by KO
agg_df = Dino_raw_df %>%
  pivot_longer(cols = starts_with("SMP"), names_to = "sample", values_to = "value") %>% 
  #dcast(ID~sample) %>% column_to_rownames("ID") %>% colSums()
  mutate(sample = str_replace(sample, "SMPier\\.", ""),
         sample = str_replace(sample, regex("\\.\\d{2}\\.\\d{2}\\.2019_S\\d+"), ""),
         sample = factor(sample, levels = c("Day1.A", "Day1.B", 
                                            "Day4.A", "Day4.B", "Day4.C", 
                                            "Day12.A", "Day12.B", 
                                            "Day18.A", "Day18.B", "Day18.C", 
                                            "Day20.A", "Day20.B", "Day20.C"))) %>% 
  group_by(ID, sample) %>%
  summarise(count = sum(value, na.rm = TRUE)) %>%
  dcast(ID~sample) #%>% 
  select(-starts_with("Day4."))

#Assign rownames as KO
rownames(agg_df) = agg_df$ID

#sanity check. These should be the same
colSums(agg_df[-1]) == colSums(Dino_raw_df[-1])

# # Normalize the agg_df dataframe with edgeR -----------------------------
#(1) create a dge object and (2) calculate normalizing factors and then (3) return the cpm as a dataframe
#Ingredients for the dge object: dataframe (counts only!), sample_list, & metadata.

#Make sure these line up with the order of the columns of the agg_df
# Create a sample list 
sample_list = c("Day1","Day4","Day12","Day18","Day20")
# Create the metadata
metadata = factor(c(rep("day1",2),
                    rep("day4",3),
                    rep("day12",2),
                    rep("day18",3),
                    rep("day20",3)),
                  levels = str_to_lower(sample_list)) # levels is important for purposefully setting the sample order.


# Use this information to create the DGEList object
dge_obj = DGEList(counts = as.data.frame.matrix(agg_df[-1]), # Counts is the actual columns of the dataframe
                  genes = agg_df[1], # This is the name of the column(similar to genes)
                  group = metadata)

#calculate the normalizing factors using TMM
dge_obj = calcNormFactors(dge_obj, method = "TMM")

keep = filterByExpr(dge_obj)
dge_obj = dge_obj[keep, ,keep.lib.sizes=FALSE]
  # setup the design matrix without an intercept as day 1!
design = model.matrix(~0+group, data = dge_obj$samples)
  
  #Set the column names for the design to match the
colnames(design) = levels(dge_obj$samples$group)
  
  # using makeContrasts() to specify the pariwise comparisons
conts = makeContrasts(day4-day1, day12-day1, day18-day1, day20-day1, levels = str_to_lower(sample_list)) #for contrasting all samples to prebloom state
  
  # estimate common and tagwise(pairwise) dispersion accrding the study design matrix
disp = estimateGLMCommonDisp(dge_obj, design = design)
disp = estimateGLMTagwiseDisp(disp, design = design)
  
  # Determine differentially expressed genes
fit = glmQLFit(disp, design = design)
edgeR_DEGs = glmQLFTest(fit, contrast = conts)
edgR_DEG_table = edgeR_DEGs$table %>% data.frame()
  
  # adjust p values due to so many comparissons to decrease chance of type 1 error.
edgR_DEG_table$P.adjusted = p.adjust(edgR_DEG_table$PValue, method = "fdr")
edgR_DEG_table$KO = rownames(edgR_DEG_table)
edgR_DEG_table_Pfilt = edgR_DEG_table %>% 
  mutate(ID = rownames(.)) %>%
  filter(P.adjusted < 0.001) %>% 
  separate(KO, into = c("tax","KO"), sep = "-")
  


#Next goal is to calculate differential expression using DESeq2 for both taxa in the way I used edgeR.
# End edgeR

# Start DESeq
# DESeq2 ------------------------------------------------------------------
library(DESeq2)

#Create data tables for the different taxa from the raw input

# Create the metadata
infodf = data.frame(sample = str_to_title(c("day1.A","day1.B",
                                            "day4.A","day4.B","day4.C",
                                            "day12.A","day12.B",
                                            "day18.A","day18.B","day18.C",
                                            "day20.A","day20.B","day20.C")),
                      condition = factor(c(rep("day1",2),rep("day4",3),rep("day12",2),rep("day18",3),rep("day20",3)), 
                                         levels = c("day1","day4","day12","day18","day20")),
                      row.names = "sample") %>% # levels is important for purposefully setting the sample order.
  mutate(condition = as.factor(condition))

#Verify sample names and column names
rownames(infodf) %in% colnames(agg_df)


#Create the DESeqDataSet for DGE analysis
#note: specifies the experimental design upon which to construct the model. 
#Its used to estimate dispersion and calcultat log fold changes.
DseqObj = DESeqDataSetFromMatrix(countData = round(agg_df[-1]),
                             colData = infodf,
                             design = ~condition)

#prefilter low counts- here genes with less than 10 reads
DseqObj = DseqObj[rowSums(counts(DseqObj)) >= 10,]

#select reference level. Comparisons of other conditions will be using this reference.
DseqObj$condition = relevel(DseqObj$condition, ref = "day1")
design(DseqObj) = ~0 + condition #set the design to be without an intercept
#perform differential gene expression analysis,
DseqObj = DESeq(DseqObj, betaPrior = FALSE) # betaPrior because no intercept

# End DESeq2

# Get both tables together.

resultsNames(DseqObj) #Get the list of result names to be contrasted

## Create a function to do this in one step
contrast_function = function(Condition1, Condition2) {
  contrast_label = str_c(Condition1,Condition2,sep = "")
  
  #Produce the DESeq table with the desired contrasts.
  res = results(DseqObj, contrast = list(c(Condition1), c(Condition2))) #get the results for the comparisons of interest in list form. Here I am contrasting day3 vs day1
  
  # Convert DESeq table into the appropriate form: splitting the taxa
  Deseq_DEG_table = res %>% 
    data.frame %>% 
    mutate(ID = rownames(.)) %>% 
    separate(ID, into = c("tax", "KO"), sep = "-", remove = F)
  
  ## Join with the edgeR table to Compare DESeq2 degs vs edgeR differential expration
  joined = edgR_DEG_table_Pfilt %>% 
    dplyr::select(logFC.day4...day1, P.adjusted, ID) %>% #These need to be manually changed to represent the day in the edgeR sample
    inner_join(., Deseq_DEG_table %>% drop_na() %>% select(log2FoldChange, padj, ID), by = "ID") %>% 
    filter(padj < 0.01)
  
  #write.csv(joined, file = "DEGS_d3VS1-joined2.csv", row.names = F)
  
  output = joined %>% 
    mutate(check = abs(logFC.day4...day1-log2FoldChange) > 1, #These need to be manually changed to represent the day in the edgeR sample
           LogFC = (logFC.day4...day1+log2FoldChange)/2, #These need to be manually changed to represent the day in the edgeR sample
           Contrast = contrast_label) %>% 
    separate(ID, into = c("tax", "KO"), sep = "-", remove = F) %>% 
    select(-contains("adj"), -contains("log", ignore.case = F)) %>% 
    mutate(Contrast = str_replace_all(contrast_label, "condition",""))
  
  return(output)
}

#Run the function for each contrast. REMEMBER TO MANUALLY CHANGE THE SAMPLES FOR THE EDGER LOGFC COLUMN.
d4vd1 = contrast_function("conditionday4","conditionday1")
d12vd1 = contrast_function("conditionday12","conditionday1")
d18vd1 = contrast_function("conditionday18","conditionday1")
d20vd1 = contrast_function("conditionday20","conditionday1")

DEG_table_joined = bind_rows(d4vd1, d12vd1, d18vd1, d20vd1) %>% 
  mutate(Contrast = str_replace_all(Contrast, "day1$",""))
save(DEG_table_joined, file = "DEG_table_d1_joined-2019.rda") #9May
load("DEG_table_d1_joined-2019.rda")
## What are the questions?
sample_list
#1. What is the number of genes that illustrated differential expression? And How did the number change with time?
#barplot
All_gymnodegs_2019 = DEG_table_joined %>% 
  group_by(Contrast) %>% 
  dplyr::count(KO) %>% 
  summarise(n_genes = sum(n)) %>% 
  ungroup() %>% 
  mutate(Year = "2019",
         Contrast = factor(Contrast, levels = c("day4","day12","day18","day20"))) %>% 
  ggplot(., aes(x = Year, y = n_genes, fill = Contrast)) +
  geom_bar(stat = "identity", color = "white") +
  scale_fill_manual(name = NULL,
                    values = jisho_picker("la_lakers")[c(3,1,2,6)]) +
  labs(x = "",
       y = "Number of DEGs") +
  theme_minimal()
ggsave(last_plot(), filename = "n_degs_d1_w4_stackbar-2019.pdf", height = 5, width = 4)

gymnoPos_neg_bars_2019 = DEG_table_joined %>% 
  mutate(Year = "2019",
         Contrast = factor(Contrast, levels = c("day4","day12","day18","day20")),
         Direction = ifelse(LogFC > 0, "Positive", "Negative"),
         Direction  = factor(Direction, levels = c("Positive", "Negative"))) %>% 
  group_by(Direction, Contrast) %>% 
  dplyr::count(Contrast) %>% 
  ungroup() %>% 
  ggplot(aes(x = Contrast, y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", color = "white") +
  scale_fill_manual(name = "",
                    values = jisho_picker('la_lakers')) +
  labs(y = "Number of DEGs",
       x = "") +
  theme_minimal() 
All_gymnodegs_2019 + Dir_waffle + gymnoPos_neg_bars_2019 +  patchwork::plot_layout(ncol = 2)
ggsave(last_plot(), filename = "Figure_S1.pdf", height = 5, width = 9)

## Table to look at
gymno2019_degstats =DEG_table_joined %>% 
  mutate(Year = "2019",
         Contrast = factor(Contrast, levels = c("day4","day12","day18","day20")),
         Direction = ifelse(LogFC > 0, "Positive", "Negative"),
         Direction  = factor(Direction, levels = c("Positive", "Negative"))) %>% 
  group_by(Contrast, Direction) %>% 
  dplyr::count(KO) %>% 
  summarise(n_genes = sum(n)) %>% 
  ungroup() %>% 
  mutate(Total = sum(n_genes)) %>% 
  group_by(Contrast) %>% 
  mutate(proportion = n_genes/Total,
         dayProportion = sum(proportion),
         Contrast = factor(Contrast, levels = c("day4","day12","day18","day20"))) 
write.csv(gymno2019_degstats, file = "gymn2019_degstats.csv")

Dir_waffle = DEG_table_joined %>% 
  mutate(Year = "2019",
         Contrast = factor(Contrast, levels = c("day4","day12","day18","day20")),
         Direction = ifelse(LogFC > 0, "Positive", "Negative"),
         Direction  = factor(Direction, levels = c("Positive", "Negative"))) %>% 
  group_by(Direction, Contrast) %>% 
  dplyr::count(Contrast) %>% 
  group_by(Direction) %>% 
  ggplot(., aes(fill = Direction, values = n)) +
  geom_waffle(make_proportional = FALSE, color = "white") +
  scale_fill_manual(name = "",
                    values = jisho_picker('la_lakers')) +
  labs(x = "",
       y = "Transcripts") +
  theme_void()
ggsave(last_plot(), filename = "Figure_S2c.pdf", height = 4.5, width = 6)

save.image("QuicDEGstats.Rdata")
#waffle plot
DEG_table_joined %>% 
  mutate(Contrast = str_to_title(Contrast)) %>% 
  group_by(Contrast) %>% 
  dplyr::count(KO) %>% 
  summarise(n_genes = sum(n)) %>% 
  ungroup() %>% 
  mutate(total = sum(n_genes),
         Contrast = factor(Contrast, levels = c("Day4","Day12","Day18","Day20"))) %>% 
  ggplot(., aes(x = Contrast, y = n_genes)) +
  geom_bar(fill = "grey30", stat = "identity", width = .5) + 
  #scale_y_continuous(breaks = c(0,seq(100,500,100))) +
  labs(x = "",
       y = "Transcripts") +
  theme_minimal()
ggsave(last_plot(), filename = "n_degs_d1_w4dodgebar-2019.pdf", height = 4, width = 5)

## Waffle plot?
library(waffle)

DEG_table_joined %>% 
  #mutate(Contrast = str_to_title(str_replace_all(Contrast, "day4$",""))) %>% 
  group_by(Contrast) %>% 
  dplyr::count(KO) %>% 
  summarise(n_genes = sum(n)) %>% 
  ungroup() %>% 
  mutate(total = sum(n_genes),
         Contrast = factor(Contrast, levels = str_to_title(sample_list))) %>%
  arrange(Contrast) %>% #because the factor statement didn't work to arrange things.
  ggplot(., aes(fill = Contrast, values = n_genes)) +
  geom_waffle(make_proportional = FALSE) +
  scale_fill_manual(name = NULL,
                    values = jisho_picker("los_angeles")[c(3,1,4,5)]) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks
ggsave(last_plot(), filename = "n_degs_d1_w4waffle-2019.pdf", height = 5, width = 6)

# 

#1. What are the identities of the genes with the highest log fold changes?
#1: I think the identities with the higest log fold change will be related to nutrient uptake.

Counts_to_fx_table = function(DEG_TABLE) {
  #Purpose: add functional annotations to a numerica table using the KO ids while cutting down on reundancy through manual curation
  table_wcat = DEG_TABLE %>%  
    left_join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T), by = "KO") %>%
    left_join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t"), by = "KO")      
  
  ##Create vectors containing all KOs that corresponded for a given pathway
  #glycolysis
  Glycolysis_KOs = table_wcat %>% 
    filter(str_detect(Category, regex("glycol", ignore_case = T))) %>% 
    pull(KO)
  
  #Calvin cycle 
  Calvin_KOs = table_wcat %>% 
    filter(str_detect(Category, regex("calvin", ignore_case = T))) %>% 
    pull(KO)
  
  #Fatty acid breakdown
  Fat_breakdown_KOs = table_wcat %>% 
    filter(str_detect(Category, regex("fatty acid breakdown", ignore_case = T))) %>% 
    pull(KO)
  
  #Fatty acid boisynthesis
  Fat_biosynth_KOs = table_wcat %>% 
    filter(str_detect(Category, regex("fatty acid biosynthesis", ignore_case = T))) %>% 
    pull(KO)
  
  # Annotate the categories according to their pathway where possible
  table_wcat = table_wcat %>% 
    mutate(Category = ifelse(str_detect(Gene_identification, regex(" npt", ignore_case = T)), "P metabolism", Category),
           case_when(str_detect(Category, regex("(Glyoxylate)|(TCA)")) ~"TCA-Glyoxylate cycle",
                     KO %in% Glycolysis_KOs ~"Glycolysis-gluconeogenesis",
                     KO %in% Calvin_KOs ~"Calvin cycle",
                     KO %in% Fat_biosynth_KOs ~"Fatty acid biosynthesis",
                     KO %in% Fat_breakdown_KOs ~"Fatty acid breakdown",
                     TRUE ~ Category)) %>%
    mutate(Category = case_when(KO %in% Calvin_KOs[which(Calvin_KOs %in% Glycolysis_KOs)] ~"Glycolysis-Calvin",
                                KO %in% Fat_biosynth_KOs[which(Fat_biosynth_KOs %in% Fat_breakdown_KOs)] ~"Lipid biosynthesis-breakdown",
                                TRUE ~ Category),
           Category = ifelse(Category == "Gluconeogenesis", "Glycolysis-gluconeogenesis", Category)) #%>% 
  #dplyr::select(Day,KO,sum_TPM,B,Target_2,Category) 
  
  ## Now i want to substitute the empty categories for Target_2
  table_wcat[is.na(table_wcat)] = "XXX" # trick to convert all NA values to a string.
  table_wcat_clean = table_wcat %>% 
    mutate(Category = ifelse(Category == "XXX", Target_2, Category)) %>% 
    mutate(Category = ifelse(Category == "XXX", B, Category),
           Category = ifelse(Category == "P uptake", "P metabolism", Category)) %>% 
    dplyr::select(-B, -Target_2) %>% 
    group_by(tax, KO, Contrast, Category) %>% 
    summarise(LogFC = mean(LogFC)) %>% 
    ungroup() %>% 
    filter(Category != "XXX")
  
  return(table_wcat_clean)
  
}

library(viridis)
#Add functional categories to the KO terms of the DEG table.
dino_table_wcat_clean = Counts_to_fx_table(DEG_table_joined)

save(dino_table_wcat_clean, file = "dinotablewcatclean_withd4-2.rda")#9May

catOrder = c("Actin polymerization", "additional breakdown", "Environmental information processing", 
             "Genetic information processing", "phosphatidylinositol", "por", "AMT", "Chitinase", "GS/GOGAT",
             "Inorganic N uptake and assimilation", "Nitrate reduction (assimilatory)", "NRT", "Nucleotide and amino acid metabolism",
             "Organic N uptake and assimilation", "P metabolism", "Urea cycle", "Energy metabolism", "Calvin cycle", "Carbohydrate and lipid metabolism",
             "Fatty acid biosynthesis", "Fatty acid breakdown", "Glycolysis", "Glycolysis-Calvin", "Glycolysis-gluconeogenesis", "Glyoxylate cycle",
             "Lipid biosynthesis-breakdown", "Metabolism", "PDH", "Photosynthesis", "TCA cycle", "Triacylglycerol biosynthesis",
             "Endocytosis", "Lysosome binding and processing", "Phagosome maturation","Motility and prey recognition", "SNARE complex", "V-type ATPase")
length(catOrder)
###try pheatmap
library(pheatmap)
dino_for_plot = dino_table_wcat_clean %>% 
  mutate(Category = factor(Category, levels = catOrder),
         Contrast = str_replace_all(Contrast, "day1$","")) %>% 
  dcast(Category~Contrast, fun.aggregate = mean) %>% 
  select(Category, day4, day12, day18, day20) 

dino_for_plot[is.na(dino_for_plot)] = 0

save(dino_for_plot, file = "dino_forPheatmap.rda")
#creating a secondary annotation on the pheatmap called process
fx_annotation_4rows = data.frame(Process = factor(c(rep("Other", 7), 
                                                    rep("Nutrient processing", 9), 
                                                    rep("Metabolism", 15), 
                                                    rep("Other catabolism", 6))))
rownames(fx_annotation_4rows) = c("Actin polymerization", "additional breakdown", "Environmental information processing", 
                                  "Genetic information processing", "Nucleotide and amino acid metabolism", "phosphatidylinositol", "por", "AMT", "Chitinase", "GS/GOGAT",
                                  "Inorganic N uptake and assimilation", "Nitrate reduction (assimilatory)", "NRT", 
                                  "Organic N uptake and assimilation", "P metabolism", "Urea cycle", "Energy metabolism", "Calvin cycle", "Carbohydrate and lipid metabolism",
                                  "Fatty acid biosynthesis", "Fatty acid breakdown", "Glycolysis", "Glycolysis-Calvin", "Glycolysis-gluconeogenesis", "Glyoxylate cycle",
                                  "Lipid biosynthesis-breakdown", "Metabolism", "PDH", "Photosynthesis", "TCA cycle", "Triacylglycerol biosynthesis",
                                  "Endocytosis", "Lysosome binding and processing", "Phagosome maturation", "Motility and prey recognition", "SNARE complex", "V-type ATPase")
fx_annotation_4rows

#read in table with environment data
env_data = read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Env_tables/Env_tables.csv")

#Creating a chlorophyll concentration annotation
Chl_annotation = env_data %>% 
  filter(Year == 2019,
         Type != "RelAbuncance") %>% 
  select(Day, Type, Value) %>% 
  pivot_wider(names_from = Type, values_from = Value)
Chl_annotation = Chl_annotation[c(1,3:5),] %>% data.frame() %>% select(-Day)
rownames(Chl_annotation) = c("day1","day12","day18","day20")

#Creating a relative abundance annotation
#Creating a Nitrate concentration annotation
#Creating a Phosphate concentration annotation
Annotation_colors = list(Process = c("Other" = "black", "Nutrient processing" = "grey40", "Metabolism" = "white", "Other catabolism" = "grey70"),
                         #here I think I'll have to use the cut() function for continuous variables. Or just a low and high
                         Chlorophyll = rev(c("#00204DFF", "#002A63FF", "#00326FFF", "#1D3B6DFF", "#32456BFF", "#424E6BFF", "#50576CFF", "#5C616EFF", "#686A71FF","#777776FF", "#838079FF", "white")), 
                         Nitrogen = rev(c("#00204DFF", "#002A63FF", "#00326FFF", "#1D3B6DFF", "#32456BFF", "#424E6BFF", "#50576CFF", "#5C616EFF", "#686A71FF","#777776FF", "#838079FF", "white")), 
                         Phosphorus = rev(c("#00204DFF", "#002A63FF", "#00326FFF", "#1D3B6DFF", "#32456BFF", "#424E6BFF", "#50576CFF", "#5C616EFF", "#686A71FF","#777776FF", "#838079FF", "white")))

#the color gradient for the annotations are incorrect because they are scaling from vlaues aken on day 4, not day 1.
#I'm going to add a dummy column to the for_plot table for day 1 to ameliorate and clip manually in illustrator.
for_plot$day1 = 0
for_plot = for_plot %>% select(day1, day12, day18, day20)
DegPheatmap = pheatmap(for_plot, angle_col = 0,
         cluster_cols = F,
         cluster_rows = F,
         color = rev(c(rev(magma(300)[1:85]),"white", cividis(300)[210:300])), 
         annotation_row = fx_annotation_4rows, 
         annotation_colors = Annotation_colors,
         annotation_col = Chl_annotation)

ggsave(DegPheatmap, filename = "DEGheatmap_gymnodiniacead1_v2-2019.pdf", width = 7, height = 6)

DEG_to_fx(DEG_table_joined)
# Which Contrast contains the most DEGS? It appears that Day11 contains the most DEGS and is more than 7X Day3, which had the fewewt, more than Day3 has the least (less than half)
dino_table_wcat_clean %>% 
  mutate(Direction = ifelse(LogFC > 0, "pos", "neg")) %>% 
  group_by(Direction) %>% 
  dplyr::count(Contrast) %>% 
  ungroup() %>% 
  mutate(Contrast = factor(Contrast, levels = str_to_lower(sample_list))) %>% 
  ggplot(aes(x = Contrast, y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "log FC",
                    values = jisho_picker('la_lakers'),
                    labels = c("negative", "positive")) +
  labs(y = "Number of DEGs",
       x = "") +
  theme_minimal() 

ggsave(last_plot(), filename = "n_degs_gymnodiacead1-2019.pdf", width = 6, height = 5)
# Phosphorus uptake becoms increasingly important for Pseudo with time
table_wcat_clean %>% 
  filter(Category == "P metabolism") 
# Jitter plot to examine differential espression across time
DegJitter = DEG_table_joined %>% 
  #filter(str_detect(tax, regex("thalas", ignore_case = T))) %>% 
  mutate(Contrast = factor(Contrast, levels = str_to_lower(sample_list)),
         hilite = abs(LogFC)>2) %>% 
  filter(hilite) %>% 
  ggplot(., aes(x = Contrast, y = LogFC)) +
  geom_jitter(color = "black", fill = "#552583", shape = 21, alpha = 0.6) +
  labs(y = "Log Fold Change",
       x = "") +
  theme_minimal()

ggsave(DegJitter, filename = "DegJitter_gymnodiniacea_d1.pdf", width = 6, height = 5)
#Whats the functional enrichment of the different taxa after the death phase?
library(clusterProfiler)
#produce the enrichment table using the enrichKegg function for positive and negative DEGs
Enrichment_DF = DEG_table_joined %>% 
  mutate(Direction = ifelse(LogFC > 0, "positive", "negative")) %>% 
  #filter(LogFC > 0) %>% 
  unite(split, Direction, Contrast, sep = "-", remove = F) %>% 
  split(.$split) %>% 
  map(., ~pull(., KO)) %>% 
  map(., ~enrichKEGG(gene = .,
                     organism = 'ko',
                     pvalueCutoff = 0.05)) %>% 
  map(., ~data.frame(.)) 

names(Enrichment_DF)

Enrich_table = bind_rows("day20-neg" = Enrichment_DF$`negative-day20`,
                         "day18-neg" = Enrichment_DF$`negative-day18`,
                         "day12-neg" = Enrichment_DF$`negative-day12`,
                         "day20-pos" = Enrichment_DF$`positive-day20`,
                         "day18-pos" = Enrichment_DF$`positive-day18`,
                         "day12-pos" = Enrichment_DF$`positive-day12`,
                         .id = "Contrast") %>% 
  dplyr::select(Contrast, Description, GeneRatio, Count, p.adjust) %>%
  separate(Contrast, into = c("day", "dir"), sep = "-") %>% 
  separate(GeneRatio, into = c("num", "denom"), sep = "/") %>% 
  mutate(GeneRatio = as.integer(num)/as.integer(denom)) %>% 
  #filter(GeneRatio > 0.025) %>% 
  dplyr::select(-num, - denom)

Human_diseases = c("Epithelial cell signaling in Helicobacter pylori infection","Human papillomavirus infection","Rheumatoid arthritis","Spinocerebellar ataxia", "Alzheimer disease", "Amyotrophic lateral sclerosis", "Chemical carcinogenesis - reactive oxygen species", "Epstein-Barr virus infection", "Fanconi anemia pathway", "Huntington disease", "Insulin signaling pathway", "Non-alcoholic fatty liver disease", "Parkinson disease", "Prion disease", "Pathways of neurodegeneration - multiple diseases", "Diabetic cardiomyopathy","Vibrio cholerae infection")
Amino_acid_metab = c("Alanine, aspartate and glutamate metabolism", "Cysteine and methionine metabolism", "")
Amino_acid_breakdown = c("Valine, leucine and isoleucine degradation")
Nucleotide_breakdown = c("RNA degradation")
Nucleotide_metabolism = c("Purine metabolism", "Pyrimidine metabolism")
Auto_phagy = c("Autophagy - yeast", "Autophagy - animal")
Fat_metab = c("Fatty acid metabolism", "Fatty acid degradation", "Fatty acid biosynthesis")

Enrich_plot = Enrich_table %>% 
  #distinct(Description) %>% arrange(desc(.))
  #mutate statements are manually aggregating and renaming for ease of interpretation.
  mutate(FxTerm = case_when(str_detect(Description, regex("autoph", ignore_case = T)) ~"Autophagy",
                            str_detect(Description, regex("COVID", ignore_case = T)) ~"Ribosome",
                            str_detect(rownames(.), "(ko00062)|(ko01212)|(ko0414[02])|(ko04966)|(ko05016)|(ko05415)") ~"Lysosome binding and processing",
                            str_detect(rownames(.), "(ko00020)|(ko00220)|(ko00250)") ~"Organic Nuptake and assimilation",
                            str_detect(rownames(.), "(ko012[034]0)") ~"Organic N and P metabolism",
                            str_detect(rownames(.), "(ko04721)|(ko05110)") ~"Oxidative phosphorylation - Lysosome V-type ATPase",
                            Description %in% Human_diseases ~"Human Disease Pathways",
                            Description %in% Nucleotide_metabolism ~"Nucleotide metabolism",
                            Description == "Fatty acid metabolism" ~"Fatty acid metabolism-other",
                            #Description %in% Nucleotide_breakdown ~"Nucleotide breakdown",
                            str_detect(Description, "Cell cycle") ~"Cell cycle",
                            str_detect(Description, regex("(proteolysis)|(proteosome)", ignore_case = T)) ~"Proteolysis")) %>% 
  mutate(FxTerm = ifelse(is.na(FxTerm), Description, FxTerm)) %>% 
  mutate(FxTerm = str_replace(FxTerm, "Biosynthesis of amino acids", "Amino acid biosynthesis"),
         FxTerm = str_replace(FxTerm, "Ribosome$", "Ribosomal protein - Translation"),
         FxTerm = str_replace(FxTerm, "Ribosome biogenesis .+", "Ribosome biogenesis"),
         FxTerm = str_replace(FxTerm, "Meiosis .+", "Meiosis"),
         FxTerm = str_replace(FxTerm, "Proteasome", "Proteolysis"),
         FxTerm = str_replace(FxTerm, "Alanine, aspartate and glutamate metabolism", "GS/GOGAT"),
         day = factor(day, levels = c("day4", "day12", "day18", "day20"))) %>% 
  group_by(day, dir, FxTerm) %>% 
  summarise(GeneRatio = sum(GeneRatio, na.rm = T)) %>% 
  filter(GeneRatio > 0.025) %>% 
  #filter(FxTerm %in% hilight) %>% 
  ggplot(., aes(y = FxTerm, x = day, size = GeneRatio, fill = dir)) +
  geom_point(shape = 21) +
  scale_size_continuous(breaks = c(0.03, 0.14, 0.28, 0.56),
                        labels = c(0.03, 0.14, 0.28, 0.56)) +
  scale_fill_manual(name = "",
                    values = jisho_picker("la_lakers"),
                    labels = c("negative", "positive")) +
  facet_wrap(~dir) +
  labs(x = "",
       y = "") +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = "solid"),
        panel.grid.major.y = element_line(linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))
ggsave(Enrich_plot, filename = "enrichplot_d1-2019.pdf", width = 6, height = 5)

#End