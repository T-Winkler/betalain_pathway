library(GenomicRanges)
library(wesanderson)
library(seqinr)
library(ggbio)
library(tidyverse)
library(msa)
library(ape)
library(ggtree)
library(poppr)
library(gggenes)
library(tidyverse)

setwd("~/Documents/Betalain_Pathway/")
# import and format BLAST results (-outfmt 7)
data <- read.table("Output/BLAST_results/BvMYB1_protein.txt", 
                   row.names = NULL, 
                   fill = T)
colnames(data) <- c("query",
                 "scaffold",
                 "identity",
                 "length",
                 "missmatches",
                 "gaps",
                 "start_q",
                 "end_q",
                 "start_s",
                 "end_s",
                 "e_value",
                 "score")
# sort based on score
data_sorted <- data %>%
  arrange(-score)

# different sorting methods, try to find the best

data_sorted <- data %>%
  filter(score > 40) %>%
  group_by(query) %>%
  arrange(start_s)


# keep best hits from each protein and save in different variables
# do only for genes with only one exon
# genes with multiple exons have to be assembled manually
cDOPA5GT <- head(data_sorted, 1)
cDOPA5GT$query <- "cDOPA5GT"
cDOPA5GT$hit <- 1

Betanidin5GT <- head(data_sorted, 8)
Betanidin5GT$query <- "Betanidin5GT"
Betanidin5GT$hit <- c(1:8)

Betanidin6GT <- head(data_sorted, 2)
Betanidin6GT$query <- "Betanidin6GT"
Betanidin6GT$hit <- c(1:2)
# keep only the best exon to ease visualisation
#CYP76AD2 <- head(data_sorted, 2)
#CYP76AD2$query <- "CYP76AD2"
#CYP76AD2_edit <- CYP76AD2 %>%
  #filter(score > 1200)
#DODA1 <- head(data_sorted, 3)
#DODA1$query <-"DODA1"  
#DODA1_edit <- DODA1 %>%
  #filter(score > 700)
#DODA2 <- head(data_sorted, 3)
#DODA2$query <- "DODA2"
#DODA2_edit <- DODA2 %>%
  #filter(score > 800)

# Load in manually corrected file for genes with exons
man_exons <- read.csv2("Output/genes_manual_corrected_with_MYB1.csv")

man_exons2 <- man_exons %>%
  filter(query %in% c("CYP76AD","DODA1","DODA2","BvMYB1-like"))

man_exons2$hit <- c(1,2,3,2,3,1,4,3,4,2,1,1,2) 

# combine all dataframes and save as file
# Make sure start is alwas smaller than end
# this does not take into account the gene orientation but produces
# error messages otherwise during intersect step
combined <- rbind(cDOPA5GT,
                  Betanidin5GT,
                  Betanidin6GT,
                  man_exons2)  
temp_min <- pmin(combined$start_s, combined$end_s)
combined$end_s <- pmax(combined$start_s, combined$end_s)
combined$start_s <- temp_min
#write.csv(combined, "Data/combined_best_exons.csv", row.names = T, quote = F)
#prepare BED-file output to match with gff (obtain gene names)
bed <- combined %>%
  select(query,
         scaffold,
         start_s,
         end_s)
bed <- bed[, c(2,3,4,1)]
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name")
# start is supposed to always be smaller than end
temp_min <- pmin(bed$chromStart, bed$chromEnd)
bed$chromEnd <- pmax(bed$chromStart, bed$chromEnd)
bed$chromStart <- temp_min
write.table(bed, "Output/betalain_genes.bed", quote=F, row.names = F, col.names = F, sep="\t")

# use bedtools intersect -a bed -b gff -wa -wb
# read in intersect outcome and join with combined table to obtain gene names and gene positions
intersect <- read.table("Output/intersect.txt")
intersect <- intersect %>%
  filter(V7=="gene")
intersect <- intersect[ -c(4:7, 10:12)]
# keep only the gene name
intersect <- intersect %>%
  mutate(V13 = substr(intersect$V13, 4, 11))
# similar colnames ease joining of tables
colnames(intersect) <- c("scaffold","start_s","end_s","start_a","end_a","gene")

# join both the combined table with the blast results and the prepared intersect table
joined_combined <- left_join(combined, intersect)
# save results, table of all identified genes
write.table(joined_combined, "Output/gene_table.txt", quote=F, row.names = F)
# add locus column for plotting
joined_combined$locus <- c("chr4_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_1",
                           "chr1_2",
                           "chr1_2",
                           "chr3_1",
                           "chr3_1",
                           "chr16_2",
                           "chr1_3",
                           "chr1_4",
                           "chr16_2",
                           "chr6_1",
                           "chr6_1",
                           "chr6_1",
                           "chr16_2",
                           "chr6_1",
                           "chr6_1",
                           "chr6_1",
                           "chr10_1",
                           "chr16_1")
# set levels for plot order
joined_combined$locus <- factor(joined_combined$locus, levels = c("chr1_1",
                                                                  "chr1_2",
                                                                  "chr1_3",
                                                                  "chr1_4",
                                                                  "chr3_1",
                                                                  "chr4_1",
                                                                  "chr6_1",
                                                                  "chr10_1",
                                                                  "chr16_1",
                                                                  "chr16_2"))


# plot with facet_wrap, think about gene orientation                           
ggplot(data=joined_combined, 
       aes(xmin = start_s, 
           xmax = end_s, 
           y = query,
           fill = query)) +
  geom_gene_arrow(aes(alpha=1/hit)) +
  facet_wrap(~ locus, 
             scales = "free", 
             ncol=1) +
  guides(alpha=F) +
  labs(y = "Locus", fill = "Pathway gene")




