library(clusterProfiler)
library(ggplot2)
library(reshape2)
library(msigdbr)
library(corrgram)
library(RColorBrewer)
library(matrixStats)
setwd("/media/barbitoff/DATA/Working issues/WES/BRK/Pregnancy_loss/pregnancy_loss_genes/")

#GO(BP,MF,CC)
variant_data = read.table('variant_full_list.txt', sep='\t', header=T)
colnames(variant_data)
#Data <- read.csv('unique_gene.txt', header = TRUE, sep = "\t")

# Trimester proportions
table(ceiling(as.numeric(variant_data$Gestational_age_weeks)/13))

genes_plot <- enrichGO(gene = unique(variant_data$Gene),
                       OrgDb = "org.Hs.eg.db",
                       keyType = "SYMBOL",
                       ont           = "BP",
                       #ont           = "CC",
                       #ont           = "MF",
                       pAdjustMethod = "BH",
                       #pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
goplot(genes_plot)
barplot(genes_plot)
write.table(genes_plot,"Genes_GO_BP.csv",quote = F,sep = "\t",row.names = F)


#C2-Pathways,C5-HPO
m_df <- msigdbr(species = "Homo sapiens")
head(m_df) %>% as.data.frame
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO") %>% 
  dplyr::select(gs_name,gene_symbol)
em <- enricher(unique(variant_data$Gene), TERM2GENE=m_t2g)
write.table(em,"Enrich_MSGDIB_C5_HPO_SAB.csv",quote = F,sep = "\t",row.names = F)
barplot(em)

#C2
m_df <- msigdbr(species = "Homo sapiens")
head(m_df) %>% as.data.frame
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name,gene_symbol)
em <- enricher(unique(variant_data$Gene), TERM2GENE=m_t2g)
write.table(em,"Enrich_MSGDIB_C2_SAB.csv",quote = F,sep = "\t",row.names = F)
barplot(em)
heatplot(em, showCategory = 15, foldChange = NULL)+ geom_point() + scale_colour_gradient(low = "blue", high = "red")


#Evolutionary constraint analysis (pLI groups, LOEUF) - Gnomad tables
#pLI groups
Data_PL <- data.frame(gene = unique(variant_data$Gene))
Data_RPL <- data.frame(gene = unique(variant_data[variant_data$Recurrence == 'RPL', ]$Gene))
Data_No_RPL <- data.frame(gene = unique(variant_data[(is.na(variant_data$Recurrence) | 
                                      variant_data$Recurrence == 'NO'), ]$Gene))
Data_Miscarriage <- data.frame(gene = unique(variant_data[variant_data$Outcome == 'Miscarriage', ]$Gene))
Data_Termination <- data.frame(gene = unique(variant_data[variant_data$Outcome == 'Termination', ]$Gene))

Gnomad_table <- read.csv('gnomad.v2.1.1.lof_metrics.by_gene.txt', header = TRUE, sep = "\t")
Gnomad_table$Type <- 'All genes'
Gnomad_table$Category <- 'All'

Gnomad_table_cut <- Gnomad_table[,c(1,21,30,78,79)]
Gnomad_table_cut$Group[Gnomad_table_cut$pLI < 0.1 | Gnomad_table_cut$pLI == 0.1] <- "Tolerant"
Gnomad_table_cut$Group[Gnomad_table_cut$pLI > 0.9 | Gnomad_table_cut$pLI == 0.9] <- "Intolerant"
Gnomad_table_cut$Group[Gnomad_table_cut$pLI > 0.1 & Gnomad_table_cut$pLI < 0.9] <- "Intermediade"


PL_table <- merge(Data_PL, Gnomad_table_cut, by='gene', all.x=T)
PL_table$Type <- 'All PL'
PL_table$Category<- 'All'

RPL_table <- merge(Data_RPL, Gnomad_table_cut, by='gene', all.x=T)
RPL_table$Type <- 'RPL'
RPL_table$Category <- 'Recurrence'

NO_RPL_table <- merge(Data_No_RPL, Gnomad_table_cut, by='gene', all.x=T)
NO_RPL_table$Type <- 'non-RPL'
NO_RPL_table$Category <- 'Recurrence'

Miscarriage_table <- merge(Data_Miscarriage,Gnomad_table_cut, by='gene', all.x=T)
Miscarriage_table$Type <- 'Miscarriage'
Miscarriage_table$Category <- 'Spontaneity'

Termination_table <- merge(Data_Termination,Gnomad_table_cut, by='gene', all.x=T)
Termination_table$Type <- 'Termination'
Termination_table$Category <- 'Spontaneity'

Gnomad_SAB <- rbind(PL_table, RPL_table, NO_RPL_table, 
                    Miscarriage_table, Termination_table, 
                    Gnomad_table_cut)
genes_by_pli <- Gnomad_SAB
genes_by_pli$Group <- factor(genes_by_pli$Group, levels = c("Intolerant", "Intermediade",
                          "Tolerant", ordered = TRUE))
genes_by_pli$Type <- factor(genes_by_pli$Type, levels = c("All genes", "All PL", 
                          "RPL","non-RPL", "Miscarriage", "Termination", ordered = TRUE))

genes_by_pli <- aggregate(Group~Type+Category, data = genes_by_pli, function(x) table(x)/length(x), simplify = FALSE)
genes_by_pli$Tolerant = sapply(genes_by_pli$Group, function(x) x['Tolerant'])
genes_by_pli$Intolerant = sapply(genes_by_pli$Group, function(x) x['Intolerant'])
genes_by_pli$Intermediate = sapply(genes_by_pli$Group, function(x) x['Intermediade'])

pli_props <- melt(genes_by_pli, id.vars=c('Type', 'Category'), 
                  measure.vars = c('Tolerant', 'Intermediate', 'Intolerant'))

ggplot(pli_props, aes(x = Type, y = value, fill = variable)) +
  geom_col(position = "fill")+facet_wrap(vars(Category), scales = "free_x")


#LOEUF
ggplot(Gnomad_SAB, aes(x=Type, y=oe_lof_upper, fill=Type)) +  
  geom_boxplot() + scale_fill_brewer(palette = "Set1") + facet_wrap(vars(Category), scales = "free_x")


# Gene expression level (several metrics)
#a - Median
GTEX_table <- read.csv('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', header = TRUE, sep = "\t")
GTEX_table$Median = rowMedians(as.matrix(GTEX_table[,c(3:56)]))
GTEX_table$TPM10 <- apply(GTEX_table[, 3:56], 1, function(elt) sum(elt > 10))
GTEX_table$Type <- 'All genes'
GTEX_table$Category <- 'All'

PL_table <- merge(PL_table, GTEX_table[, 1:58], by='gene', all.x=T)
RPL_table <- merge(RPL_table, GTEX_table[, 1:58], by="gene", all.x=T)
NO_RPL_table <- merge(NO_RPL_table, GTEX_table[, 1:58], by="gene", all.x=T)
Miscarriage_table <- merge(Miscarriage_table, GTEX_table[, 1:58], by="gene", all.x=T)
Termination_table <- merge(Termination_table, GTEX_table[, 1:58], by="gene", all.x=T)

eplot_cols = c('gene', 'Type', 'Category', 'Median', 'TPM10')
GTEX_SAB <- rbind(GTEX_table[, eplot_cols], PL_table[, eplot_cols], 
                  RPL_table[, eplot_cols], NO_RPL_table[, eplot_cols], 
                  Miscarriage_table[, eplot_cols], Termination_table[, eplot_cols])
GTEX_SAB$Type <- factor(GTEX_SAB$Type, levels = c("All genes", "All PL", "RPL","non-RPL", 
                            "Miscarriage", "Termination", ordered = TRUE))
ggplot(GTEX_SAB, aes(x=Type, y=Median, fill=Type)) +  
  geom_boxplot(outlier.shape = NA) + scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits=c(0, 60))+facet_wrap(vars(Category), scales = "free_x")  

ggplot(GTEX_SAB, aes(x=Type, y=TPM10, fill=Type)) +  
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(vars(Category), scales = "free_x")  


#Prediction
# LOEUF filtering
loeuf_cutoff = median(PL_table$oe_lof_upper, na.rm = T)
#loeuf_cutoff = quantile(PL_table$oe_lof_upper, probs=0.75, na.rm = T)
print(loeuf_cutoff)
Gnomad_Loeuf_filter <- Gnomad_table_cut[Gnomad_table_cut$oe_lof_upper < loeuf_cutoff, ]
Gnomad_Loeuf <- na.omit(Gnomad_Loeuf_filter) 

#Median expression filter
median_exp_cutoff = median(PL_table$Median, na.rm = T)
#median_exp_cutoff = quantile(PL_table$Median, probs=0.25, na.rm = T)
print(median_exp_cutoff)
GTEX_Gnomad_Loeuf <-merge(Gnomad_Loeuf, GTEX_table, by="gene")
Loeuf_Median_filter <- GTEX_Gnomad_Loeuf[GTEX_Gnomad_Loeuf$Median > median_exp_cutoff,]
Loeuf_Median_filter <- Loeuf_Median_filter [, c('gene', 'oe_lof_upper', 'Median')]

#BIOGRID
BIOGRID<- read.csv('BIOGRID-ORGANISM-Homo_sapiens-4.4.225.tab3.txt', 
                   header = TRUE, sep = "\t")
DF_A_B <- BIOGRID[(BIOGRID$Official.Symbol.Interactor.A %in% Data_PL$gene) | 
                    (BIOGRID$Official.Symbol.Interactor.B %in% Data_PL$gene) ,]
Interactors <- unique(c(DF_A_B$Official.Symbol.Interactor.A, 
                        DF_A_B$Official.Symbol.Interactor.B))
Interacting_genes <- Loeuf_Median_filter[Loeuf_Median_filter$gene %in% Interactors, ]

table(Data_PL$gene %in% Interacting_genes$gene)

#C2 pathways filter
C2_all <- read.csv('c2.all.tsv', header = FALSE, sep = "\t")
ALL_SAB_C2 <- read.csv('Enrich_MSGDIB_C2_SAB.csv', header = TRUE, sep = "\t")
ALL_SAB_C2_TOP <- head(ALL_SAB_C2, 10)
C2_top <- C2_all[C2_all$V1 %in% ALL_SAB_C2_TOP$ID, ]
C2_enriched <- C2_all[C2_all$V1 %in% ALL_SAB_C2$ID, ]

UNIQUE_PR_GENES <- unique(unlist(sapply(c(C2_top$V3), function(x) strsplit(x, ';')[[1]])))
Predicted_genes <- Interacting_genes[Interacting_genes$gene %in% UNIQUE_PR_GENES, ]
table(Data_PL$gene %in% Predicted_genes$gene)

PR_GENES_ALL <- unique(unlist(sapply(c(C2_enriched$V3), function(x) strsplit(x, ';')[[1]])))
Broad_Predicted_genes <- Interacting_genes[Interacting_genes$gene %in% PR_GENES_ALL, ]
table(Data_PL$gene %in% Broad_Predicted_genes$gene)


write.table(Predicted_genes, file = 'prediction_results_narrow.tsv', 
            sep='\t', row.names=F)
write.table(Broad_Predicted_genes, file = 'prediction_results_broad.tsv', 
            sep='\t', row.names=F)

# Checking disease associations for the predicted gene list
Monogen <- read.csv('HPO_InheritanceInfo.csv', header = TRUE, sep = "\t")
Monogen_predicted_narrow <- merge(Predicted_genes, Monogen, 
                                  by.x ="gene", by.y ="gene_symbol")
table(sapply(Predicted_genes$gene, function(x) x %in% Monogen_predicted_narrow$gene))

Monogen_predicted_broad <- merge(Broad_Predicted_genes, Monogen, 
                                  by.x ="gene", by.y ="gene_symbol")
table(sapply(Broad_Predicted_genes$gene, 
             function(x) x %in% Monogen_predicted_broad$gene))

