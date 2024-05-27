#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("maftools")
#install.packages("R.utils")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("NMF")
#BiocManager::install("pheatmap")
#BiocManager::install("barplot3d")

# https://bioconductor.statistik.tu-dortmund.de/packages/3.5/bioc/vignettes/maftools/inst/doc/maftools.html#oncoplots-aka-waterfall-plots

# 27may2024

library("maftools")
library("R.utils")
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
library('NMF')
library("barplot3d")

setwd("processed/")

# this is the functional mutation dataset
combined = read.maf(maf = "cohort_10x_10p_vaf.maf")

# mutational load needs to be on the full dataset
laml.mutload = tcgaCompare(maf = combined, cohortName = 'WES-OSCC', logscale = TRUE, capture_size = 50)

# full previous genes list
genes =c("AJUBA", "ANO1", "APC", "APOBEC2", "ARRDC1", "ATG7", "BRAF", "CCND1", "CCND3", "CCNE1", 
         "CCNF", "CDK6", "CDKN2A", "CHEK2", "CREBBP", "CTNNB1", "CTTN", "CUL3", "DCHS1", "DDX5", 
         "EGFR", "EP300", "EPHA7", "ERBB2", "ERRFI1", "FAT1", "FAT2", "FAT3", "FAT4", "FBXW7", 
         "FGF19", "FGF2", "FGFR1", "FGFR2", "FOXA1", "IGF1R", "JAG1", "KDM5A", "KDM6A", "KEAP1", 
         "KIT", "KMT2C", "KMT2D", "KRAS", "LRP5", "MAP3K13", "MET", "MNT", "MUC2", "MYC", 
         "NACA", "NCOR1", "NFE2L2", "NOTCH1", "NOTCH2", "NOTCH3", "NPRL2", "P63", "PAX9", "PBRM1", 
         "PIK3CA", "PIK3R1", "PTCH1", "PTEN", "RASA3", "RB1", "RUNX2", "SAV1", "SERPINB4", "SMAD2", 
         "SMAD4", "SMARCA4", "SMARCB1", "SOX2", "SOX5", "SOX8", "SOX9", "TGFBR2", "TP53", "TP63", 
         "VEGFA", "ZFP36L2", "ZNF750")

# subset - 10x, 10% vaf: 23 genes
genes_10x =c("ZNF750", "TGFBR2", "RASA3", "NOTCH1", "NOTCH2", "NOTCH3", "LRP5", "JAG1", "FGFR2", "FBXW7", 
             "FAT3", "ERRFI1", "EGFR", "CDKN2A", "NCOR1", "NACA", "FAT4", "FAT2", "FAT1", "KMT2C", 
             "MUC2", "KMT2D", "TP53")

# sample order
sample =c("LF312", "LF367", "LF374", "LF294", "LF320", "LF357", "LF289", "LF293", "LF141", "LF306", 
          "LF175", "LF317", "LF280", "LF415", "LF298")

# pathways
oncoplot(maf = combined, draw_titv = TRUE, top = 50, gene_mar=20,legendFontSize = 1.5,annotationFontSize = 5, anno_height = 2, pathways = 'auto')

# previous genes
oncoplot(maf = combined, draw_titv = TRUE, genes=genes_10x, gene_mar=7,legendFontSize = 1,annotationFontSize = 1, anno_height = 1, fontSize = 0.5, sampleOrder = sample)

# top 22
oncoplot(maf = combined, draw_titv = TRUE, top = 22, gene_mar=7,legendFontSize = 1,annotationFontSize = 1, anno_height = 2, fontSize = 0.5, sampleOrder = sample)

# mutational signatures

#end

getSampleSummary(combined)

getGeneSummary(combined)

getFields(combined)

combined.titv = titv(maf = combined, plot = FALSE, useSyn = TRUE)
plotTiTv(res = combined.titv, sampleOrder = sample)

data.titv = titv(maf = combined, plot = FALSE, useSyn = TRUE)

plotTiTv(res = data.titv)

# vaf plot
vafPlot = plotVaf(maf = combined, vafCol = 'vaf', flip = TRUE)
