#!/usr/bin/env Rscript

library(maftools)

files <- list.files(path=".", pattern=".oncefiltered.funcotated.vep.maf$")
mafs <- merge_mafs(files)
write.table(mafs@data[,-1], file="COHORT.somatic.maf", sep="\t", row.names=FALSE, na="", quote=FALSE)

laml <-  mafs

getSampleSummary(laml)
getGeneSummary(laml)
getClinicalData(laml)
getFields(laml)

write.mafSummary(maf = laml, basename = 'COHORT.summary')

pdf('COHORT_maf_summary.pdf', paper='a4r', onefile=FALSE)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf('COHORT_oncoplot.pdf', paper='a4r', onefile=FALSE)
oncoplot(maf = laml, top = 10, draw_titv = TRUE)
dev.off()

pdf('COHORT_mutation_load.pdf', paper='a4r', onefile=FALSE)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'COHORT', logscale = TRUE, capture_size = 50)
dev.off()

pdf('COHORT_somatic_interactions.pdf', paper='a4r', onefile=FALSE)
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

pdf('COHORT_oncogenic_pathways.pdf', paper='a4r', onefile=FALSE)
OncogenicPathways(maf = laml)
dev.off()
