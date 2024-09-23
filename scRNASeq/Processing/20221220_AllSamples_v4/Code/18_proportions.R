# Compositional analysis for single-cell and FACS data

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "18_proportions"
dir.create(outDir)

#### Load packages ####
library(corrplot)
library(foreign)
library(limma)
library(readxl)
library(Seurat)
# Installed the speckle package from https://github.com/phipsonlab/speckle
library(speckle)
# Documentation in the vignette is helpful
# vignette("speckle")
library(tidyr)

#### scRNA-seq propeller analysis: test for difference in proportions #####

# Read the 10X data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"
# Why does line 72-73 in 14_libra_deg_20230426.r state that names(matrices)
# have only 10 different values? Seems to be due to default argument min_cells = 3
# of to_pseudobulk

# Basic analysis
res <- propeller(clusters = seurat$integrated_snn_res.1, sample = seurat$orig.ident, 
          group = seurat$Condition)
res
write.table(res, file=paste0(outDir, "/proportion_analyis_cluster_all.txt"), 
            sep = "\t", col.names = NA)

plotCellTypeProps(clusters = seurat$integrated_snn_res.1, sample = 
                    seurat$orig.ident)

# Only include clusters 0, 1, 4 and 8
seurat.subset <- subset(seurat, subset = integrated_snn_res.1 %in% c(0,1,4,8))
seurat.subset$integrated_snn_res.1 <- factor(seurat.subset$integrated_snn_res.1)
res <- propeller(clusters = seurat.subset$integrated_snn_res.1, sample = seurat.subset$orig.ident, 
          group = seurat.subset$Condition)
res
write.table(res, file=paste0(outDir, "/proportion_analyis_cluster_0148.txt"), sep = "\t", col.names = NA)

plotCellTypeProps(clusters = seurat.subset$integrated_snn_res.1, sample = 
                    seurat.subset$orig.ident)

# Do the same thing but 'manually'
props <- getTransformedProps(seurat.subset$integrated_snn_res.1, seurat.subset$orig.ident, 
                    transform="logit")

# Make a 'targets' table at the sample level
colData <- seurat.subset@meta.data[, c("orig.ident","Condition","subType")]
colData <- unique(colData)
rownames(colData) <- colData$orig.ident
# And put it in the right order
colData <- colData[colnames(props$Counts),]

group <- colData$Condition
design <- model.matrix(~ 0 + group)
colnames(design) <- sub("group", "", colnames(design))
# Results are identical to the basic propeller analysis above
propeller.anova(prop.list=props, design=design, coef = 1:3, 
                robust=TRUE, trend=FALSE, sort=TRUE)

# ANCA vs. NC
contrasts <- c(1,-1, 0)
res <- propeller.ttest(prop.list = props, design=design, contrasts=contrasts, robust=TRUE,
                trend=FALSE, sort=TRUE)
res
write.table(res, file=paste0(outDir, "/proportion_analyis_cluster_0148_ANCA_vs_NC.txt"), 
            sep = "\t", col.names = NA)

# PR3 vs. MPO
group <- factor(paste(colData$Condition, colData$subType, sep = "_"))
group <- relevel(group, ref = "ANCA_PR3")
design <- model.matrix(~ 0 + group)
colnames(design) <- sub("group", "", colnames(design))
propeller.anova(prop.list=props, design=design, coef = 1:4, 
                robust=TRUE, trend=FALSE, sort=TRUE)

contrasts <- c(1,-1, 0, 0)
res <- propeller.ttest(prop.list = props, design=design, contrasts=contrasts, robust=TRUE,
                       trend=FALSE, sort=TRUE)
res
write.table(res, file=paste0(outDir, "/proportion_analyis_cluster_0148_ANCA_PR3_vs_ANCA_MPO.txt"), 
            sep = "\t", col.names = NA)


#### propeller analysis: continuous variable ####

# Read the patient data
patient.data <- read.spss("../../../Data/Dataset_1/Meta/230526_single-cell_patients_voorPerry.sav",
                          to.data.frame = TRUE)
# Gives some warning messages related to long string value/missing labels and 
# duplicated levels 

# Just have a look at some of the variables
colnames(patient.data)
patient.data[,grepl("percentage", colnames(patient.data))]
# Put it in the right order
rownames(patient.data) <- paste0("MOMA", patient.data$incl_subj_number)
patient.data <- patient.data[rownames(colData), ]

# Let's first have a look at an example
# According to Yosta the IFTA score is an indicator of chronic damage
design <- model.matrix(~ 1 + patient.data$IFTA_score)
fit <- lmFit(props$Proportions, design)
fit <- eBayes(fit, robust=TRUE)
topTable(fit,coef=2)

# Plot the regression lines
par(mfrow=c(1,4))
for(i in seq(1,4,1)){
  plot(patient.data$IFTA_score, props$Proportions[i,], main = rownames(props$Proportions)[i], 
       pch=16, cex=2, xlab =  "IFTA_score", ylab="Proportions", cex.lab=1.5, cex.axis=1.5,
       cex.main=2)
  abline(a=fit$coefficients[i,1], b=fit$coefficients[i,2], col=4, 
         lwd=2)
}
# This seems to be in line with Yosta's hypotheses

# However, for this analysis it makes more sense to exclude the NC (MOMA57) and 
# SLE (MOMA68) samples 
index.anca <- colData$Condition == "ANCA"
design <- model.matrix(~ 1 + patient.data[index.anca, ]$IFTA_score)
fit <- lmFit(props$Proportions[, index.anca], design)
fit <- eBayes(fit, robust=TRUE)
topTable(fit,coef=2)

# Generic function for fitting a linear model with a continuous covariate
propeller.continuous <- function(var, index = 1:nrow(patient.data)){
  # var: character string that should correspond to a column name in patient.data
  index <- index[!is.na(patient.data[index, var])]
  design <- model.matrix(~ 1 + patient.data[index, var])
  fit <- lmFit(props$Proportions[ , index], design)
  fit <- eBayes(fit, robust=TRUE)
  fit
}

# Generic function for plotting regression lines per cluster for a linear model 
# with a continuous covariate for two regression models (fit1, fit2)
propeller.regline <- function(var, fit1, fit2, index){
  # var: character string that should correspond to a column name in patient.data
  par(mfrow=c(2,nrow(props$Proportions)))
  tt <- topTable(fit1, coef = 2, sort.by = "none")
  for(i in seq(1,nrow(props$Proportions),1)){
    main.txt <- paste0(rownames(props$Proportions)[i],": p = ", round(tt[i,"P.Value"], 3))
    plot(patient.data[, var], props$Proportions[i,], 
         main = main.txt, pch=16, cex=2, xlab = var, 
         ylab="Proportions", cex.lab=1.5, cex.axis=1.5, cex.main=2)
    abline(a=fit1$coefficients[i,1], b=fit1$coefficients[i,2], col=4, 
           lwd=2)
  }
  tt <- topTable(fit2, coef = 2, sort.by = "none")
  for(i in seq(1,nrow(props$Proportions),1)){
    main.txt <- paste0(rownames(props$Proportions)[i],": p = ", round(tt[i,"P.Value"], 3))
    plot(patient.data[index, var], props$Proportions[i, index], 
         main = main.txt, pch=16, cex=2, xlab = var, 
         ylab="Proportions", cex.lab=1.5, cex.axis=1.5, cex.main=2)
    abline(a=fit2$coefficients[i,1], b=fit2$coefficients[i,2], col=4, 
           lwd=2)
  }
  mtext(var, side = 3, line = -1.3, outer = TRUE)
}
fit1 <- propeller.continuous("IFTA_score")
fit2 <- propeller.continuous("IFTA_score", index.anca)
tt1 <- topTable(fit1, coef=2, sort.by = "none")
tt2 <- topTable(fit2, coef=2, sort.by = "none")
tt1
tt2
propeller.regline("IFTA_score", fit1, fit2, index.anca)

# Make two vectors with the covariates of interest
# See ProjectDocumentation/BackgroundDocumentation/Associaties met histologie.docx
# The ones that are not yet properly coded as a numeric vector have been
# commented out
acute.damage <- c("perc_cellular_cresc_totalglom",
          # "Cellularcresc_score",
          "eGFR",
          "Kreatinine",
          "Monocytes",
          "BVAS",
          "perc_normal_glom",
          "bx_tin")

pdf(file = paste0(outDir,"/acute_damage_reg.pdf"), width = 14)
for (var in acute.damage){
  fit1 <- propeller.continuous(var)
  fit2 <- propeller.continuous(var, index.anca)
  propeller.regline(var, fit1, fit2, index.anca)
}
dev.off()

chronic.damage <- c("IFTA_score",
          "bx_ifta",
          "NIH_chronicity",
          #"Fibrouscresc_score",
          "perc_fibrous_cresc_total",
          #"Glomsclerosis_score",
          "perc_glomscl",
          "perc_fibrosed_or_fibrocellular_cresc",
          "perc_normal_glom"
)

pdf(file = paste0(outDir,"/chronic_damage_reg.pdf"), width = 14)
for (var in chronic.damage){
  fit1 <- propeller.continuous(var)
  fit2 <- propeller.continuous(var, index.anca)  
  propeller.regline(var, fit1, fit2, index.anca)
}
dev.off()

#### FACS data propeller analysis ####
# - New file on 231114

# Import the FACS data and select only relevant columns
#dat <- read_excel("../Data/scRNASeq/Meta/Data-restructured_231013_v3.xlsx", 
#           sheet = 1)[,1:12]
dat <- read_excel("../Data/scRNASeq/Meta/Data_R_restructured_231114.xlsx", 
                  sheet=1)[,1:12]

# Strip the first column to proper patient identifiers
name.split <- sapply(dat$Name, strsplit , " ") 
dat$Name <- sapply(name.split, "[", 2)
dat$Prop_subset <- dat$Perc_subset/100
  
# Pivot to wide format
dat.wide <- pivot_wider(dat, id_cols = c(Name, Groep),
  names_from = Monocyte_subset,
  values_from = c(Prop_subset),
)
# And take the transpose
dat.wide.tr <- t(dat.wide[, -(1:2)])
colnames(dat.wide.tr) <- dat.wide$Name

barplot(dat.wide.tr,col=ggplotColors(nrow(dat.wide.tr)),
        ylab="Monocyte subset proportions", xlab="Samples", las=2)

#### Cross-sectional ####
# Use arcsin square root transformation since it is insensitive to 
# the choice of scale.fac (see vignette, section 13)
prop.list <- convertDataToList(dat.wide.tr,data.type="proportions", 
                               transform="asin", scale.fac = 1)
design <- model.matrix(~0 + dat.wide$Groep)
colnames(design) <- c("predneg","predpos","hc","remission")

res <- propeller.anova(prop.list = prop.list, design = design, coef = 1:3,
                robust=TRUE,trend=FALSE,sort=TRUE)
res
#write.table(res, file=paste0(outDir, "/proportion_analyis_facs_anova.txt"), 
#            sep = "\t", col.names = NA)
write.table(res, file=paste0(outDir, "/proportion_analyis_facs_anova.231114.txt"), 
            sep = "\t", col.names = NA)

# prednison + versus -
mycontr <- makeContrasts(predpos - predneg,levels=design)
res <- propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                robust=TRUE,trend=FALSE,sort=TRUE)
res
# write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predpos_vs_predneg.txt"), 
#             sep = "\t", col.names = NA)
write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predpos_vs_predneg.231114.txt"), 
            sep = "\t", col.names = NA)

# prednison + versus healthy control
mycontr <- makeContrasts(predpos - hc, levels=design)
res <- propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                robust=TRUE,trend=FALSE,sort=TRUE)
res
# write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predpos_vs_hc.txt"), 
#             sep = "\t", col.names = NA)
write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predpos_vs_hc.231114.txt"), 
            sep = "\t", col.names = NA)

# prednison - versus healthy control
mycontr <- makeContrasts(predneg - hc, levels=design)
res <- propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                robust=TRUE,trend=FALSE,sort=TRUE)
res
# write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predneg_vs_hc.txt"), 
#             sep = "\t", col.names = NA)
write.table(res, file=paste0(outDir, "/proportion_analyis_facs_predneg_vs_hc.231114.txt"), 
            sep = "\t", col.names = NA)

#### Longitudinal ####
# Select only the individuals with two measurements
dat.wide$Indiv <- sub("v1|v2","",dat.wide$Name)
dat.wide$Indiv[dat.wide$Indiv=="53"] <- "053"
indx <- dat.wide$Indiv %in% dat.wide$Indiv[duplicated(dat.wide$Indiv)]
dat.wide.ltd <- dat.wide[indx,]
dat.wide.tr.ltd <- dat.wide.tr[,indx]
dat.wide.ltd$Groep <- sub(" \\+ pred| \\- pred","",dat.wide.ltd$Groep)

prop.list <- convertDataToList(dat.wide.tr.ltd,data.type="proportions", 
                               transform="asin", scale.fac = 1)
design <- model.matrix(~0 + dat.wide.ltd$Groep + dat.wide.ltd$Indiv)
colnames(design)[1:2] <- c("active","remission")
colnames(design) <- make.names(colnames(design))

# remission versus active
mycontr <- makeContrasts(remission - active,levels=design)
res <- propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                robust=TRUE,trend=FALSE,sort=TRUE)
res
# write.table(res, file=paste0(outDir, "/proportion_analyis_facs_remission_vs_active.txt"), 
#             sep = "\t", col.names = NA)
write.table(res, file=paste0(outDir, "/proportion_analyis_facs_remission_vs_active.231114.txt"), 
            sep = "\t", col.names = NA)

# remission versus active (without pairing, just for the record)
design <- model.matrix(~0 + dat.wide.ltd$Groep )
colnames(design) <- c("active","remission")

mycontr <- makeContrasts(remission - active,levels=design)
propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                robust=TRUE,trend=FALSE,sort=TRUE)

#### MPO vs. PR3 ####
# mail Yosta 231114 - toevoegen of er een  verschil is tussen MPO en PR3 actieve samples
# Are the samples with missing Serology Healthy Controls?? Yes!!
dat$Serology <- ifelse(is.na(dat$Serology), "HC", dat$Serology)

# But, select only the active samples (whether they are on prednison or not ..)
dat.subset <- dat[grepl("Actief", dat$Groep), ]

# Pivot to wide format
dat.wide.serology <- pivot_wider(dat.subset, id_cols = c(Name, Serology),
                        names_from = Monocyte_subset,
                        values_from = c(Prop_subset),
)
# And take the transpose
dat.wide.serology.tr <- t(dat.wide.serology[, -(1:2)])
colnames(dat.wide.serology.tr) <- dat.wide.serology$Name

prop.list <- convertDataToList(dat.wide.serology.tr, data.type="proportions", 
                               transform="asin", scale.fac = 1)
design <- model.matrix(~0 + dat.wide.serology$Serology)
colnames(design) <- c("MPO","PR3")

mycontr <- makeContrasts(MPO - PR3,levels=design)
res <- propeller.ttest(prop.list = prop.list, design = design, contrasts = mycontr,
                       robust=TRUE,trend=FALSE,sort=TRUE)
res

write.table(res, file=paste0(outDir, "/proportion_analyis_facs_activeSamples_MPO_vs_PR3.231114.txt"), 
            sep = "\t", col.names = NA)
