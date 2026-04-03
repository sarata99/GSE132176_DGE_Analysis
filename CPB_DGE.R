# set working directory
setwd("~/Documents/smith_lab/CPB_3pod_analysis")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("GEOquery")
install.packages("tidyverse")
library(tidyverse)
BiocManager::install("GEOquery", force = TRUE)
library(GEOquery)
BiocManager::install("oligo")
library(oligo)
library(limma)
install.packages("edgeR")
library(edgeR)

# Download geo dataset
geoData <- getGEO("GSE132176", GSEMatrix = TRUE)
metadata <- pData(geoData[[1]])

# Read files
celFiles <- list.celfiles(path = "./GSE132176_RAW", full.names = TRUE)
rawData <- read.celfiles(celFiles)

# Use RMA to do background subtraction, normalization, and summarization
eset <- rma(rawData)

# Get ExpressionFeatureSet
data <- exprs(eset)

# Lable genes
BiocManager::install("hthgu133pluspm.db", ask = FALSE, update = FALSE)
library(hthgu133pluspm.db)

geneSymbols <- mapIds(hthgu133pluspm.db,
                       keys = rownames(data),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Filter out probes that don't map to any known gene (NAs)
validProbes <- !is.na(geneSymbols)
exprMatrix <- data[validProbes, ]
mappedSymbols <- geneSymbols[validProbes]

# Assign the actual Gene Symbols and collapse duplicates
rownames(exprMatrix) <- mappedSymbols
dataAnnotated <- avereps(exprMatrix, ID = rownames(exprMatrix))

#----------- Determining sex -----------#
# Define sex markers, XIST = female only; KDM5D, RPS4Y1, ElF1AY, NLGN4Y = male only
sexGenes <- c("XIST", "KDM5D", "RPS4Y1", "EIF1AY", "NLGN4Y") 

# Only pull the genes that successfully mapped in annotation step
availableSexGenes <- sexGenes[sexGenes %in% rownames(dataAnnotated)]
maleGenes <- availableSexGenes[availableSexGenes != "XIST"]

sexMatrix <- dataAnnotated[availableSexGenes, ]

# Transpose into a dataframe
sexDf <- as.data.frame(t(sexMatrix)) %>%
  rownames_to_column("sampleID")

# Calculate the mean male signal and apply the threshold
sexDf <- sexDf %>%
  mutate(sampleID = gsub("^(GSM[0-9]+)_.*$", "\\1", basename(sampleID))) %>% # Clean sample ID
  rowwise() %>% # Compute row-by-row for each patient
  mutate(
    meanMaleSignal = mean(c_across(all_of(maleGenes)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    transcriptomicSex = case_when(
      XIST > meanMaleSignal ~ "Female",
      meanMaleSignal > XIST ~ "Male",
      TRUE ~ "Unknown"
    )
  )

# Plot XIST against the Mean Male Signal
ggplot(sexDf, aes(x = XIST, y = meanMaleSignal, color = transcriptomicSex, label = sampleID)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Transcriptomic Sex Stratification",
       x = "XIST Log2 Intensity (Female Marker)",
       y = "Mean Y-Chromosome Log2 Intensity (KDM5D, RPS4Y1, NLGN4Y)") +
  scale_color_manual(values = c("Female" = "#D55E00", "Male" = "#0072B2"))

# Grab just the essential columns and sort them so all Females/Males are grouped together
finalSexList <- sexDf %>%
  dplyr::select(sampleID, transcriptomicSex, XIST, meanMaleSignal) %>%
  arrange(transcriptomicSex, sampleID)
print(finalSexList, n = Inf)

# Merge sex labels into metadata
metadata$sampleID <- rownames(metadata)
metadata <- metadata %>%
  left_join(sexDf %>% dplyr::select(sampleID, transcriptomicSex), 
            by = c("geo_accession" = "sampleID"))

# Create group factor with clean renaming (ex: TOF_Pre_Male, ASD_Post_Female)
metadata <- metadata %>%
  mutate(
    Time = factor(ifelse(grepl("before", `time:ch1`), "Pre", "Post"), levels = c("Pre", "Post")),
    Disease = factor(`disease state:ch1`),
    PatientID = factor(`individual:ch1`),
    Sex = factor(transcriptomicSex)
  )

group <- factor(paste(metadata$Disease, metadata$Time, metadata$Sex, sep = "_"))

# Set up design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# duplicateCorrelation estimates the internal consistency of each child's samples
corfit <- duplicateCorrelation(dataAnnotated, design, block = metadata$PatientID)

fit <- lmFit(dataAnnotated, design, block = metadata$PatientID, correlation = corfit$consensus)

# Create contrasts
contrastMatrix <- makeContrasts(
  TOF_Male_Response = TOF_Post_Male - TOF_Pre_Male,
  TOF_Female_Response = TOF_Post_Female - TOF_Pre_Female,
  ASD_Male_Response = ASD_Post_Male - ASD_Pre_Male,
  ASD_Female_Response = ASD_Post_Female - ASD_Pre_Female,
  levels = design
)

# Final Statistics
fit2 <- contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2, robust = TRUE)

# View Summaries
resultsSummary <- decideTests(fit2)
print(summary(resultsSummary))
vennDiagram(resultsSummary, cex = c(1, 0.6, 1))

# Iterate through results and create a CSV for each entry
results <- list()
for(i in colnames(contrastMatrix)){
  results[[i]] <- topTable(fit2, coef = i, number = Inf) %>% rownames_to_column("Symbol")
}

for(i in names(results)) {
  columns = results[[i]][, c("Symbol","logFC", "P.Value", "adj.P.Val")]
  write.csv(columns, file = paste0(i, ".csv"), row.names = FALSE)
}
