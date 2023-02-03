################################################################################
## Data import and compositional data transformations
################################################################################
library(zCompositions)
library(compositions)
library(dplyr)
library(phyloseq)
library(ggplot2)

# Read Archaea ASV table (noRares)
ASVarch <- read.csv("./Data/ArchaeaASVnoRares.csv")

# Assign ASV column to row names
rownames(ASVarch) <- as.character(unlist(ASVarch[, 1]))
ASVarch <- ASVarch[, -1]
# Imputation of zero values using Geometric Bayesian multiplicative replacement
arch.asv.ZeroRepl <- cmultRepl(t(ASVarch), label = 0, method = "GBM", 
                               output = "p-counts")
# Centred log ratio transformation
arch.asv.clrTrans <- apply(t(arch.asv.ZeroRepl), 2, function(x) {
  log(x) - mean(log(x))
})
# Transpose for WGCNA
arch.asv.clean = t(arch.asv.clrTrans)

# Import taxonomy table and create phyloseq objects 
ARCH_ASV_tax <- read.csv("./Data/ArchaeaASVtaxonomy.csv") 
rownames(ARCH_ASV_tax) <- as.character(unlist(ARCH_ASV_tax[, 1]))
ARCH_ASV_tax <- as.matrix(ARCH_ASV_tax[, -1])

# Combine ARCH 16S data into phyloseq object
ASVarch = otu_table(arch.asv.clrTrans, taxa_are_rows = TRUE)
TAXarch = tax_table(ARCH_ASV_tax)
ARCHseq_asv = phyloseq(ASVarch, TAXarch)

# Subset out Thaumarchaeal (AOA) reads
nspmls_asv <- subset_taxa(ARCHseq_asv, Phylum=="Thaumarchaeota")
nspmls_asv

# Covert back to data matrix and transpose for MixOmics
nspmls.asv.clean = t(as(otu_table(nspmls_asv), "matrix"))

# Read Bacteria ASV table (noRares)
ASVbact <- read.csv("./Data/BacteriaASVnoRares.csv")

# Assign ASV column to row names
rownames(ASVbact) <- as.character(unlist(ASVbact[, 1]))
ASVbact <- ASVbact[, -1]
# Imputation of zero values using Geometric Bayesian multiplicative replacement
bact.asv.ZeroRepl <- cmultRepl(t(ASVbact), label = 0, method = "GBM", 
                               output = "p-counts")
# Centred log ratio transformation
bact.asv.clrTrans <- apply(t(bact.asv.ZeroRepl), 2, function(x) {
  log(x) - mean(log(x))
})
# Transpose for WGCNA
bact.asv.clean = t(bact.asv.clrTrans)

# ASV-level Clustering 
BACT_ASV_tax <- read.csv("./Data/BacteriaASVtaxonomy.csv") 
rownames(BACT_ASV_tax) <- as.character(unlist(BACT_ASV_tax[, 1]))
BACT_ASV_tax <- as.matrix(BACT_ASV_tax[, -1])

# Combine BACT 16S data into phyloseq object
ASVbact = otu_table(bact.asv.clrTrans, taxa_are_rows = TRUE)
TAXbact = tax_table(BACT_ASV_tax)
BACTseq_asv = phyloseq(ASVbact, TAXbact)
BACTseq_asv <- subset_taxa(BACTseq_asv, Order!="Chloroplast")

bact.asv.clean <- t(as(otu_table(BACTseq_asv), "matrix"))

# Select all SUP05, Nitrospina (NOB) and SAR11 reads
bactSelect_asv <- subset_taxa(BACTseq_asv, Genus=="SUP05_cluster" |
                                Family=="Nitrospinaceae" | Order=="SAR11_clade")

select.asv.clean = t(as(otu_table(bactSelect_asv), "matrix"))

# Read Saanich Inlet Metadata
## Samples MUST be in the same order as in your OTU/ASV tables
saanich <- read.csv("./Data/SaanichN2OMetadata.csv", 
                    header = TRUE, sep = ",")
rownames(saanich) <- saanich$Sample_ID

# Data column for combined nitrate/nitrite
saanich$NO3_NO2 <- saanich$NO3Mol + saanich$NO2Mol 
# Calculate molar N2O yield
saanich$N2O_yield <- (saanich$NH4_Prod) / (saanich$Nitrification+1e-05)*100 

# Select relevant subset of variables for PLSR and WGCNA
saanich.trim <- saanich %>% 
  dplyr::select(O2Mol, delta_N2O, NH4Mol, NO3_NO2, NH4_Prod, NO3_Prod, 
                Nitrification, N2O_yield)

# Subset Bacteria / Metadata to match samples with Archaeal reads (n=18)
common <- intersect(rownames(saanich.trim), rownames(nspmls.asv.clean))
saanich.trim.common <- subset(saanich.trim, rownames(saanich.trim) %in% common)

nspmls.asv.clean.common <- subset(nspmls.asv.clean, rownames(nspmls.asv.clean) 
                                %in% common)
bact.asv.clean.common <- subset(bact.asv.clean, rownames(bact.asv.clean) 
                                %in% common)
select.asv.clean.common <- subset(select.asv.clean, rownames(select.asv.clean) 
                                %in% common)

# Combine archaea and bacteria datasets
saanich.select.asv.comb <- cbind(nspmls.asv.clean, select.asv.clean.common)
saanich.asv.comb <- cbind(arch.asv.clean, bact.asv.clean.common)

################################################################################
## Weighted Gene Correlational network analyses
################################################################################
library(WGCNA)

data = saanich.asv.comb

# Take a quick look at what is in the data set:
dim(data);
names(data);

# Assign data to new working variable to streamline code
datExpr0 = as.data.frame(data);

# Check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# Cluster samples to check for obvious outliers.
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# No outliers apparent so no removal steps necessary
datExpr = datExpr0
# Assign trimmed metadata to sample traits dataframe
datTraits = saanich.trim.common
datTraits.scale <- as.data.frame(scale(datTraits, center = TRUE, scale = TRUE))

dim(datTraits.scale)
names(datTraits.scale)

# We now have the expression data in the variable datExpr, and the scaled 
# environmental traits in the variable datTraits. Before we continue with 
# network construction and module detection, we visualize how the sample traits
# relate to the sample dendrogram.

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, 
# grey means missing entry
traitColors = numbers2colors(datTraits.scale, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits.scale),
                    main = "Sample dendrogram and trait heatmap")  
# In the plot, white means a low value, red a high value, and 
# grey a missing entry.The last step is to save the relevant expression and 
# trait data for use in the next steps of the tutorial.
save(datExpr, datTraits.scale, file = "./WGCNA_dataInput.RData")

# Load the data saved in the first part
lnames = load(file = "WGCNA_dataInput.RData");
# The variable lnames contains the names of loaded variables.

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Set the soft power threshold
softPower = 8;
adjacency = adjacency(datExpr, power = softPower, type = "signed");

# Turn adjacency into topological overlap 
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#  Set the minimum module size relatively high: 20 taxa
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Set similarity threshold for combining modules
MEDissThres = 0.30
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, 
                          verbose = 3)
# The merged module colors
mergedColors = merge$colors;
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNA-networkConstruction-auto.RData")

# Load network data saved in the second part.
lnames = load(file = "WGCNA-networkConstruction-auto.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Extract the taxa-specific correlation coefficients for each sample trait
# Isolate correlation coefficients for Nitrification
trait.ammox = as.data.frame(datTraits$Nitrification);
names(trait.ammox) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), 
                                          nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.ammox, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), 
                                          nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.ammox), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.ammox), sep="");

module.member.traits.ammox <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                    GSPvalue, moduleColors)

# Isolate coefficients for NO3->N2O
trait.reduc1 = as.data.frame(datTraits$NO3_Prod);
names(trait.reduc1) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.reduc1, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.reduc1), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.reduc1), sep="");

module.member.traits.reduc1 <- cbind(taxaModuleMembership, taxaTraitSignificance,
                                     GSPvalue, moduleColors)

# Isolate coefficients for NH4->N2O
trait.oxProd = as.data.frame(datTraits$NH4_Prod);
names(trait.oxProd) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.oxProd, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.oxProd), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.oxProd), sep="");

module.member.traits.oxProd <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                     GSPvalue, moduleColors)

# Isolate coefficients for N2O yields
trait.yield = as.data.frame(datTraits$N2O_yield);
names(trait.yield) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.yield, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.yield), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.yield), sep="");

module.member.traits.yield <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                    GSPvalue, moduleColors)

# Isolate coefficients for N2O yields
trait.delta = as.data.frame(datTraits$delta_N2O);
names(trait.delta) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.delta, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.delta), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.delta), sep="");

module.member.traits.delta <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                    GSPvalue, moduleColors)

# Now we are going to construct a singal data frame that contains correlation
# coefficients, p-values, connectivity measures and taxonomic information for 
# each ASV contained in the WGCNA analyses
Alldegrees1=intramodularConnectivity(adjacency, mergedColors)
head(Alldegrees1)

connectivity <- cbind(Alldegrees1, mergedColors,taxaModuleMembership)
connectivity <- tibble::rownames_to_column(connectivity, "ASV")
connectivity <- connectivity %>% 
  rename(Subnetwork = mergedColors)
# Rename moduleColors to Subnetworks
connectivity$Subnetwork[connectivity$Subnetwork == "brown"] <- "SNET1"
connectivity$Subnetwork[connectivity$Subnetwork == "blue"] <- "SNET2"
connectivity$Subnetwork[connectivity$Subnetwork == "turquoise"] <- "SNET3"

connectivity.traits <- cbind(connectivity, 
                             module.member.traits.ammox$GS.Trait,
                             module.member.traits.ammox$p.GS.Trait,
                             module.member.traits.delta$GS.Trait,
                             module.member.traits.delta$p.GS.Trait,
                             module.member.traits.reduc1$GS.Trait,
                             module.member.traits.reduc1$p.GS.Trait,
                             module.member.traits.oxProd$GS.Trait,
                             module.member.traits.oxProd$p.GS.Trait,
                             module.member.traits.yield$GS.Trait,
                             module.member.traits.yield$p.GS.Trait)

connectivity.traits <- connectivity.traits %>% 
  dplyr::rename(
    "GS.Ammox" = "module.member.traits.ammox$GS.Trait",
    "p.Ammox" = "module.member.traits.ammox$p.GS.Trait",
    "GS.Delta" =  "module.member.traits.delta$GS.Trait",
    "p.Delta" =  "module.member.traits.delta$p.GS.Trait",
    "GS.Reduc" = "module.member.traits.reduc1$GS.Trait",
    "p.Reduc" = "module.member.traits.reduc1$p.GS.Trait",
    "GS.OxProd" = "module.member.traits.oxProd$GS.Trait",
    "p.OxProd" = "module.member.traits.oxProd$p.GS.Trait",
    "GS.Yield" = "module.member.traits.yield$GS.Trait",
    "p.Yield" = "module.member.traits.yield$p.GS.Trait"
    )

# Include taxonomy information in data frame
# Match datasets by rows
tax <- as.data.frame(rbind(BACT_ASV_tax, ARCH_ASV_tax))
tax <- tibble::rownames_to_column(tax, "ASV")

## Create Data Vector of taxonomies for heatmap labelling
tax <-rbind(ARCH_ASV_tax, BACT_ASV_tax)
tax <- tibble::rownames_to_column(as.data.frame(tax), "ASV")

tax.mod <- tax %>% 
  mutate(Tax = ifelse(grepl("SUP05_cluster", Genus), "SUP05", Genus))
tax.mod <- tax.mod %>% 
  mutate(Tax = ifelse(grepl("uncultured", Genus), Family, Tax))
tax.mod <- tax.mod %>% 
  mutate(Tax = ifelse(grepl("uncultured", Order), Class, Tax))

tax.mod$Tax <- gsub("_unclassified","",as.character(tax.mod$Tax))
tax.mod$Tax <- gsub("__ge","",as.character(tax.mod$Tax))
tax.mod$Tax <- gsub("_ge","",as.character(tax.mod$Tax))
tax.mod$Tax <- gsub("Candidatus_","",as.character(tax.mod$Tax))
tax.mod$Tax <- gsub("marine_group","MG",as.character(tax.mod$Tax))

connectivity.traits.tax <- merge(connectivity.traits, tax.mod[, c("ASV", "Phylum", "Class","Order",
        "Family", "Genus", "Tax")], by="ASV")
rownames(connectivity.traits.tax) <- as.character(unlist(connectivity.traits.tax[, 1]))
#write.csv(connectivity.traits.tax, "TaxaTraitRelationships.csv")

# GGPLOTs for Connectivity figures (Fig.3)
snet1.col <- RColorBrewer::brewer.pal(1, "Dark2")[1]
snet2.col <- RColorBrewer::brewer.pal(2, "Dark2")[2]
snet3.col <- RColorBrewer::brewer.pal(3, "Dark2")[3]

# Select all SUP05 and NOB reads
sup05 <- subset_taxa(BACTseq_asv, Genus=="SUP05_cluster")
sup05 = as(otu_table(sup05), "matrix")

sar11 <- subset_taxa(BACTseq_asv, Order=="SAR11_clade")
sar11 = as(otu_table(sar11), "matrix")

names1 = rownames(sup05)
names2 = rownames(sar11)
names3 = colnames(nspmls.asv.clean)

# Create connectivity plots
# summary(lm(GS.Ammox ~ MMblue, data=connectivity.traits.tax))
plot1 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMblue, GS.Ammox), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMblue, GS.Ammox, colour=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(.4, 4),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (Nitrification)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(face = "plain"),
        legend.position = c(0.775,0.175), legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm")) 
plot1 <- plot1 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, GS.Ammox), shape = 0, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, GS.Ammox), shape = 2, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMblue, GS.Ammox), shape = 4, size = 2.5) +
  guides(color="legend", size = "none")
plot1


plot2 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMbrown, GS.Reduc), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMbrown, GS.Reduc, colour=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  labs(y = expression("ASV importance (N"*O[3]^"-"%->%N[2]*"O)"), x = NULL) +
  scale_size_continuous(range = c(.4, 4),
                        name=expression(bold("K"["in"]))) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(face = "bold"),
        legend.position = "none")
plot2 <- plot2 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMbrown, GS.Reduc), shape = 0, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMbrown, GS.Reduc), shape = 2, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMbrown, GS.Reduc), shape = 4, size = 2.5) 
plot2

# summary(lm(GS.OxProd ~ MMbrown, data=connectivity.traits.tax))
plot3 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMbrown, GS.OxProd), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMbrown, GS.OxProd, colour=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(.4, 4),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (N"*H[4]^"+"%->%N[2]*"O)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(face = "bold"),
        legend.position = "none")
plot3 <- plot3 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMbrown, GS.OxProd), shape = 0, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMbrown, GS.OxProd), shape = 2, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMbrown, GS.OxProd), shape = 4, size = 2.5) 
plot3

# summary(lm(GS.Yield ~ MMbrown, data=connectivity.traits.tax))
plot4 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMbrown, GS.Yield), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMbrown, GS.Yield, colour=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(.4, 4),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance ("*N[2]*"O"~"yield)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(face = "bold"), 
        legend.position = "none")
plot4 <- plot4 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMbrown, GS.Yield), shape = 0, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMbrown, GS.Yield), shape = 2, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMbrown, GS.Yield), shape = 4, size = 2.5) 
plot4


plot5 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMblue, GS.Delta), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMblue, GS.Delta, colour=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(.4, 4), 
                        name=expression("K"["in"])) +
  labs(y = expression("ASV importance ("*Delta*N[2]*"O)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(face = "plain"),
        legend.position = c(0.875,0.25), legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm"))  
plot5 <- plot5 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, GS.Delta), shape = 0, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, GS.Delta), shape = 2, size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMblue, GS.Delta), shape = 4, size = 2.5) +
  guides(size = "legend", color="none")
plot5

plot6 <- ggplot(connectivity.traits.tax, aes(MMblue, MMturquoise)) + 
  geom_smooth(aes(MMblue, MMturquoise), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(color=Subnetwork,
                 size=kWithin), alpha = 0.45) +
  scale_size_continuous(range = c(.4, 4), 
                        name=expression(bold("K"["in"]))) +
  labs(x=NULL, y="SNET3 Membership") +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2",
                     direction = -1) +
  theme_bw() 

plot6 <- plot6 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, MMturquoise, shape = "SUP05"), size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, MMturquoise, shape = "SAR11"),  size = 2.5) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(MMblue, MMturquoise, shape = "AOA"), size = 2.5) +
  scale_shape_manual(name = "Taxa", values = c('SUP05' = 0, 'SAR11' = 2, 'AOA' = 4)) +
  theme(panel.grid = element_blank(), legend.position = c(0.83,0.16),
        legend.background = element_rect(color="black"),
        legend.key.size = unit(1, "mm"), legend.title = element_blank()) +
  guides(shape = "legend", size = "none", color="none")
plot6


rates <- cowplot::plot_grid(plot3 + theme(axis.title.x = element_blank()),
                           plot5 + theme(axis.title.x = element_blank()),
                           plot4 + theme(axis.title.x = element_blank()), 
                           plot1 + theme(axis.title.x = element_blank()),
                           plot2 + scale_x_continuous(name = "SNET1 Membership", 
                                                      position = "bottom"),
                           plot6 + scale_x_continuous(name = "SNET2 Membership", 
                                                        position = "bottom"),
                           labels= c("a)", "b)", "c)", "d)", "e)", "f)"),
                           label_fontface="plain",
                           label_x = 0.2, label_y = 0.955,
                           heights = c(1,1,1.05), 
                           common.legend = TRUE,
                           nrow=3, ncol=2)
rates

#ggsave(rates, filename = "Jameson_et_al_Fig4.pdf",
#        height = 7.87, width = 6.06, dpi = 300)

################################################################################
## Propr Network analyses
################################################################################
library(propr)

arch.asv.ZeroRepl.pr <- subset(arch.asv.ZeroRepl, rownames(arch.asv.ZeroRepl) %in% common)
bact.asv.ZeroRepl.pr <- subset(bact.asv.ZeroRepl, rownames(bact.asv.ZeroRepl) %in% common)

arch.bact.clr <- cbind(arch.asv.ZeroRepl.pr, bact.asv.ZeroRepl.pr)

arch.bact.pr <- propr(arch.bact.clr, # rows as samples, like it should be
                      metric = "rho", # or "phi", "phs", "cor", "vlr"
                      ivar = 'clr', # or can use "iqlr" instead
                      alpha = NA, # use to handle zeros
                      p = 1000) # used by updateCutoffs

NoRaresCorrALL <- getResults(arch.bact.pr)

# Match datasets by rows
test <- tibble::rownames_to_column(connectivity.traits.tax, "Partner")

tax.source <-  tax.mod %>% dplyr::rename(Partner = ASV)
tax.target <- tax.mod %>% dplyr::rename(Pair = ASV)

NoRaresCorrALL <- merge(NoRaresCorrALL, tax.source[, c("Partner", "Tax")], by="Partner")
NoRaresCorrALL <- merge(NoRaresCorrALL, tax.target[, c("Pair", "Tax")], by="Pair")

NoRaresCorrALL <- merge(NoRaresCorrALL, test[, c("Partner", "kWithin")], by="Partner")
NoRaresCorrALL <- merge(NoRaresCorrALL, test[, c("Partner", "Subnetwork")], by="Partner")

# Subset network to include rho values > |0.60|
NoRaresCorrALL.subset <- subset(NoRaresCorrALL, propr >= 0.60 | propr <= -0.60)

#write.csv(NoRaresCorrALL.subset, "arch.bact.propr.network.60.csv")

# Import outputs of Cytoscape network analysis
node <- read.csv("./Data/ProprNetworkAnalysis.csv")

matchtax <- c("ARCH28", "BACT54", "BACT121", "ARCH37", "ARCH40", "BACT260", "BACT1644", "ARCH4", "BACT211",
              "BACT110", "BACT5", "BACT19", "BACT94", "BACT52", "BACT15", "BACT10", "BACT2", "BACT25")
splsr.taxa<-node[match(matchtax, node$ASV, nomatch=0),]


# Create important taxa layers
SAR11 <- node %>% subset(Tax == "SAR11")
SUP05 <- node %>% subset(Tax == "SUP05")
Rhodo <- node %>% subset(Tax == "Rhodobacteraceae" | Tax == "Amylibacter" |
                           Tax == "Roseibacillus" | Tax == "Planktomarina" |
                           Tax == "Tateyamaria")
verruco <- node %>% subset(Tax == "Verrucomicrobiales" | Tax == "Verrucomicrobiae" | Tax == "MB11C04_MG")
flavo <- node %>% subset(Tax == "Flavobacteriaceae" | Tax == "Flavobacteriales")
marini <- node %>% subset(Tax == "Marinimicrobia")
desulfo <- node %>% subset(Tax == "Desulfobacteraceae")
ecto <- node %>% subset(Tax == "Ectothiorhodospiraceae")
bacter <- node %>% subset(Tax == "Bacteroidales" | Tax == "Bacteroidetes_BD2-2")
nitro <- node %>% subset(Tax == "Nitrospina")
pelag <- node %>% subset(Tax == "Nitrosopelagicus")
pumil <- node %>% subset(Tax == "Nitrosopumilus")
greys <- node %>% subset(Tax == "Synechococcus_CC9902")

keystone <- ggplot() + 
  geom_point(data=node,aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour="darkgrey", alpha =0.75) +
  geom_point(data=marini, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#8C510A", alpha =0.75) +
  geom_point(data=desulfo, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#FBF500", alpha =0.75) +
  geom_point(data=ecto, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#00E3A1", alpha =0.75) +
  geom_point(data=bacter, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#0000FF", alpha =0.75) +
  geom_point(data=nitro, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#B3DE69", alpha =0.75) +
  geom_point(data=pelag, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#0B3B3B", alpha =0.75) +
  geom_point(data=pumil, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#CB181D", alpha =0.75) +
  geom_point(data=SAR11, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#2B8CBE", alpha =0.75) +
  geom_point(data=Rhodo, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#E6AB02", alpha =0.75) +
  geom_point(data=SUP05, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour = "#88419D", alpha =0.75) +
  geom_point(data=greys, aes(Degree, ClosenessCentrality, size = BetweennessCentrality, 
        shape=Subnetwork), colour="darkgrey") +
  geom_point(data=splsr.taxa, aes(Degree, ClosenessCentrality), 
        colour = "black", size=2, shape=8) +
  scale_size_continuous(range = c(1, 6), name="Betweenness centrality") +
  scale_shape_manual(values = c(15,17,16)) +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  labs(x="Node degree", y="Closeness centrality") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = c(0.60,0.2),
        legend.key.size = unit(1, 'lines'), legend.title = element_text(size=10),
        legend.text = element_text(size=8), legend.box = "horizontal")
keystone 

#ggsave(rates, filename = "Jameson_et_al_Fig3b.pdf",
#        height = 4, width = 4, dpi = 300)
################################################################################
## Partial Least Squares Regression
################################################################################
library(mixOmics)

# PLSR with 5 latent components on the full microbial ASV table
spls.asv.env <- spls(saanich.asv.comb, saanich.trim.common, 
                     mode = "regression", ncomp = 5)

# Evaluate performance of model components
perf.pls = perf(spls.asv.env, validation="loo", folds = 18, 
                progressBar = FALSE, nrepeat = 1)
perf.pls$measures$Q2.total
perf.pls$measures$R2
perf.pls$measures$RMSEP
plot(perf.pls, criterion = 'Q2.total')
plot(perf.pls, criterion = 'R2')
plot(perf.pls, criterion = 'RMSEP')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(15, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(2:8) 


# Perform variable selection on X and Y dataframe 
tune.spls<- tune.spls(saanich.asv.comb, saanich.trim.common, ncomp = 2,
                      test.keepX = list.keepX,
                      test.keepY = list.keepY,
                      nrepeat = 1, folds = 18,
                      validation="loo",
                      mode = 'regression', measure = 'cor') 
plot(tune.spls)

optimal.keepX <- tune.spls$choice.keepX # extract optimal number of X variables
optimal.keepY <- tune.spls$choice.keepY # extract optimal number of Y variables
optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

# Re-fit PLSR model using tuned parameters
spls.asv.env <- spls(saanich.asv.comb, saanich.trim.common,
                     mode = "regression", ncomp = optimal.ncomp,
                     keepX = optimal.keepX, keepY = optimal.keepY)

## Create Data Vector of taxonomies for heatmap labelling
taxa <- as.data.frame(spls.asv.env$names$colnames$X)
names(taxa)[1] <- "ASV"

taxmatch <- merge(taxa, tax.mod[, c("ASV", "Tax")], by="ASV")
taxmatch.reorder <- taxmatch[match(taxa$ASV, taxmatch$ASV), ]

head(saanich.trim.common)
# Clean up text for sample trait names
sample.traits <- c(expression(italic("In situ")~O[2]),
                   expression(Delta*"N"[2]*"O"),
                   expression("[ NH"[4]^"+"*" ]"),
                   expression("[ NO"[3]^"-" * "+" * "NO"[2]^"-"*" ]"),
                   expression("N"*H[4]^"+"%->%N[2]*"O"),
                   expression("N"*O[3]^"-"%->%N[2]*"O"),
                   expression("Nitrification"),
                   expression("N"[2]*"O yield")
)

# Create taxonomy table mapped to spls object
a<-as.data.frame(BACTseq_asv@tax_table)
b<-as.data.frame(ARCHseq_asv@tax_table)
t<-rbind(a,b)
t <- tibble::rownames_to_column(as.data.frame(t), "ASV")
t <- merge(t, connectivity.traits[, c("ASV", "Subnetwork")], by="ASV")
t <- t[match(taxa$ASV, t$ASV), ]


# Assign color key to distinguish bacteria from archaea
tax.col = color.mixo(as.factor(t$Subnetwork))
unique(tax.col)
tax.col<- replace(tax.col, tax.col=="#388ECC", "#7570b3")
tax.col<- replace(tax.col, tax.col=="#F68B33", "#d95f02")
tax.col<- replace(tax.col, tax.col=="#C2C2C2", "#1b9e77")

dev.off()
# Clustered coefficient heatmap (Fig. 1)
cim.select <- cim(spls.asv.env, 
                  row.names = taxmatch.reorder$Tax, 
                  symkey = TRUE,
                  col.names = sample.traits, 
                  cutoff = 0.2,
                  row.sideColors = tax.col,
                  row.cex = 0.60,
                  legend=list(legend = c("SNET1", "SNET2", "SNET3"), 
                              col = c("#7570b3", "#d95f02", "#1b9e77"),
                              cex=0.75, vjust=1), keysize = c(1,1),
                  margins = c(6,6),
                  save = 'pdf'
)

# Print pairwise correlation matrix
plsr.correlations <- cim.select$mat
plsr.correlations

plsr.correlations <- tibble::rownames_to_column(as.data.frame(plsr.correlations), "ASV")
plsr.taxmatch <- merge(plsr.correlations, saanich.tax.comb[, 
                       c("ASV","Phylum", "Class","Order", "Family", "Genus")], 
                       by="ASV")

#write.csv(plsr.taxmatch, "plsr.correlations.csv")
#
################################################################################
## PICRUSt2 nosZ predictions
################################################################################

# Read bacterial gene predictions
GENEbact <- read.csv("./PICRUSt2/BacteriaKEGGabund.csv")
rownames(GENEbact) <- as.character(unlist(GENEbact[, 1]))# Change header names to top row names
GENEbact <- GENEbact[, -1]

gene.path.ZeroRepl <- cmultRepl(t(GENEbact), label = 0, method = "GBM", output = "p-counts")

gene.path.clrTrans <- apply(t(gene.path.ZeroRepl), 2, function(x) {
  log(x) - mean(log(x))
})

gene.path.clrTrans.trans <- t(gene.path.clrTrans)


# Combine gene abundances and metadata
GenePathComp <-cbind(gene.path.clrTrans.trans, saanich)
GenePathComp$Month <- recode(GenePathComp$Month,
                        "4" = "Apr 2018",
                        "6" = "Jun 2018",
                        "8" = "Aug 2018",
                        "10" = "Oct 2018",
)

# Create nosZ scatterplot
x1_title <- expression(paste("Predicted"~italic("nosZ"), "gene abundance (clr-transformed)"))
nosPlot <- ggplot(GenePathComp, aes(K00376, delta_N2O)) + 
  geom_point(aes(shape=Month), size=3) +
  geom_smooth(method = "lm", se=FALSE, size=0.5, colour="black") +
  labs(x=x1_title, y = expression(Delta*N[2]*"O (nmol L"^-1*")")) +
  scale_y_continuous(limits=c(-14,20)) +
  scale_x_continuous(limits=c(-0.5,3.0)) +
  scale_shape_manual(values = c(22,1,2,9), 
                     breaks = c("Apr 2018", "Jun 2018", "Aug 2018", "Oct 2018")) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size=10),
        axis.text = element_text(size=8), legend.position = c(0.75, 0.75),
        legend.background = element_rect(colour = "black", size =0.25)) +
  annotate("text", x = 0, y = -12, label = expression("R"^2~"= 0.59"), size=3) +
  annotate("text", x = 0, y = -14, label = expression(italic("p =")~"9.6 x 10"^-6), size=3)
nosPlot

#ggsave(rates, filename = "Jameson_et_al_Fig3b.pdf",
#        width = 4.25, height=4, dpi = 300)