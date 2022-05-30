################################################################################
## Data import and compositional data transformations
################################################################################
library(zCompositions)
library(compositions)
library(dplyr)
library(phyloseq)

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
## Partial Least Squares Regression
################################################################################
library(mixOmics)

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

# Creat colour pallette to distinguish taxa
a<-as.data.frame(nspmls_asv@tax_table)
b<-as.data.frame(bactSelect_asv@tax_table)
t<-rbind(a,b)
tax.col = color.mixo(as.factor(t$Genus))
unique(tax.col)

# PLSR with 3 latent components
spls.select.asv.env <- pls(saanich.select.asv.comb, saanich.trim.common, 
                           mode = "regression", ncomp = 3)

# Determine appropriate number of model components
perf.pls = perf(spls.select.asv.env, validation="Mfold", folds = 10, 
                progressBar = FALSE, nrepeat = 100)
perf.pls$measures$Q2.total
perf.pls$measures$R2
plot(perf.pls, criterion = 'Q2.total')
plot(perf.pls, criterion = 'R2')

# Re-fit PLSR model with 1 component
spls.select.asv.env <- pls(saanich.select.asv.comb, saanich.trim.common, 
                           mode = "regression", ncomp = 1)

dev.off()
# Clustered coefficient heatmap (Fig. 1)
cim.select <- cim(spls.select.asv.env, row.sideColors = tax.col, 
                  row.names = TRUE, symkey = TRUE,
                  col.names = sample.traits, 
                  legend=list(legend = c("Nitrosopumilus", "Nitrosopelagicus", 
                                         "Nitrosopumilaceae", "SUP05", "SAR11", 
                                         "Nitrospina"), 
                  col = c("#F68B33", "#388ECC", "#C2C2C2", 
                          "#F0E442", "#CC79A7", "#009E73"),
                  cex=0.75, vjust=1), keysize = c(1,1))#, save = 'pdf'
)

plsr.correlations <- cim.select$mat
#write.csv(plsr.correlations, "plsr.correlation.coeffs.csv")

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
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
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
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

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
connectivity.traits.tax <- merge(connectivity.traits, tax[, c("ASV", "Phylum", "Class","Order",
        "Family", "Genus")], by="ASV")
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

plot1 <- ggplot(connectivity.traits.tax) +
  geom_point(aes(GS.Ammox, GS.Delta, colour = Subnetwork, size = kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"), 
                     name = "Subnetwork", palette = "Dark2") +
  labs(x = expression("ASV importance (Nitrification)"), y = NULL) +
  theme_bw()

plot1 <- plot1 + 
  scale_size_continuous(name="Connectivity") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(GS.Ammox, GS.Delta), shape=0, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(GS.Ammox, GS.Delta), shape=2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(GS.Ammox, GS.Delta), shape=4, size = 3) 
plot1

plot2 <- ggplot(connectivity.traits.tax) +
  geom_point(aes(GS.OxProd, GS.Delta, colour = Subnetwork, size = kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2") +
  labs(x = expression("ASV importance (N"*H[4]^"+"%->%N[2]*"O)"), y = NULL) +
  theme_bw()

plot2 <- plot2 + 
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_size_continuous(name="Connectivity") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(GS.OxProd, GS.Delta), shape = 0, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(GS.OxProd, GS.Delta), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(GS.OxProd, GS.Delta), shape = 4, size = 3) 
plot2

plot3 <- ggplot(connectivity.traits.tax) +
  geom_point(aes(GS.Reduc, GS.Delta, colour = Subnetwork, size = kWithin), alpha = 0.45) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3"),
                     name = "Subnetwork", palette = "Dark2") +
  labs(x = expression("ASV importance (N"*O[3]^"-"%->%N[2]*"O)"), y = NULL) +
  theme_bw()

plot3 <- plot3 + 
  scale_size_continuous(name="Connectivity") +
  theme(panel.grid = element_blank(), legend.position = c(0.85, 0.75), 
        legend.background = element_rect(linetype = 2, size = 0.5, colour = 1),
        legend.key.size = unit(0.1, 'mm')) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(GS.Reduc, GS.Delta), shape = 0, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(GS.Reduc, GS.Delta), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names3), ],
             aes(GS.Reduc, GS.Delta), shape = 4, size = 3) 
plot3

# Creat text for common y-axis
library(grid)
library(gridExtra)
y.grob <- textGrob(expression("ASV importance ("*Delta*N[2]*"O)"), 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot=90)

# Combine and save connectivity plots
connectivityPlots <- ggpubr::ggarrange(plot1, plot2,  plot3,
                        labels = "AUTO", # labels
                        align = "h", # Align them both, horizontal and vertical
                        nrow = 1, ncol = 3)  # number of rows

# Add labels
connectivityPlots <- grid.arrange(arrangeGrob(connectivityPlots, left = y.grob))
#ggsave(filename= "Connectivity_ASVsig.pdf", plot = connectivityPlots, 
#       width = 12.5, height = 4.5)


#####################################################################################################
## Propr Network analyses
#####################################################################################################
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
tax.source <-  tax %>% dplyr::rename(Partner = ASV)
tax.target <- tax %>% dplyr::rename(Pair = ASV)

NoRaresCorrALL <- merge(NoRaresCorrALL, tax.source[, c("Partner", "Genus")], by="Partner")
NoRaresCorrALL <- merge(NoRaresCorrALL, tax.target[, c("Pair", "Genus")], by="Pair")
NoRaresCorrALL <- merge(NoRaresCorrALL, test[, c("Partner", "kWithin")], by="Partner")
NoRaresCorrALL <- merge(NoRaresCorrALL, test[, c("Partner", "Subnetwork")], by="Partner")

NoRaresCorrALL.subset <- subset(NoRaresCorrALL, propr >= 0.60 | propr <= -0.60)

# write.csv(NoRaresCorrALL.subset, "arch.bact.propr.network.60.csv")