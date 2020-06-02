# specify working directory
setwd("C:/Users/pakpoomsu/Desktop/Banana_Metagenomics/")
# specify data file (csv)
rawData = read.csv('HN00126595_Report_otu_table_shared_with_tax_assignment.xlsx - Sheet1.csv', sep = ',', header = T)
# specify which taxonomic level to create stack plot (Phylum, Class, Order, Family, Genus, Species)
OTUlevel = "Phylum"

# select only relevant data 
otumat = as.matrix(rawData[,15:dim(rawData)[2]])
otumat = otumat[,1:15]  # use only banana data, exclusing zerg data
rownames(otumat) = rawData$Group
taxmat = as.matrix(rawData[,3:9])
rownames(taxmat) = rawData$Group

# scaling plot margins
sc=3

relOTU = TRUE


###########################################
####### load libraries ###################

library(phyloseq)
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("scales")
packageVersion("scales")

library("grid")
packageVersion("grid")

theme_set(theme_bw())



######################################################
########### Processing data before plotting #########

## convert data into phyloseq compatible forms
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX) # initial phyloseq object from raw data 

phyloGlom = tax_glom(physeq, OTUlevel) # group data by pre-specified OTUlevel

## for taxonomic matrix, exclude data at levels below OTUlevel
TAXmt = tax_table(phyloGlom)
indTaxLevel = match(OTUlevel,colnames(TAXmt))
TAXm = TAXmt[, 1:indTaxLevel]
## get otu table (after grouping)
OTUm = otu_table(phyloGlom, taxa_are_rows = TRUE)

# convert absolute to relative count as specified
if (relOTU) {OTUm = OTUm / rep(colSums(OTUm), each = nrow(OTUm))}

physeqm = phyloseq(OTUm, TAXm) ## phyloseq object after grouping 

# generate abundance bar plots
plot_bar(physeqm, fill = OTUlevel) + theme(plot.margin = margin(6*sc,1*sc,6*sc,1*sc,"cm")) +
  theme(legend.position="right") + guides(fill=guide_legend(ncol=5))
par(pin=c(1.9,1.9))  
ggsave(filename = "myplot.png", plot = last_plot(), 
       width=80, height=80, unit="cm")



### ref. https://github.com/joey711/phyloseq/issues/616 
