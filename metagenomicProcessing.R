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

# do not show data on stack bar if count (%) of a taxa is blow this CUTOFF
CUTOFF = 0.01

# scaling plot margins
sc=2

relOTU = TRUE

# order of bar plot
desired_order = c("C1", "C2", "C3", "F1", "F2", "F3", 
              "D1", "D2", "D3", "X0036.1" , "X0036.2", "X0036.3", 
              "X1887.1", "X1887.2", "X1887.3")


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

library(matrixStats)

theme_set(theme_bw())

library("wesanderson")
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
if (relOTU) {
  OTUm = OTUm / rep(colSums(OTUm), each = nrow(OTUm))
  # screen out taxa with abundance lower than CUTOFF
  idxSel  = rowMaxs(as.matrix(OTUm)) > CUTOFF
  OTUm = OTUm[idxSel, ]
  TAXm = TAXm[idxSel, ]
  }

physeqm = phyloseq(OTUm, TAXm) ## phyloseq object after grouping 


set.seed(123458)
colorset = sample(wes_palette(length(rownames(OTUm)), name = "Darjeeling1", type = "continuous"))


# generate abundance bar plots
p <- plot_bar(physeqm, fill = OTUlevel)  + theme(plot.margin = margin(6*sc,1*sc,6*sc,1*sc,"cm")) +
  theme(legend.position="right") + guides(fill=guide_legend(ncol=1)) +
  theme(text = element_text(size=20)) + scale_fill_manual(values = colorset)

# take care of plotting order
pd <- p$data
pd$Sample <- factor(pd$Sample, levels = desired_order)
p$data <- pd
print(p)

ggsave(filename = "myplot.png", plot = last_plot(), 
       width=40, height=40, unit="cm")


### ref. https://github.com/joey711/phyloseq/issues/616 
