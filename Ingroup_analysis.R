#############################################################################

# In group analysis Grevillea project   - Started 07/03/2024

#############################################################################
library(dartR)

gl.install.vanilla.dartR()
library(ape)
library(plotly)
library(StAMPP)
library(hierfstat)
library(directlabels)
library(pheatmap)
library(rrBLUP)
library(plotrix)
library(graph4lg)
library(gdata)
library(poppr)
library(RColorBrewer)
library(gplots)

sessionInfo()

gl <- gl.read.dart(filename="Report_DGr23-8087_SNP_2.csv",
                       ind.metafile="metadata_05032024_wo_clones.csv")

gl                                  # gives overview of the genlight object

# ind.metafile ids not matched were:
# [1] "OUT1-01" "OUT1-02" "OUT1-03" "OUT1-04" "OUT1-05" "OUT2-01" "OUT2-02"
# [8] "OUT2-03" "OUT2-04" "OUT3-03" "OUT3-04" "OUT3-05" "OUT4-01" "OUT4-02"
# [15] "OUT4-03" "OUT4-04"
# 
# DArT file ids not matched were: (CLONES)
#   [1] "AUS1-02" "AUS1-03" "BRO4-01" "BLA9-02" "YEL1-01" "BRO4-03" "TYP6-01"
# [8] "YEL1-02" "TYP6-02" "YEL1-03" "YEL1-04" "GRA4-05" "BLA9-05" "GRA2-06"
# [15] "YEL1-05" "BLA9-06" "TYP6-03" "YEL1-06" "TYP6-04" "YEL1-08"


# 146 genotypes,  175,577 SNPs , size: 164.6 Mb
# 
# missing data: 8533497 (=33.29 %) scored as NA

class(gl)<- "genlight"  
gl.report.rdepth(gl)

# No. of loci = 175577 
# No. of individuals = 146 
# Minimum      :  2.5 
# 1st quartile :  8 
# Median       :  16.2 
# Mean         :  21.41601 
# 3r quartile  :  29.8 
# Maximum      :  368.8 
# Missing Rate Overall:  0.33 

gl.smearplot(gl)

# FILTER 1: For POP GEN analysis, remove secondary (linked) snps to reduce dataset to single SNP per fragment

# This script filters out all but the first sequence tag with the same CloneID after ordering
# the genlight object on based on repeatability, avgPIC in that order (method='best') 
# or at random (method='random').

filter1 <- gl.filter.secondaries(gl, method="best", v=3)
filter2 <- gl.filter.rdepth(filter1, lower=5, upper=100, v=3)
filter3 <- gl.filter.reproducibility(filter2, t=0.99, v=3)
filter4 <- gl.filter.callrate(filter3, method="loc", threshold=0.90, v=3)
filter5 <- gl.filter.maf(filter4, threshold=0.02, verbose=3)
filter6 <- gl.filter.monomorphs(filter5, v=2)
filter7 <- gl.filter.callrate(filter6, method="ind", threshold=0.90)

filter7
# // 122 genotypes,  9,805 binary SNPs, size: 74.1 Mb
# 47006 (3.93 %) missing data

data <- gl.recalc.metrics(filter7)
# save genlight object (.Rdata) and/or csv file of filtered genotypes and locus metadata if wanted
gl.save(data, "Grevillea_9kloci.Rdata")

pcoa <- gl.pcoa(data, nfactors=5)

#2 axes looks good

gl.pcoa.plot(pcoa, data, pop.labels="pop", xaxis=1, yaxis=2, zaxis=NULL, interactive=TRUE)

###### SNMF analysis ##################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

library(LEA)

strgl <- data
index <- 1
pops <- sort(unique(strgl$pop))
for (pop in pops) {
  levels(strgl$pop)[match(pop, levels(strgl$pop))] <- index
  index <- index + 1
}

gl2structure(strgl, outfile = "9k_Grevillea.str", addcolumns = strgl$pop, exportMarkerNames = FALSE, outpath = ".")

#Edit in Textpad8

filename3 <- "9k_Grevillea_.str"
filename3
struct2geno(filename3, ploidy = 2, FORMAT=2, extra.col=2)

geno <- "9k_Grevillea.str.geno"

# need to set upper limit of K values and number of repetitions (100-1000 depending on size of dataset)
# all other parameters can be left as they are
project <- snmf(input.file = geno,
                K = 1:10,
                ploidy = 2,
                repetitions = 100,
                entropy = TRUE,
                percentage = 0.2,
                project = "new")

# indicate max K value and plot the cross-entropy criterion as a function of K to see where it plateaus (= optimal K)
kmax <- 10
cross <- sapply(1:kmax, function(k) { LEA::cross.entropy(project, run = 1, K = k) })
kplot <- plot(1:kmax, cross, main = "Cross Entropy", type="b", xlab = "K", ylab = "cross entropy")

print(k_plot)
png(filename = "Plots/kplot_9k.png")

#looks like it is plateauing between 4 - 8 

# for chosen K value and run, choose # colours needed and plot barchart 
# can also export the Q values if prefer to plot in excel

kchoice=6
mycol=palette.colors(kchoice)
par(mfrow = c(7,1),  mai = c(0.2, 0.3, 0.1, 0.1))

kchoice=2
p2 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol , border=TRUE, xaxt="n", ylab = "")

kchoice=3
p3 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xaxt="n", ylab = "")

kchoice=4
p4 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xaxt="n", ylab = "")

kchoice=5
p5 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xaxt="n", ylab = "")

kchoice=6
p6 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xaxt="n", ylab = "")

kchoice=7
p7 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xlab = "Individuals", ylab = "")

kchoice=8
p8 <- barchart(project, K=kchoice, run=65, sort.by.Q = FALSE, lab = FALSE, col=mycol, border=TRUE, xlab = "Individuals", ylab = "Admixture coefficients")


# save data for chosen k value/s
qvalues<-LEA::Q(project, K=kchoice, run=65)
write.table(qvalues, "Grevillea_9k_k_values.txt", quote=FALSE)


# Graphic now exported to Inkscape and populations added manually by counting the number of
# individuals for each pop and drawing a red line through all seven column graphs.

#####################################################################################

#                        DIVERSITY + F-sTATS for Grevillea 11/03/2024

# Note: Pop-level diversity measures analysed using automated script known as Frankenstein 

######################################################################################

# OPTIONAL: if needing to change the pop labels to another ind metrics slot, e.g. lineage, species, sex etc. Just don't forget to change back!
pop(data) <- data@other$ind.metrics$pop
pop(data) <- data@other$ind.metrics$species

# estimates standard He, Ho and Fis and standard error per population (can also estimate per individual, method="ind")
# ignores estimates of adjusted Ho/He, see next section if want to estimate these parameters
het_stats_pop <- gl.report.heterozygosity(data, method="pop")
het_stats_loci <- utils.basic.stats(data)

# There are populations with one individual. Please remove populations with one individual or merged them with other populations for his function to work
# Error in colSums(apply(sgl_mat[[y]], 2, is.na)) : 
#   'x' must be an array of at least two dimensions

table(pop(data))  

# AUS1 AUS2 AUS3 BLA1 BLA2 BLA3 BLA4 BLA5 BLA6 BLA7 BLA8 BLA9 BRO1 BRO2 BRO3 BRO4 GRA1 GRA2 GRA3 GRA4 TYP1 
# 6    4    6    6    5    4    4    5    5    4    6    3    6    7    5    3    6    3    6    4    5 
# TYP2 TYP3 TYP4 TYP5 TYP6 
# 5    3    5    5    1 

#So, we have to remove Type6, which only has one individual for the population.

# will delete all individuals in all populations except those listed.
data_new <- gl.drop.pop(data, pop.list=c("BLA9", "BRO4","GRA2","TYP3","TYP6"))

# estimate pairwise fst values
fst <- gl.fst.pop(data_new, nboots=1000, percent=0.95)
pairwise_fst <- as.matrix(fst$Fsts)
pairwise_fst
final_fst<- pairwise_fst

# custom re-ordering of pairwise fst matrix for table or heatmap
sym_fst <- as.matrix(as.dist(fst$Fsts))
View(sym_fst)

table(pop(data_new))


order <- c("AUS1","AUS2","AUS3",                               #Australis - South
           "BLA1","BLA3","BLA4", "BLA5","BLA7","BLA8",         #Australis - North
           "GRA1","GRA3","BLA2","TYP1",                        #Grandis
           "GRA4","BLA6", "TYP2","TYP4","TYP5",                #Brachystylis
           "BRO2", "BRO1","BRO3")                              #Bronweneae

final_fst <- reorder_mat(mat=sym_fst, order=order)
final_fst

# export fst values to csv file
write.csv(final_fst, "Grevillea_pop_fst_Final.csv")

# basic heatmap of fst values - not included in manusciprt
npop <- nPop(data_new)
mycol <- brewer.pal(n=9, name="Reds")
heatmap_def <- heatmap.2(final_fst, notecol="black", density.info="none", trace="none", Rowv=NA, Colv=NA, dendrogram=c("none"), col=mycol, colsep=c(0:npop), rowsep=c(0:npop), sepcolor="black", sepwidth=c(0.02, 0.02))

# calculate global fst/fis values (overall values at bottom of table)
gl.basic.stats(data_new)

# overall
# Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.0437 0.0640 0.1414 0.0773 0.1453 0.0812 0.5470 0.5591 0.3174 0.0868 

#-----------------------------------------------------------------------------------------
##12/03/2024

# next, we need to work out what samples belong in what lineage.
# To do this, I am going to use another program called Splitstree. The following code is
# to export the SNP data into a compatible format that plays nicely with this software.

# export distance matrix for SPLITSTREE input
# generates fasta file, converts to DNAbin object, calculates distance metric and exports in the correct format for splitstree
library(pofadinr)

install.packages("devtools")

#And then, install the pofadinr package from the github repository.

library(devtools)

install_github("simjoly/pofadinr")

#^ Doesn't work

isfar <- get(load('c:/Desktop/Stump/004_Data/001_SNP_Geno_G_brachystylis/Grevillea_9kloci.Rdata'))
class(a)<- "genlight"  
gl.report.rdepth(a)

data_fasta <- gl2fasta(data2, method=3, outfile="Grev_9k.fas", outpath=getwd())
filename = "Grev_9k.fas"

install_github("simjoly/pofadinr")

data_dnabin <- fasta2DNAbin(filename)

library(pofadinr)

genpofad_dist <- dist.snp(data_dnabin, model = "GENPOFAD", pairwise.deletion = TRUE, as.matrix=TRUE)  # can also use model="MATCHSTATES"
outfile <- file("Grev_9k.txt", open = "w")
writeLines(as.character(nInd(data2)), con = outfile)
write.table(genpofad_dist, file = outfile, sep = " ", append=TRUE, col.names=FALSE, quote=FALSE)

close(outfile)




