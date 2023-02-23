# Load packages and data
library(reshape2)
library(multtest)
library(lattice)
library(devtools)
library(ggbiplot)
library(plyr)
library(scales)
library(ggplot2)
library(grid)
library(tidyverse)
library(cowplot)
library(ggthemes)
library(broom)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(gplots)

# Link to the dada2 `strain_table.tsv` output file
infile <- "../01_dada2_out/strain_table.tsv"

# Let's make some output directories
dir.create("Indexed")
dir.create("Read_counts")
dir.create("Plots/Strains", recursive = TRUE)

# Make a vector list of all input file names
myData <-  read.table(file = infile, sep = '\t', header = TRUE)
sample.names <- sapply(strsplit(colnames(myData), "[.]"), `[`, 1)
sample.names <- sample.names[-c(1)]

# And read in the get Indexing file
# WARNING, certain Indexing names in Alleles will not comply with differential equation function (notably x and df). 
Alleles <- read.delim("../Allele_indexing.txt", sep = '\t', header=TRUE)
Indexing <- read.delim("../Indexing.txt", header=TRUE)
#sapply(Alleles, class)

# First we setup the differential equation
solve_allele_equation <- function(data, n_obs, allele_col, equation) {

  # create a data.frame to work on
  df <- data

  # assign the number of observations to the correct allele
  for (j in 1:nrow(df)) {
    assign(x = df[[allele_col]][j], value = df[[n_obs]][j] )
  }

  # solve the equations
  x <- sapply(df[[equation]], function(x) eval(parse(text = x)) )
  names(x) <- NULL

  # add the equation solutions to the original data.frame
  df[["Total_C"]] <- x

  return(df)
}


# Let's loop the analysis
# First make a file to append to.
All_reads <- data.frame(matrix(NA, ncol=26, nrow=0))[-1]

for (i in sample.names) {

	# Input the data, subselect, and adapt naming
	myData <- read.table(file = infile, sep = '\t', header = TRUE)
	myData <- myData[c("strain",i)]
	colnames(myData) <- c("Barcode","Total")

	# Change Indexing of strains and alleles
	myData2 <- gsub('_C12W1_', ',', myData$Barcode)
	myData <- cbind(myData,myData2)
	myData <- myData %>% separate(myData2, c("ID.1", "allel"), sep = ",",remove = TRUE, convert = FALSE)
	myData3 <- join(Alleles, myData, by = "Barcode")
	#sapply(myData3, class)

	# Solve differential equation
	myData4 <- solve_allele_equation(data = myData3, 
	n_obs = "Total", 
	allele_col = "Diffrential_ID", 
	equation = "Differential_equation")

	# The equation produces NA for 0/0, change to 0 again
	myData4$Total_C[is.na(myData4$Total_C)] <- 0

	# Generate sample ID for indexing experimental samples
	Name <- i
	ID <- gsub('.tsv', '', Name)
	N <- nrow(myData4)
	Sample <- rep(ID, each = N)
	Sample <- as.data.frame(Sample)
	colnames(Sample) <- c("ID")
	Index <- join(Sample, Indexing, by = "ID")

	# Fuse the sample indexing with data
	myData5 <- cbind(Index, myData4)

	# And normalize amplicon counts on trimmed data
	Nreads <- sum(myData5$Total_C)
	myData5$Relative_abundance <- (myData5$Total_C/Nreads)
	Read_counts <- cbind(ID, Nreads)
	#sapply(myData5, class)

	# Remove * on allele so it can be used for indexing, and second ID column
	myData5$allel <- gsub('*', '', myData5$allel)

	# We'll save the header and export data
	output <- paste("Indexed/",ID,".csv", sep="")
	output2 <- paste("Read_counts/",ID,".csv", sep="")
	output3 <- paste("Plots/",ID,".pdf", sep="")

	header <- colnames(myData5)
	write.table(header, "Indexed/1_header.tmp", sep=",",  col.names=FALSE, row.names =FALSE)
	write.table(myData5, output, sep=",",  col.names=FALSE, row.names =FALSE)
	write.table(Read_counts, output2, sep=",",  col.names=FALSE, row.names =FALSE)

	# Append data to dataframe All_reads
	All_reads <- rbind(All_reads, myData5)

	# Let's also make some plots
	Titel <- unique(myData5$Sample)

	Fig1 <- ggplot(myData5, stats = "identity", aes(Strain, fill = allel)) +
	  geom_col(data=myData5, aes(Strain, Relative_abundance)) +
	  geom_hline(yintercept = 0.017) +
	  background_grid(major = "none", minor = "none") + # add thin horizontal lines
	  panel_border(colour = "black", size = 1) + # and a border around each panel
	  labs (title = Titel) +
	  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
	  theme(panel.spacing = unit(0.1, "lines")) +
	  theme(legend.title=element_blank()) +
	  theme(legend.text=element_text(size=10)) +
	  theme(text=(element_text(size=10))) +
	  theme(axis.text=(element_text(size=5))) +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	  theme(panel.background = element_blank()) +
	  theme(legend.key = element_blank()) +
	  theme(legend.text = element_text(face = "italic")) +
	  theme(aspect.ratio=0.25) +
	  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
	  theme(legend.position ="top",
	        legend.margin = margin(2, 2, 2, 2),
	        legend.box.background = element_rect(fill='white'),
	        legend.background = element_blank(),
	        legend.key = element_rect(fill = NA, color = NA),
	        legend.spacing.x=unit(0, "cm"),
	        legend.spacing.y=unit(0, "cm"))

	pdf(output3)
	print(Fig1)
	dev.off()
}

header <- read.csv("Indexed/1_header.tmp", header = F)
header <- t(header)
colnames(All_reads) <- (header)
write.table(All_reads, "Indexed/All_reads.csv", sep=",",  col.names=TRUE, row.names =FALSE)

# Count the occurrences of the wrong strain in the wrong population (False positives)
VGs <- subset.data.frame(All_reads, grepl('VG', All_reads$Strain))
VGs <- subset.data.frame(VGs, grepl('GP', VGs$Experiment))
VG_FP <- length(VGs$Population[VGs$Total_C > 0])

GPs <- subset.data.frame(All_reads, grepl('GP', All_reads$Strain))
GPs <- subset.data.frame(GPs, grepl('VG', GPs$Experiment))
GP_FP <- length(GPs$Population[GPs$Total_C > 0])

# The number of indisputable false positives is
#VG_FP+GP_FP

# Graphical interpretation of data####
# Let's make a graph of only the MMs abundances.
MMsamples <- subset.data.frame(All_reads, grepl("0", All_reads$Timepoint))
MMsamples <- MMsamples[order(MMsamples$Strain),]
MMsamples$Relative_abundance[MMsamples$Relative_abundance == 0] <- NA
#sapply(MMsamples, class)
AveragesMM <- ddply(MMsamples, c("Experiment", "Strain", "allel"), summarise,
                     mean = mean(Relative_abundance), sd = sd(Relative_abundance))

# This way hard to make errorbars right.
Fig4 <- ggplot(AveragesMM, stats = "identity", aes(Strain, fill = allel)) +
  geom_bar(data=AveragesMM, aes(Strain, mean), position = "stack", stat = "identity") +
  facet_grid(rows = vars(Experiment)) +
  geom_hline(yintercept = 1/29) +
  geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=0.2), stat = "identity") + #position = position_dodge(0.3
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  scale_fill_manual(values = c("forestgreen", "darkolivegreen2", "chartreuse", "blue3", "dodgerblue2", "deepskyblue1", "darkorchid1", "grey")) +
  labs (x="Strain", y="Fraction of amplicons") +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position= c(0.83, 0.82),
        legend.margin = margin(2, 2, 2, 2),
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

pdf("Plots/Fig4.pdf")
print(Fig4)
dev.off()

# This way hard to remove borders around bars, and keep errorbbars.
Fig5 <- ggbarplot(MMsamples, x = "Strain", y = "Relative_abundance", add = "mean_sd", fill = "allel", size = 0,  facet.by = "Experiment") + #color=NA removes border but also error bars, ar
  facet_grid(rows = vars(Experiment)) +
  geom_hline(yintercept = 1/29) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  scale_fill_manual(values = c("forestgreen", "darkolivegreen2", "chartreuse", "blue3", "dodgerblue2", "deepskyblue1", "darkorchid1", "grey")) +
  labs (x="Strain", y="Fraction of amplicons") +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position= c(0.83, 0.82),
        legend.margin = margin(2, 2, 2, 2),
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

pdf("Plots/Fig5.pdf")
print(Fig5)
dev.off()

# This code can be used to pull out any single strains alleles by name
#subset.data.frame(MMsamples, grepl("GP2-4_52", MMsamples$Strain))

# Let's pick one strain at a time and plot them individualy (Then Loop this)
df_uniq <- unique(All_reads$Strain)
colourCount = length(unique(mtcars$hp))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Let's make an Dataframe to save all ratios in
AlleleRatios <- data.frame(ID=character(),
                           Barcode1=character(),
                           Barcode1_RA=numeric(),
                           Total1=numeric(),
                           Barcode2=character(),
                           Barcode2_RA=numeric(),
                           Total2=numeric(),
                           AllelRatio=numeric(),
                 stringsAsFactors=FALSE)

#quartz()

# THIS LOOP KEEPS FAILING DUE TO GGPLOT ERRORWarningsS (in last FigStrains part)
# Stop loop before plots to finish rest of script

for (i in df_uniq) {

	OneStrain <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
	PIPT <- subset.data.frame(OneStrain, grepl(pattern = 'GP|VG', OneStrain$Experiment))
	PIPT <- subset.data.frame(PIPT, grepl(pattern = 'GP|VG', PIPT$Population))
	Strain2 <- subset.data.frame(OneStrain, grepl('Strain', OneStrain$Experiment))
	output <- paste("Plots/Strains/",i,".pdf", sep="")
	output2 <- paste("Plots/Strains/",i,"_Strains.pdf", sep="")

	# Let's plot how the allele ratios vary with abundance. Should be no correlation and close to 1
	OneStrain2 <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
	Barcode1 <- subset.data.frame(OneStrain2, grepl("C12W1_1", OneStrain2$Barcode))
	Barcode1 <- select(Barcode1, c('ID',  "Barcode", "Relative_abundance", "Total_C"))
	Barcode2 <- subset.data.frame(OneStrain2, grepl("C12W1_2", OneStrain2$Barcode))
	Barcode2 <- select(Barcode2, c('ID',  "Barcode", "Relative_abundance", "Total_C"))
	OneStrain3 <- join(Barcode1, Barcode2, by = "ID")
	colnames(OneStrain3) <- (c('ID', "Barcode1", "Barcode1_RA", "Total1", "Barcode2", "Barcode2_RA", "Total2"))
	OneStrain3$AllelRatio <- OneStrain3$Barcode1_RA/OneStrain3$Barcode2_RA
	OneStrain3$ReadSum <- OneStrain3$Total1+OneStrain3$Total2

	# Append results here
	AlleleRatios <- rbind(AlleleRatios, OneStrain3)
	#End below here to avoid crash on figures
#}

	# This needs to be here for some annoying reason, loop terminates otherwise
	#quartz()

	PIPT$Replicate <- as.character(as.integer(PIPT$Replicate))
	suppressWarnings(FigRA <- ggplot(PIPT, aes(Timepoint, Relative_abundance), message=FALSE) +
                   geom_point(mapping = aes(shape = Barcode, color = Replicate), size = 2) +
                   facet_wrap(~Experiment + Treatment) +
                   theme(strip.text.x = element_text(size=12, angle=0, face = "italic"),
                         strip.background = element_rect(colour="white", fill="white")) +
                   geom_path(mapping = aes(color = Replicate, linetype = Barcode), size = 0.5) +
                   scale_color_manual(values = getPalette(colourCount)) +
                   coord_cartesian(ylim=c(0.0001, 1), xlim=c(-1, 44), expand = F) + #ylim=c(-0, 12.5)
                   scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                 labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                   background_grid(major = "none", minor = "none") + # add thin horizontal lines
                   panel_border(colour = "black", size = 1) + # and a border around each panel
                   labs (x=("Time (days)"), y=("Relative abundance"), title = i) +
                   theme(panel.spacing = unit(0.1, "lines")) +
                   theme(legend.title=element_text(size=12)) +
                   theme(legend.text=element_text(size=6)) +
                   theme(text=(element_text(size=12))) +
                   theme(axis.text=(element_text(size=12))) +
                   theme(panel.background = element_blank()) +
                   theme(legend.text = element_text(face = "italic")) +
                   theme(aspect.ratio=1))

	suppressWarnings(pdf(output))
	suppressWarnings(print(FigRA))
	suppressWarnings(dev.off())

	# End here to also make figures (sometimes crashes the loop unpredictably for me)
#}
	#quartz()
	# This figure is only relevant for plotting the occurences of alleles in the singe genotype samples (P21502_101-164)
	suppressWarnings(FigStrains <- ggplot(Strain2, aes(Sample, Relative_abundance)) + #, message=FALSE
                   geom_point(mapping = aes(color = Barcode), size = 2, shape = 1) +
                   coord_cartesian(ylim=c(0.0001, 1), expand = F) + #ylim=c(-0, 12.5)
                   scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                 labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                   background_grid(major = "none", minor = "none") + # add thin horizontal lines
                   panel_border(colour = "black", size = 1) + # and a border around each panel
                   labs (x=("Strain Genotype sample"), y=("Relative abundance"), title = i) +
                   theme(panel.spacing = unit(0.1, "lines")) +
                   theme(legend.title=element_text(size=12)) +
                   theme(legend.text=element_text(size=6)) +
                   theme(text=(element_text(size=10))) +
                   theme(axis.text=(element_text(size=6))) +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                   theme(panel.background = element_blank()) +
                   theme(legend.text = element_text(face = "italic")) +
                   theme(aspect.ratio=0.25) +
                   theme(legend.position ="top"))

	suppressWarnings(pdf(output2))
	suppressWarnings(print(FigStrains))
	suppressWarnings(dev.off())

	# End here to also make genotype figures (crashes the loop unpredictably for me)
}

# Allele ratios and correlations####

# Change 0 and inf to NA to fix plotting of ratios
AlleleRatios2 <- join(AlleleRatios, Indexing, by = "ID")
AlleleRatios2$Total1 <- as.numeric(as.integer(AlleleRatios2$Total1))
AlleleRatios2$Total2 <- as.numeric(as.integer(AlleleRatios2$Total2))
AlleleRatios2$Timepoint <- as.numeric(as.integer(AlleleRatios2$Timepoint))
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == 0] <- NA
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == "inf"] <- NA
AlleleRatios2$AllelRatio[sapply(AlleleRatios2$AllelRatio, is.infinite)] <- NA
AlleleRatios2 <- na.omit(AlleleRatios2)


# Let's look at the distribution of some values and parameters
#min(AlleleRatios2$AllelRatio)
#max(AlleleRatios2$AllelRatio)

pdf("Plots/Strains/HistRatios.pdf")
hist(log10(AlleleRatios2$AllelRatio), breaks=80, main="AllelRatio") #xlim=c()
dev.off()

#min(AlleleRatios2$Total1)
#max(AlleleRatios2$Total1)

#min(AlleleRatios2$Total2)
#max(AlleleRatios2$Total2)

#min(AlleleRatios2$ReadSum)
#max(AlleleRatios2$ReadSum)
#median(AlleleRatios2$ReadSum)

pdf("Plots/Strains/Barcode1_RA_hist.pdf")
hist(log10(AlleleRatios2$ReadSum), breaks=80, main="BarcodeBoth_abundance") #xlim=c()
dev.off()

#unique(AlleleRatios2$Population)

# May need to remove RO5AC samples for this
AlleleRatios2 <- subset.data.frame(AlleleRatios2, grepl(pattern = 'GP|VG', AlleleRatios2$Population))

colourCount = 6
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

FigRatio <- ggplot(AlleleRatios2, aes(Barcode1_RA, AllelRatio), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE)), shape=1, size = 1) +
  facet_wrap(~Barcode1) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta")) +
  coord_cartesian(ylim=c(0.1, 10), xlim=c(0.0001, 1), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 2^x),
                labels = trans_format("log10", math_format(2^.x))) + #change the scale on y axis
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=("Relative abundance allel 1)"), y=("Ratio allel 1 to 2")) + #title = i
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.text=element_text(size=5)) +
  theme(text=(element_text(size=5))) +
  theme(axis.text=(element_text(size=5))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="top")

pdf("Plots/Strains/FigRatios.pdf")
print(FigRatio)
dev.off()

# This we can then do to remove 0+0 observations of alleles for correlation (Outside loop)
AlleleRatios3 <- join(AlleleRatios, Indexing, by = "ID")
AlleleRatios3$ReadSum[AlleleRatios3$ReadSum == 0] <- NA
AlleleRatios3 <- AlleleRatios3 %>% drop_na(ReadSum)
AlleleRatios3$Total1 <- as.numeric(as.integer(AlleleRatios3$Total1))
AlleleRatios3$Total2 <- as.numeric(as.integer(AlleleRatios3$Total2))
#sapply(AlleleRatios3, class)

# Remove Strain observations
AlleleRatios3 <- subset.data.frame(AlleleRatios3, grepl(pattern = 'GP|VG', AlleleRatios3$Experiment))

# Look at distribution
#max(AlleleRatios3$Total1)
#min(AlleleRatios3$Total2)

pdf("Plots/Strains/Allel1_counts.pdf")
hist(log10(AlleleRatios3$Total1), breaks=80, main="Allel1_counts") #xlim=c()
dev.off()

pdf("Plots/Strains/Allel2_counts.pdf")
hist(log10(AlleleRatios3$Total2), breaks=80, main="Allel2_counts") #xlim=c()
dev.off()

# Need to remove RO5AC samples
AlleleRatios3 <- subset.data.frame(AlleleRatios3, grepl(pattern = 'GP|VG', AlleleRatios3$Population))

# Look for a barcode
#OneStrain <- subset.data.frame(AlleleRatios3, grepl("P21502_101_C12W1_1*", AlleleRatios3$Barcode1))

StrainIndex <- subset(Alleles, select = c(1,6))
colnames(StrainIndex) <- c("Barcode1", "Strain")
AlleleRatios3 <- inner_join(AlleleRatios3, StrainIndex, by = "Barcode1")
#length(unique(AlleleRatios3$Strain))

FigReg <- ggplot(AlleleRatios3, aes(Total1+1, Total2+1), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE), shape = Treatment), size = 1.2) +
  facet_wrap(~Strain) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = c("black", "orange", "red", "blue", "cyan", "magenta")) +
  scale_shape_manual(values=c(1, 2)) +
  coord_cartesian(ylim=c(0.5, 28000), xlim=c(0.5, 28000), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=("Read count allele 1"), y=("Read count allele 2")) + #title = i
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.text=element_text(size=6)) +
  theme(text=(element_text(size=6))) +
  theme(axis.text=(element_text(size=6))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="top")


pdf("Plots/Strains/FigReg.pdf")
print(FigReg)
dev.off()

#Clustering of samples####

# We also need a dataframe for cbinding Matrix data
MatrixAllR <- data.frame(unique(All_reads$Sample))
colnames(MatrixAllR) <- (c('Sample'))
#MatrixAllR
#df_uniq

for (i in df_uniq) {
  #Lets also make a matrix format for further analysis
  Reads <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain), select = c("Sample","Barcode", "Relative_abundance"))
  Barcodes <- unique(Reads$Barcode)
  Barcodes <- gsub("\\*", "", Barcodes)
  Barcodes <- gsub("\\?", "", Barcodes)
  
  Read1 <- subset.data.frame(Reads, grepl(paste(Barcodes[1]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read1) <- (c("Sample", paste("",i,"_",paste(Barcodes[1],""), sep = "")))
  
  Read2 <- subset.data.frame(Reads, grepl(paste(Barcodes[2]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read2) <- (c("Sample", paste("",i,"_",paste(Barcodes[2],""), sep = "")))
  Reads_matrix <- join(Read1, Read2, by = "Sample")
  
  Read3 <- subset.data.frame(Reads, grepl(paste(Barcodes[3]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read3) <- (c("Sample", paste("",i,"_",paste(Barcodes[3],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read3, by = "Sample")
  
  Read4 <- subset.data.frame(Reads, grepl(paste(Barcodes[4]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read4) <- (c("Sample", paste("",i,"_",paste(Barcodes[4],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read4, by = "Sample")
  
  # Remove NA columns since not all strains have 4 alleles
  Reads_matrix <- Reads_matrix[ , apply(Reads_matrix, 2, function(x) !any(is.na(x)))]
  
  # Then we just join this with growing matrix
  MatrixAllR <- join(MatrixAllR, Reads_matrix, by = "Sample")
  
}

#sapply(MatrixAllR, class)

# Let's change the header now
MatrixHeadR <- colnames(MatrixAllR)
MatrixHeadR <- gsub("_P21502_[0-9]{3}_C12W1", "", MatrixHeadR)
MatrixAll2R <- MatrixAllR
colnames(MatrixAll2R) <- MatrixHeadR


# Reformat to Matrix
V <- as.vector(MatrixAll2R$Sample)
MatrixAll2R <- as.matrix(MatrixAll2R)
class(MatrixAll2R) <-"numeric"
MatrixAll2R <- MatrixAll2R[,-1]
row.names(MatrixAll2R)=V
MatrixAll2R <- t(MatrixAll2R)

# I would like to blank out 0 values but I can't get the heatmap code to work then
#MatrixAll2R[MatrixAll2R == 0] <- NA

MatrixAll2R_cF <- scale(MatrixAll2R, scale = T, center = F) #center = F, scale = T)
MatrixAll2R_cT <- scale(MatrixAll2R, scale = T, center = T) #center = T, scale = T)
MatrixAll2R_log <- log10(MatrixAll2R+1)
MatrixAll2R_cFt <- t(MatrixAll2R_cF)

# Make a heatmap without clustering
colRamp <- colorRampPalette(c("red", "black", "green"), space="rgb")(64)
pdf("Plots/Strains/Heat1.pdf")
heatmap(MatrixAll2R, col = colRamp, na.rm = TRUE, distfun = dist, cexRow=0.3, cexCol=0.5) #Rowv = NA, Colv = NA
dev.off()

#2#Apply clustering (Eucladian distances)
#heatmap(MatrixAll2R, col = colRamp, distfun = dist, revC = T, #margins = c(8,6)
#        xlab = "Sample", ylab = "Allel", scale = "none")

# This prioritizes strain first
pdf("Plots/Strains/Heat_cT_Eucl.pdf")
heatmap(MatrixAll2R_cT, col = colRamp, cexRow=0.3, cexCol=0.5, distfun = dist, revC = T, margins = c(8,6),
        xlab = "Sample", ylab = "Allel", scale = "none")
dev.off()

heatmap.2(MatrixAll2R, col = colRamp, scale="none", revC = F, margins = c(8,6),
          xlab = "Sample", ylab = "Allel", trace="none", cexRow=0.4)

breaks <- seq(min(MatrixAll2R, na.rm = T), max(MatrixAll2R, na.rm = T), length.out = 65)

pdf("Plots/Strains/Heat_cT_Green_row.pdf")
heatmap.2(MatrixAll2R_cF, col = colRamp, scale="none", cexRow=0.3, cexCol=0.5, revC = F, margins = c(8,6),
          xlab = "Sample", ylab = "Allel", trace="none", na.color = "white", breaks=breaks)
dev.off()

# Let's make a PCA and then add indexing
PCA1 <- prcomp(na.omit(MatrixAll2R_cT), center = TRUE, scale = TRUE)
Indexing2 <- as.data.frame(Barcode1$ID)
colnames(Indexing2) <- c("ID")
#sapply(Indexing2, class)
Indexing3 <- join(Indexing2, Indexing, by = "ID")

#PCA1
# Plot the eigenvalues of the PC components
#plot(PCA1)
#summary(PCA1)

# Export PCA results
PCAdata <- PCA1$rotation
#class(PCAdata)
#dim(PCAdata)
#plot(PCAdata[,1])

# Add Population /species/Strain) vectors
PCAdata_df <- as.data.frame(cbind(PCAdata, Indexing3))
#dim(PCAdata_df)
#head(PCAdata_df)
PCAdata_df$PC1 <- as.numeric(as.character(PCAdata_df$PC1))
PCAdata_df$PC2 <- as.numeric(as.character(PCAdata_df$PC2))
PCAdata_df$PC3 <- as.numeric(as.character(PCAdata_df$PC3))
#unique(PCAdata_df$Population)
#unique(PCAdata_df$Population)

# Make plots
plot(PCAdata_df$PC1, PCAdata_df$PC2, col = PCAdata_df$Timepoint, xlab = "PC1 (25%)", ylab = "PC2 (15%)",
     cex.lab = 1.3, cex.axis=1.3, pch = c(PCAdata_df$Population), cex =  1.5)
# "'legend' is of length 0"?
#legend("bottomleft", cex = 1.2, legend=levels(PCAdata_df$Bottel), pch = c(1,2,3,4,5), col = c(1,2,3,4,5))
text(PCAdata, labels = ID , pos = 1, cex = 0.7)


PCA2 <- prcomp(na.omit(t(MatrixAll2R_cT)), scale = TRUE)
#dim(PCA2)
#head(PCA2)

#ggbiplot(PCA1, ellipse=TRUE)	# Can't test this
#ggbiplot(PCA2, ellipse=TRUE)	# Can't test this

#StrainCounts#####
# Now I will modify the analysis to make quantitative predictions about strain density using barcodes

# Let's loop the analysis
# First make a file to append to
All_Strains <- data.frame(matrix(NA, ncol=26, nrow=0))[-1]
header <- read.csv("Indexed/1_header.tmp", header = F)
header <- t(header)
colnames(All_reads) <- (header)

# The loop below hasn't been updated yet
for (i in InputFiles) {
  # Import data (too be looped)
  input <- paste("Input/",i,"", sep="")
  myData <- read.table(file = input, sep = '\t', header = TRUE)
  
  # Change Indexing of strains and alleles
  myData2 <- gsub('_C12W1_', ',', myData$Barcode)
  myData <- cbind(myData,myData2)
  myData <- myData %>% separate(myData2, c("ID.1", "allel"), sep = ",",remove = TRUE, convert = FALSE)
  myData3 <- join(Alleles, myData, by = "Barcode")
  #sapply(myData3, class)
  
  # Solve differential equation
  myData4 <- solve_allele_equation(data = myData3, 
                                   n_obs = "Total", 
                                   allele_col = "Diffrential_ID", 
                                   equation = "Differential_equation_strains")
  
  # The equation produces NA for 0/0, change to 0 again
  myData4$Total_C[is.na(myData4$Total_C)] <- 0
  
  # Now we can remove the second alleles
  myData4 <- subset.data.frame(myData4, grepl('YES|NO', myData4$Modified_allelcounts))

  # Generate sample ID for indexing experimental samples
  Name <- i
  ID <- gsub('.tsv', '', Name)
  N <- nrow(myData4)
  Sample <- rep(ID, each = N)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- c("ID")
  Index <- join(Sample, Indexing, by = "ID")
  
  # Fuse the sample indexing with data
  myData5 <- cbind(Index, myData4)
  
  # And normalize amplicon counts on trimmed data
  Nreads <- sum(myData5$Total_C)
  myData5$Relative_abundance <- (myData5$Total_C/Nreads)
  Read_counts <- cbind(ID, Nreads)
  #sapply(myData5, class)
  
  # Remove * on allel so it can be used for indexing, and second ID column
  myData5$allel <- gsub('*', '', myData5$allel)
  
  # We'll save the header and export data
  header <- colnames(myData5)
  
  # Append data to dataframe All_reads
  All_Strains <- rbind(All_Strains, myData5)

}

write.table(All_Strains, "Indexed/All_strains.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#Check the results####

#Subselect C samples for VG
VG_C <- subset.data.frame(All_Strains, grepl('Control', All_Strains$Treatment))
VG_C <- subset.data.frame(VG_C, grepl('VG', VG_C$Population))

FigPie_StrainCounts <- ggplot(VG_C, aes(x = "", y = Relative_abundance, fill = Strain)) +
  geom_col() +
  coord_polar(theta = "y") +
  facet_grid(rows = vars(Timepoint), cols = vars(Replicate)) +
  theme_void() 

pdf("Plots/FigPie_StrainCounts")
print(FigPie_StrainCounts)
dev.off()

##############################

#sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS  10.16
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base
# other attached packages:
#  [1] ggbiplot_0.55      devtools_2.4.5     usethis_2.1.6      digest_0.6.31      dada2_1.16.0       Rcpp_1.0.10        gplots_3.1.3       RColorBrewer_1.1-3
#  [9] gridExtra_2.3      ggpubr_0.6.0       broom_1.0.3        ggthemes_4.2.4     cowplot_1.1.1      forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0
# [17] purrr_1.0.1        readr_1.3.1        tidyr_1.3.0        tibble_3.1.8       tidyverse_1.3.0    ggplot2_3.4.1      scales_1.2.1       plyr_1.8.8
# [25] lattice_0.20-45    reshape2_1.4.4
# loaded via a namespace (and not attached):
#   [1] colorspace_2.1-0            ggsignif_0.6.4              deldir_1.0-6                hwriter_1.3.2.1             ellipsis_0.3.2
#   [6] XVector_0.28.0              GenomicRanges_1.40.0        fs_1.6.1                    rstudioapi_0.14             remotes_2.4.2
#  [11] fansi_1.0.4                 lubridate_1.9.2             xml2_1.3.3                  cachem_1.0.6                pkgload_1.3.2
#  [16] jsonlite_1.8.4              Rsamtools_2.4.0             dbplyr_2.3.0                png_0.1-8                   shiny_1.7.4
#  [21] compiler_4.0.2              httr_1.4.4                  backports_1.4.1             assertthat_0.2.1            Matrix_1.2-18
#  [26] fastmap_1.1.0               cli_3.6.0                   later_1.3.0                 prettyunits_1.1.1           htmltools_0.5.4
#  [31] tools_4.0.2                 gtable_0.3.1                glue_1.6.2                  GenomeInfoDbData_1.2.3      tinytex_0.44
#  [36] ShortRead_1.46.0            carData_3.0-5               Biobase_2.48.0              cellranger_1.1.0            vctrs_0.5.2
#  [41] Biostrings_2.56.0           xfun_0.37                   ps_1.7.2                    rvest_1.0.3                 mime_0.12
#  [46] timechange_0.2.0            miniUI_0.1.1.1              lifecycle_1.0.3             gtools_3.9.4                rstatix_0.7.2
#  [51] zlibbioc_1.34.0             promises_1.2.0.1            hms_1.1.2                   parallel_4.0.2              SummarizedExperiment_1.18.2
#  [56] memoise_2.0.1               latticeExtra_0.6-30         stringi_1.7.12              S4Vectors_0.26.1            caTools_1.18.2
#  [61] BiocGenerics_0.34.0         pkgbuild_1.4.0              BiocParallel_1.22.0         GenomeInfoDb_1.24.2         rlang_1.0.6
#  [66] pkgconfig_2.0.3             bitops_1.0-7                matrixStats_0.63.0          htmlwidgets_1.6.1           GenomicAlignments_1.24.0
#  [71] processx_3.8.0              tidyselect_1.2.0            magrittr_2.0.3              R6_2.5.1                    profvis_0.3.7
#  [76] IRanges_2.22.2              generics_0.1.3              DelayedArray_0.14.1         DBI_1.1.3                   pillar_1.8.1
#  [81] haven_2.5.1                 withr_2.5.0                 abind_1.4-5                 RCurl_1.98-1.10             modelr_0.1.10
#  [86] crayon_1.5.2                car_3.1-1                   interp_1.0-33               KernSmooth_2.23-20          utf8_1.2.3
#  [91] urlchecker_1.0.1            jpeg_0.1-8.1                readxl_1.3.1                callr_3.7.3                 reprex_2.0.2
#  [96] xtable_1.8-4                httpuv_1.6.8                RcppParallel_5.1.6          stats4_4.0.2                munsell_0.5.0
# [101] sessioninfo_1.2.2
