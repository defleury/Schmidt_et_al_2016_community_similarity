#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Prepare analyses shown in Figures 3 & 4
#
#=>	load HMP data (V35) for all samples, de novo 97% AL OTUs
#=> load community similarities (previously calculated in "prepare.community_similarity_data.R")
#=> plot networks and distance histograms for selected, anecdotal sample pairs
#
#
#2015-08-31
#sebastian.schmidt@imls.uzh.ch
################################################################################

################################################################################
################################################################################
# Load Packages
library("foreach")
library("doMC")
library("iterators")
library("doParallel")
library("parallel")
library("Matrix")
library("bigmemory")
library("biganalytics")
library("gRbase")
library("gplots")
library("ggplot2")
library("grid")
library("gridExtra")
library("data.table")
library("plyr")
library("ape")
library("phyloseq")
library("vegan")
library("RColorBrewer")
library("igraph")
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.input <- #PUT INPUT FOLDER HERE
	PARAM$folder.data <- #PUT DATA FOLDER HERE
	PARAM$folder.output <- #PUT RESULTS/OUTPUT FOLDER HERE
	PARAM$file.functions <- #PUT FUNCTIONS TO SOURCE HERE ("funtions.communtiy_similarity.R")
	PARAM$file.sample_list <- #PUT RAW LIST OF SAMPLES FILE HERE
	PARAM$file.sample_metadata <- #PUT RAW SAMPLE METADATA FILE HERE
	PARAM$file.otu_table <- #PUT OTU TABLE FILE HERE
	PARAM$file.tree <- #PUT OTU PHYLO TREE FILE HERE
	PARAM$file.otu_data <- #PUT OTU DATA FILE HERE
	PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$sample.steps <- c(0.01, 0.02, 0.05, 0.1, 0.15, seq(0.2, 0.9, by=0.1));
PARAM$use.cores <- 40;
###################################

###################################
#Import functions
source(PARAM$file.functions)
###################################

###################################
#Set parameters for data processing
###################################
#Ordered factor of body subsite names
order.body_subsite <- factor(c(
	"Stool",
	"Anterior_nares",
	"Left_Antecubital_fossa", "Right_Antecubital_fossa", "Left_Retroauricular_crease", "Right_Retroauricular_crease",
	"Vaginal_introitus", "Mid_vagina", "Posterior_fornix",
	"Saliva", "Subgingival_plaque", "Supragingival_plaque", "Tongue_dorsum", "Buccal_mucosa", "Hard_palate", "Attached_Keratinized_gingiva", "Throat", "Palatine_Tonsils"
));
#Minimum OTU size required across all samples
PARAM$thresh.otu_size <- 3;
#Minimum sample size required
PARAM$thresh.sample_size <- 1000;
#Minimum relative number of sequences retained
PARAM$thresh.retained_seqs <- 0.5;
#Get list of colors to use for body sites
#=> based on the R Color Brewer qualitative palette "Accent"
PARAM$colors.body_site <- c("#33a02c", "#bf5b17", "#386cb0", "#fdc086", "#beaed4");
names(PARAM$colors.body_site) <- c("Airways", "Gastrointestinal_tract", "Oral", "Skin", "Urogenital_tract");
###################################

###################################
#Read list of samples to process
samples.raw <- as.vector(read.delim(PARAM$file.sample_list, header=F, sep="\n"))$V1;
n.sample.raw <- length(PARAM$use.samples);
###################################
#Read sample raw metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, sep="\t");
rownames(sample.data.raw) <- paste("SN_", as.character(sample.data.raw$SN), sep="");
################################################################################
################################################################################


################################################################################
################################################################################
#Load data
load(file=paste(PARAM$folder.data, "hmp_samples.otu_table.RData", sep=""));
load(file=paste(PARAM$folder.data, "hmp_samples.beta_div.get_cs.RData", sep = ""));

#Get lists of communtiy similarity methods
all.methods <- sapply(get.cs, function(i) {i$name});
all.calls <- sapply(get.cs, function(i) {i$call});
################################################################################
################################################################################


################################################################################
################################################################################
#Generate analyses for "anecdotal" networks
#=> identify interesting anecdotal sets of samples
#=> generate down-sampled cooccurrence network containing anecdotal data plus context
#=> export graphs for Illustrator post-processing
################################################################################
################################################################################
#Extract distances in Jaccard and TINA space for all pairs of samples
combn.samples <- t(combn(samples, 2));
combn.jac <- get.cs[[1]]$cs[combn.samples];
combn.tina <- get.cs[[12]]$cs[combn.samples];
#Prepare sample distance histograms
samples.rand.pick <- sample(seq(1, nrow(combn.samples)), 1000000);
samples.rand.dist_frame <- data.frame(Jaccard=combn.jac, tina=combn.tina);
dens.jac <- ggplot(samples.rand.dist_frame, aes(x=Jaccard)) + geom_density(fill="#8dd3c7", alpha=0.7) + xlim(c(0,1)) + theme_bw();
dens.tina <- ggplot(samples.rand.dist_frame, aes(x=tina)) + geom_density(fill="#80b1d3", alpha=0.7) + xlim(c(0,1)) + theme_bw();
############################
#Case 1: Jaccard and TINA agree well
############################
#Find pairs which are very similar in Jaccard and TINA index
cand.case_1 <- which(combn.jac < 0.5 & combn.tina < 0.05);
#Settle for the pair with minimum Jaccard distance (most shared OTUs)
anecdotes.case_1 <- combn.samples[which.min(combn.jac), ];
otus.case_1 <- rownames(ot)[which(rowSums(ot[, anecdotes.case_1]) > 0)];
otus.case_1.s_1 <- rownames(ot)[which(ot[, anecdotes.case_1[1]] > 0)];
otus.case_1.s_2 <- rownames(ot)[which(ot[, anecdotes.case_1[2]] > 0)];
otus.case_1.shared <- intersect(otus.case_1.s_1, otus.case_1.s_2);
otus.case_1.s_1.excl <- setdiff(otus.case_1.s_1, otus.case_1.shared);
otus.case_1.s_2.excl <- setdiff(otus.case_1.s_2, otus.case_1.shared);
#Export distance histograms with current sample pair marked
curr.plot <- dens.jac + geom_vline(xintercept=get.cs[[1]]$cs[anecdotes.case_1[1], anecdotes.case_1[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_1.Jaccard.pdf", sep = ""), useDingbats=F);
curr.plot <- dens.tina + geom_vline(xintercept=get.cs[[12]]$cs[anecdotes.case_1[1], anecdotes.case_1[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_1.TINA.pdf", sep = ""), useDingbats=F);
############################
#Case 2: No shared OTUs (Jaccard distance == 1), but very high similarity in TINA
############################
#Find (pairs of) samples which share little or no OTUs (high Jaccard distance) but are highly similar in their interaction context (low TINA distance)
cand.case_2 <- which(combn.jac == 1 & combn.tina < 0.03);
#Pick a suitable pair from these
anecdotes.case_2 <- combn.samples[cand.case_2[2],]
otus.case_2 <- rownames(ot)[which(rowSums(ot[, anecdotes.case_2]) > 0)];
otus.case_2.s_1 <- rownames(ot)[which(ot[, anecdotes.case_2[1]] > 0)];
otus.case_2.s_2 <- rownames(ot)[which(ot[, anecdotes.case_2[2]] > 0)]
otus.case_2.shared <- intersect(otus.case_2.s_1, otus.case_2.s_2);
otus.case_2.s_1.excl <- setdiff(otus.case_2.s_1, otus.case_2.shared);
otus.case_2.s_2.excl <- setdiff(otus.case_2.s_2, otus.case_2.shared);
#Export distance histograms with current sample pair marked
curr.plot <- dens.jac + geom_vline(xintercept=get.cs[[1]]$cs[anecdotes.case_2[1], anecdotes.case_2[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_2.Jaccard.pdf", sep = ""), useDingbats=F);
curr.plot <- dens.tina + geom_vline(xintercept=get.cs[[12]]$cs[anecdotes.case_2[1], anecdotes.case_2[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_2.TINA.pdf", sep = ""), useDingbats=F);
############################
#Case 3: Relatively close in Jaccard space, but very high dissimilarity in TINA
############################
#Rank samples by Jaccard and TINA distances
ranks.jac <- rank(combn.jac);
ranks.tina <- rank(combn.tina);
#Find greatest discrepancy in ranks (high rank in TINA, low rank in Jaccard)
ranks.diff <- ranks.tina - ranks.jac;
cand.case_3 <- which.max(ranks.diff);
anecdotes.case_3 <- combn.samples[cand.case_3, ];
otus.case_3 <- rownames(ot)[which(rowSums(ot[, anecdotes.case_3]) > 0)];
otus.case_3.s_1 <- rownames(ot)[which(ot[, anecdotes.case_3[1]] > 0)];
otus.case_3.s_2 <- rownames(ot)[which(ot[, anecdotes.case_3[2]] > 0)]
otus.case_3.shared <- intersect(otus.case_3.s_1, otus.case_3.s_2);
otus.case_3.s_1.excl <- setdiff(otus.case_3.s_1, otus.case_3.shared);
otus.case_3.s_2.excl <- setdiff(otus.case_3.s_2, otus.case_3.shared);
#Export distance histograms with current sample pair marked
curr.plot <- dens.jac + geom_vline(xintercept=get.cs[[1]]$cs[anecdotes.case_3[1], anecdotes.case_3[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_3.Jaccard.pdf", sep = ""), useDingbats=F);
curr.plot <- dens.tina + geom_vline(xintercept=get.cs[[12]]$cs[anecdotes.case_3[1], anecdotes.case_3[2]]);
ggsave(plot=curr.plot, width=20, height=10, filename=paste(PARAM$folder.output, "hmp_samples.anecdotes.similarity_hist.case_3.TINA.pdf", sep = ""), useDingbats=F);
############################
#Plot larger network with context
############################
library("igraph");
#Add some random samples per body site to fill up network for plotting
network.samples <- c(
	anecdotes.case_1, anecdotes.case_2, anecdotes.case_3,
	sample(samples[sample.data$Body_Site == "Gastrointestinal_tract"], 2),
	sample(samples[sample.data$Body_Site == "Oral"], 2),
	sample(samples[sample.data$Body_Site == "Skin"], 2),
	sample(samples[sample.data$Body_Site == "Airways"], 2),
	sample(samples[sample.data$Body_Site == "Urogenital_tract"], 2)
);
network.samples <- unique(network.samples);
#Get current pruned OTU table
curr.ot <- ot[which(rowSums(ot[,network.samples]) > 0), network.samples];
curr.otus <- rownames(curr.ot);
############################
#Plot SparCC and derived S.SparCC networks as (subsampled) matrices as heatmaps
curr.otus.s <- sample(curr.otus, 300);
curr.S <- 2 * S.sparcc[curr.otus.s, curr.otus.s] - 1;
curr.sparcc <- global.sparcc[curr.otus.s, curr.otus.s]
#Cluster by S.sparcc
curr.clust.S <- hclust(as.dist(1-S.sparcc[curr.otus.s, curr.otus.s]));
#Plot heatmap of SparCC cooccurrences
pdf(file = paste(PARAM$folder.output, "hmp_samples.example_network.heatmap.SparCC.pdf", sep = ""), width=20, height=20, useDingbats=F);
heatmap.2(
	curr.sparcc,
	Rowv = as.dendrogram(curr.clust.S), Colv = as.dendrogram(curr.clust.S),
	dendrogram = "none",
	na.rm = T,
	revC = T,
	col=colorRampPalette(rev(brewer.pal(9, "PRGn")))(200),
	breaks = seq(-1, 1, length.out=201),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
#Plot heatmap of S.SparCC
pdf(file = paste(PARAM$folder.output, "hmp_samples.example_network.heatmap.S_SparCC.pdf", sep = ""), width=20, height=20, useDingbats=F);
heatmap.2(
	curr.S,
	Rowv = as.dendrogram(curr.clust.S), Colv = as.dendrogram(curr.clust.S),
	dendrogram = "none",
	na.rm = T,
	revC = T,
	col=colorRampPalette(rev(brewer.pal(9, "PRGn")))(200),
	breaks = seq(-1, 1, length.out=201),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
############################
#Generate network plots
############################
#Get SparCC cooc between current OTUs
#curr.sparcc <- S.sparcc[curr.otus, curr.otus];
#curr.sparcc[curr.sparcc < 0.978] <- 0;
curr.sparcc <- global.sparcc[curr.otus, curr.otus];
curr.sparcc[curr.sparcc < 0.5] <- 0;
#Coerce into igraph object
curr.graph <- graph.adjacency(curr.sparcc, mode="upper", weighted=T, diag=F);
#Remove unconnected vertices
remove.vertices <- which(degree(curr.graph) == 0);
keep.vertices <- which(degree(curr.graph) > 0);
keep.vertices <- seq(1, vcount(curr.graph));
#curr.graph <- delete.vertices(curr.graph, remove.vertices);
#Get some stats and metadata on current graph
curr.graph.dat <- list();
curr.graph.dat$size <- rowSums(ot[curr.otus[keep.vertices],]);
#Color by dominant body site
curr.graph.dat$color <- rep.int("#bdbdbd", length(curr.otus));
curr.graph.dat$color[otu.data[curr.otus[keep.vertices], "Dominant_Body_Site"] == "Gastrointestinal_tract"] <- "#fdc086";
curr.graph.dat$color[otu.data[curr.otus[keep.vertices], "Dominant_Body_Site"] == "Oral"] <- "#386cb0";
curr.graph.dat$color[otu.data[curr.otus[keep.vertices], "Dominant_Body_Site"] == "Skin"] <- "#ffff99";
curr.graph.dat$color[otu.data[curr.otus[keep.vertices], "Dominant_Body_Site"] == "Airways"] <- "#7fc97f";
curr.graph.dat$color[otu.data[curr.otus[keep.vertices], "Dominant_Body_Site"] == "Urogenital_tract"] <- "#beaed4";
#Color by taxonomy
curr.graph.dat$color.tax <- rep.int("#bdbdbd", length(curr.otus));
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "consensus_phylum"] == "Firmicutes"] <- "#a6cee3";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "consensus_phylum"] == "Bacteroidetes"] <- "#fb9a99";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "consensus_phylum"] == "Proteobacteria"] <- "#33a02c";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "consensus_phylum"] == "Actinobacteria"] <- "#b2df8a";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "consensus_phylum"] == "Fusobacteria"] <- "#1f78b4";
############################
#Mark OTUs involved in anecdotal cases
#Case 1
#Background
curr.graph.dat$color.case_1 <- rep.int("#f0f0f0", length(curr.otus));
#OTUs exclusive to sample 1
curr.graph.dat$color.case_1[curr.otus %in% otus.case_1.s_1.excl] <- "#a1dab4"; #fcc5c0; #a1dab4 #74c476
#OTUs exclusive to sample 2
curr.graph.dat$color.case_1[curr.otus %in% otus.case_1.s_2.excl] <- "#225ea8"; #ae017e;
#OTUs shared between samples
curr.graph.dat$color.case_1[curr.otus %in% otus.case_1.shared] <- "#41b6c4"; #beaed4;
#Remove frame.colors for non-focal OTUs
curr.graph.dat$frame_col.case_1 <- rep.int("grey", length(curr.otus));
curr.graph.dat$frame_col.case_1[curr.otus %in% otus.case_1] <- "black";
#Get color by taxonomy
curr.graph.dat$color.case_1.tax <- rep.int("#f0f0f0", length(curr.otus));
curr.graph.dat$color.case_1.tax[curr.otus %in% otus.case_1] <- curr.graph.dat$color.tax[curr.otus %in% otus.case_1];
#Case 2
curr.graph.dat$color.case_2 <- rep.int("#f0f0f0", length(curr.otus));
curr.graph.dat$color.case_2[curr.otus %in% otus.case_2.s_1.excl] <- "#a1dab4";
curr.graph.dat$color.case_2[curr.otus %in% otus.case_2.s_2.excl] <- "#225ea8";
curr.graph.dat$color.case_2[curr.otus %in% otus.case_2.shared] <- "#41b6c4";
curr.graph.dat$frame_col.case_2 <- rep.int("grey", length(curr.otus));
curr.graph.dat$frame_col.case_2[curr.otus %in% otus.case_2] <- "black";
curr.graph.dat$color.case_2.tax <- rep.int("#f0f0f0", length(curr.otus));
curr.graph.dat$color.case_2.tax[curr.otus %in% otus.case_2] <- curr.graph.dat$color.tax[curr.otus %in% otus.case_2];
#Case 3
curr.graph.dat$color.case_3 <- rep.int("#f0f0f0", length(curr.otus));
curr.graph.dat$color.case_3[curr.otus %in% otus.case_3.s_1.excl] <- "#a1dab4"; #ffff99;
curr.graph.dat$color.case_3[curr.otus %in% otus.case_3.s_2.excl] <- "#225ea8"; #386cb0;
curr.graph.dat$color.case_3[curr.otus %in% otus.case_3.shared] <- "#41b6c4"; #41ae76;
curr.graph.dat$frame_col.case_3 <- rep.int("grey", length(curr.otus));
curr.graph.dat$frame_col.case_3[curr.otus %in% otus.case_3] <- "black";
curr.graph.dat$color.case_3.tax <- rep.int("#f0f0f0", length(curr.otus));
curr.graph.dat$color.case_3.tax[curr.otus %in% otus.case_3] <- curr.graph.dat$color.tax[curr.otus %in% otus.case_3];
############################
#Compute a layout
curr.layout <- layout.kamada.kawai(curr.graph)
#Export graph, coloured by dominant body site
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by taxonomy
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.color_tax.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.tax, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 1
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_1.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_1, vertex.frame.color=curr.graph.dat$frame_col.case_1, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 1, by taxonomy
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_1.color_tax.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_1.tax, vertex.frame.color=curr.graph.dat$frame_col.case_1, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 2
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_2.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_2, vertex.frame.color=curr.graph.dat$frame_col.case_2, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 2, by taxonomy
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_2.color_tax.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_2.tax, vertex.frame.color=curr.graph.dat$frame_col.case_2, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 3
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_3.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_3, vertex.frame.color=curr.graph.dat$frame_col.case_3, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
#Export graph, coloured by anecdotal case 3, by taxonomy
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.full_network.highlight_case_3.color_tax.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.case_3.tax, vertex.frame.color=curr.graph.dat$frame_col.case_3, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
############################
#Save stuff
save(curr.graph, curr.graph.dat, curr.layout, curr.otus, network.samples, anecdotes.case_1, anecdotes.case_2, anecdotes.case_3, file=paste(PARAM$folder.output, "hmp_samples.anecdotes.data.RData", sep = ""));
############################
#Plot subnetwork for case 1
curr.graph.case_1 <- graph.adjacency(curr.sparcc[otus.case_1, otus.case_1], weighted=T, mode="upper", diag=F);
#Annotate
curr.graph.case_1.dat <- list();
curr.graph.case_1.dat$size <- rowSums(ot[otus.case_1,]);
curr.graph.case_1.dat$color[otus.case_1 %in% otus.case_1.s_1.excl] <- "#a1dab4";
curr.graph.case_1.dat$color[otus.case_1 %in% otus.case_1.s_2.excl] <- "#225ea8";
curr.graph.case_1.dat$color[otus.case_1 %in% otus.case_1.shared] <- "#41b6c4";
#Get layout for current subnetwork
curr.layout.case_1 <- layout.fruchterman.reingold(curr.graph.case_1, maxiter=5000);
#Plot graph
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.sub_network.case_1.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph.case_1, vertex.size=log10(curr.graph.case_1.dat$size), vertex.label=NA, vertex.color=curr.graph.case_1.dat$color, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.3), edge.curved=0.3, layout=curr.layout.case_1);
dev.off();
############################
#Plot subnetwork for case 2
curr.graph.case_2 <- graph.adjacency(curr.sparcc[otus.case_2, otus.case_2], weighted=T, mode="upper", diag=F);
#Annotate
curr.graph.case_2.dat <- list();
curr.graph.case_2.dat$size <- rowSums(ot[otus.case_2,]);
curr.graph.case_2.dat$color[otus.case_2 %in% otus.case_2.s_1.excl] <- "#a1dab4";
curr.graph.case_2.dat$color[otus.case_2 %in% otus.case_2.s_2.excl] <- "#225ea8";
curr.graph.case_2.dat$color[otus.case_2 %in% otus.case_2.shared] <- "#41b6c4";
#Get layout for current subnetwork
curr.layout.case_2 <- layout.fruchterman.reingold(curr.graph.case_2, maxiter=5000);
#Plot graph
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.sub_network.case_2.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph.case_2, vertex.size=log10(curr.graph.case_2.dat$size), vertex.label=NA, vertex.color=curr.graph.case_2.dat$color, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.3), edge.curved=0.3, layout=curr.layout.case_2);
dev.off();
############################
#Plot subnetwork for case 3
curr.graph.case_3 <- graph.adjacency(curr.sparcc[otus.case_3, otus.case_3], weighted=T, mode="upper", diag=F);
#Annotate
curr.graph.case_3.dat <- list();
curr.graph.case_3.dat$size <- rowSums(ot[otus.case_3,]);
curr.graph.case_3.dat$color[otus.case_3 %in% otus.case_3.s_1.excl] <- "#a1dab4";
curr.graph.case_3.dat$color[otus.case_3 %in% otus.case_3.s_2.excl] <- "#225ea8";
curr.graph.case_3.dat$color[otus.case_3 %in% otus.case_3.shared] <- "#41b6c4";
#Get layout for current subnetwork
curr.layout.case_3 <- layout.fruchterman.reingold(curr.graph.case_3, maxiter=5000);
#Plot graph
pdf(file = paste(PARAM$folder.output, "hmp_samples.anecdotes.sub_network.case_3.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph.case_3, vertex.size=log10(curr.graph.case_3.dat$size), vertex.label=NA, vertex.color=curr.graph.case_3.dat$color, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.3), edge.curved=0.3, layout=curr.layout.case_3);
dev.off();
################################################################################
################################################################################



q()