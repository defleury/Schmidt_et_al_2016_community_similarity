#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Prepare community similarity data for the HMP dataset
#
#=>	load HMP data (V35) for all samples, de novo 97% AL OTUs
#=>	calculate global cooccurrence network
#=>	estimate community similarity using either traditional or interaction-adjusted indices
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
#Read & pre-process data
########################
#Load raw OTU table
ot.raw <- Matrix(as.matrix(read.table(file=PARAM$file.otu_table, header=F, skip=1, sep="\t", row.names=1)), sparse=T);
colnames(ot.raw) <- paste("SN_", samples.raw, sep="");
size.sample.raw <- colSums(ot.raw);
size.otu.raw <- rowSums(ot.raw);
#Load raw OTU data
otu.data.raw <- read.table(file=PARAM$file.otu_data, header=T, sep="\t", row.names=1);
colnames(otu.data.raw) <- c("size", colnames(otu.data.raw)[2:ncol(otu.data.raw)]);
########################

########################
#Filter data
#=> by a minimum of "used" sequences per sample (relative to the raw sequence count per sample)
#=> by minimum sample size
#=> by minimum OTU size (across all samples)
########################
#Minimum relative nuber of retained sequences & minimum sample size
retained.tmp <- size.sample.raw / sample.data.raw[names(size.sample.raw), "Raw_v35"];
samples <- colnames(ot.raw)[size.sample.raw > PARAM$thresh.sample_size & retained.tmp > PARAM$thresh.retained_seqs];
ot.tmp <- ot.raw[, samples];
#Minimum OTU size
otus <- rownames(ot.tmp)[rowSums(ot.tmp) > PARAM$thresh.otu_size];
#Prune data accordingly
ot <- ot.raw[otus, samples];
size.sample <- colSums(ot); n.sample <- length(size.sample);
size.otu <- rowSums(ot); n.otu <- length(size.otu);
otu.data <- otu.data.raw[rownames(ot), ];
otu.data$size <- size.otu;
ot.rel <- t(t(ot) / size.sample);
sample.data <- sample.data.raw[samples, ];
n.otu <- nrow(ot);
n.sample <- length(samples);
########################

########################
#Coerce into phyloseq object
tmp.tree <- read_tree(PARAM$file.tree);
my.tree <- drop.tip(tmp.tree, rownames(ot.raw)[ !rownames(ot.raw) %in% rownames(ot)]);
my.tree.rooted <- root(my.tree, sample(rownames(ot), 1), resolve.root=T);
my.ps <- phyloseq(otu_table(as.matrix(ot), taxa_are_rows=T), sample_data(sample.data), phy_tree(my.tree.rooted));
#Save
save(ot, otu.data, sample.data, my.ps, PARAM, n.otu, n.sample, samples, file=paste(PARAM$folder.data, "hmp_samples.otu_table.RData", sep=""));
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate pairwise indices on OTUs, for count correction
#=> calculate "raw" pairwise similarities between OTUs (SparCC, cophenetic phylogenetic distance)
#=> for each OTU-wise matrix, get the transformed OTU association matrix S (or "C", "Phi"), based on OTU-OTU-correlation across the raw matrix
#=> in other words, pairwise OTU association is calculated as the similarity of OTUs in the original similarity space
#=> matrix S is a symmetric, transformed correlation matrix: S[i,j] = 0.5 * (1 + rho) => S scales on [0,1]
################################################################################
#SparCC correlation
if ("SparCC" == "SparCC") {
	########################
	#SparCC correlation analysis
	#=> following Friedman & Alm, PLOS CB, 2012
	#=> apparently, the pseudocount has a large influence on the results...
	#=> expects an OTU table in which OTUs are rows, samples are columns
	########################
	#Set parameters
	size.thresh <- 1;
	pseudocount <- 10^-6;
	nblocks <- 400;
	use.cores <- 40;
	########################

	########################
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	#Filter OTU table by removing all OTUs that are observed less then <size.thresh> times across all samples
	#=> their correlations to all other OTUs will be (manually) set to 0 later on
	#Add pseudocount to all observed counts (to avoid issues with zeroes in log-space)
	cat(paste("SparCC preallocation steps =>", Sys.time()), sep="\n");
	remove.otus <- rowSums(ot) < size.thresh; keep.otus <- ! remove.otus;
	#o.t <- otu.table[1:1000, ] + pseudocount;
	o.t <- ot[! remove.otus, ] + pseudocount;
	otus <- rownames(o.t);
	n.otu <- length(otus);
	########################
	#Preallocate blocks for parallel processing & Aitchinson's T matrix
	#=> based on https://gist.github.com/bobthecat/5024079
	size.split <- floor(n.otu / nblocks);
	if (size.split < 1) {size.split <- 1}
	my.split <- list(); length(my.split) <- nblocks;
	my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
	my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu;
	dat.split <- mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
	#Get combinations of splits
	my.combs <- expand.grid(1:length(my.split), 1:length(my.split));
	my.combs <- t(apply(my.combs, 1, sort));
	my.combs <- unique(my.combs);
	#Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
	mat.T <- big.matrix(nrow=n.otu, ncol=n.otu, dimnames=list(otus, otus), shared=T);
	mat.T.desc <- describe(mat.T);
	cat(paste("Done with preallocations =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Compute Aitchinson's T matrix
	#=> iterate through each block combination, calculate matrix
	#		between blocks and store them in the preallocated matrix on both
	#		symmetric sides of the diagonal
	cat(paste("Starting parallel T matrix calculation =>", Sys.time()), sep="\n");
	results <- foreach(i = 1:nrow(my.combs)) %dopar% {
		#Get current combination
		curr.comb <- my.combs[i, ];
		#Get current data
		g.1 <- my.split[[curr.comb[1]]];
		g.2 <- my.split[[curr.comb[2]]];
		dat.1 <- dat.split[[curr.comb[1]]];
		dat.2 <- dat.split[[curr.comb[2]]];
		#Get current part of Aitchinson's matrix
		curr.T <- apply(dat.1, 1, function(x) {apply(dat.2, 1, function(y) {var(log(x/y))})});
		#Store
		curr.mat.T <- attach.big.matrix(mat.T.desc);
		curr.mat.T[g.2, g.1] <- curr.T;
		curr.mat.T[g.1, g.2] <- t(curr.T);
		#Return
		TRUE
	}
	cat(paste("Done with parallel T matrix calculation =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Compute component variations ("t_i")
	cat(paste("Computing component variations =>", Sys.time()), sep="\n");
	var.t <- colsum(mat.T);
	cat(paste("Done with component variation calculations =>", Sys.time()), sep="\n");
	#Estimate component variances ("omega_i") from t_i by solving a linear equation system
	cat(paste("Estimating component variances from linear equation system =>", Sys.time()), sep="\n");
	mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu); diag(mat.a) <- n.otu-1;
	omega <- sqrt(solve(a=mat.a, b=var.t));
	cat(paste("Done with component variation estimation =>", Sys.time()), sep="\n");
	#Estimate pairwise correlations based on these values
	cat(paste("Estimating correlations =>", Sys.time()), sep="\n");
	global.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
	rownames(global.sparcc) <- colnames(global.sparcc) <- otus;
	cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Plot histogram
	curr.data <- global.sparcc[upper.tri(global.sparcc, diag=F)];
	pdf(file = paste(PARAM$folder.output, "histogram.SparCC.pdf", sep = ""), width=10, height=10);
	ggplot(data.frame(data=sample(curr.data, 100000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise SparCC Correlations") + ylab("Density");
	dev.off();
	########################
	
	########################
	#Tidy up
	save(global.sparcc, file=paste(PARAM$folder.output, "hmp_samples.SparCC_global.RData", sep=""));
	rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
	########################
}
########################

############################
#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
cat(paste("Done =>", Sys.time()), sep="\n");
S.sparcc <- 0.5 * (tmp.S + 1);
#Plot histogram
curr.data <- S.sparcc[upper.tri(S.sparcc, diag=F)];
pdf(file = paste(PARAM$folder.output, "histogram.S_SparCC.pdf", sep = ""), width=10, height=10);
ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise SparCC Correlations") + ylab("Density");
dev.off();
#Save correlations & tidy up
save(S.sparcc, file=paste(PARAM$folder.data, "hmp_samples.S_SparCC_global.RData", sep=""));
rm(tmp.S, curr.data);
############################

############################
#Get cophenetic distance matrix for current tree
#=> based on the "cophenetic.phylo" function from the ape package
my.tree <- phy_tree(my.ps);
my.tree$node.label <- c("Root", paste("Node", 2:my.tree$Nnode, sep = "_"));
global.cophenetic_tree <- fast.cophenetic.phylo(my.tree, use.cores=40);
#Reorder to match order in OTU table
global.cophenetic_tree <- global.cophenetic_tree[rownames(ot), rownames(ot)];
#Save
save(global.cophenetic_tree, file=paste(PARAM$folder.data, "hmp_samples.pairwise_cophenetic_tree_global.RData"));
############################
#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
cat(paste("Done =>", Sys.time()), sep="\n");
S.phylo <- 0.5 * (tmp.S + 1);
#Plot histogram
curr.data <- S.phylo[upper.tri(S.phylo, diag=F)];
curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise Phylogenetic Correlations") + ylab("Density");
ggsave(curr.plot, width=10, height=10, filename=paste(PARAM$folder.output, "histogram.S_phylo.pdf", sep = ""), useDingbats=F);
#Save correlations & tidy up
save(S.phylo, file=paste(PARAM$folder.data, "hmp_samples.S_phylo_global.RData", sep=""));
rm(tmp.S, curr.data);
############################
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate pairwise community similarities between all HMP samples
#=> classical indices (un-corrected)
#=> corrected indices
################################################################################
#Preallocate template sample similarity matrix
tmp.cs <- matrix(data=NA, nrow=n.sample, ncol=n.sample);
rownames(tmp.cs) <- colnames(tmp.cs) <- samples;

############################
#Preallocate classical indices to calculate
get.cs <- list(); t <- 1;
#Jaccard index, classical, unweighted
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical";
get.cs[[t]]$call <- "jaccard";
t <- t + 1;
#Jaccard index, classical, weighted
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical, weighted, no fractions";
get.cs[[t]]$call <- "jaccard.abd";
t <- t + 1;
#Jaccard index, classical, weighted, fraction-based
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical, weighted";
get.cs[[t]]$call <- "jaccard.abd.frac";
t <- t + 1;
#Jaccard index, classical, weighted, fraction-based, Chao's version
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical, Chao fraction-based";
get.cs[[t]]$call <- "jaccard.abd.frac_chao";
t <- t + 1;
#Jaccard index, classical, weighted, alternative formula
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical, alternative formula";
get.cs[[t]]$call <- "jaccard.abd.alt";
t <- t + 1;
#Jaccard index, classical, Chao' abundance correction
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Jaccard, classical, Chao";
get.cs[[t]]$call <- "jaccard.abd.chao";
t <- t + 1;
#Bray-Curtis, classical
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Bray-Curtis";
get.cs[[t]]$call <- "bray_curtis";
t <- t + 1;
#Morisita-Horn, classical
get.cs[[t]] <- list();
get.cs[[t]]$name <- "Morisita-Horn";
get.cs[[t]]$call <- "morisita_horn";
t <- t + 1;
############################


############################
#Iterate through (classical) indices and calculate pairwise community similarities
for (t in seq(1, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise distances
	#=> as "1-similarity"
	my.cs <- community.similarity.par(ot, distance=get.cs[[t]]$call, use.cores=PARAM$use.cores);
	#Manually correct for rounding errors
	my.cs[my.cs < 0] <- 0;
	#Export histogram
	curr.data <- my.cs[upper.tri(my.cs)];
	pdf(file = paste(PARAM$folder.output, "hmp_samples.similarity_hist.", get.cs[[t]]$name, ".pdf", sep = ""), width=10, height=10);
	curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + xlim(0,1) + ggtitle(paste("Distribution of Pairwise", get.cs[[t]]$name, "Distances")) + ylab("Density");
	print(curr.plot);
	dev.off();
	#Store current similarities
	get.cs[[t]]$cs <- my.cs;
	#Report
	print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
}
############################

############################
#Calculate UniFrac distances / similarities
#=> modified from https://github.com/joey711/phyloseq/blob/6eeb569025d330c5b1b709c103075b37b2feeaff/R/distance-methods.R => "fastUniFrac"
#=> see also https://github.com/joey711/phyloseq/issues/524
########################
#Unweighted UniFrac
t <- 9;
get.cs[[t]] <- list();
get.cs[[t]]$name <- "UniFrac, unweighted";
get.cs[[t]]$call <- "UniFrac_unweighted";
print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
get.cs[[t]]$cs <- fasterUniFrac(ot, my.tree, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
t <- t + 1;
#Weighted UniFrac
get.cs[[t]] <- list();
get.cs[[t]]$name <- "UniFrac, weighted";
get.cs[[t]]$call <- "UniFrac_weighted";
print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
get.cs[[t]]$cs <- fasterUniFrac(ot, my.tree, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
t <- t + 1;
########################
#Plot histograms
########################
#Unweighted UniFrac
curr.data <- as.dist(get.cs[[9]]$cs);
pdf(file = paste(PARAM$folder.output, "hmp_samples.similarity_hist.UniFrac_unweighted.pdf", sep = ""), width=10, height=10);
curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + xlim(0,1) + ggtitle("Distribution of Pairwise unweighted UniFrac Distances") + ylab("Density");
print(curr.plot);
dev.off();
#Weighted UniFrac
curr.data <- as.dist(get.cs[[10]]$cs);
pdf(file = paste(PARAM$folder.output, "hmp_samples.similarity_hist.UniFrac_weighted.pdf", sep = ""), width=10, height=10);
curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + xlim(0,1) + ggtitle("Distribution of Pairwise weighted UniFrac Distances") + ylab("Density");
print(curr.plot);
dev.off();
########################


############################
#Preallocate corrected indices to calculate
t <- start.t <- length(get.cs) + 1;
#Jaccard index, SparCC-corrected, unweighted, normalized
#=> unweighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected, weighted, normalized
#=> weighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, unweighted, normalized
#=> unweighted PINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, weighted, normalized
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
############################

############################
#Iterate through (corrected) indices and calculate pairwise community similarities
for (t in seq(start.t, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise similarities
	if (get.cs[[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
		curr.cs <- community.similarity.corr.par(ot, S=S.sparcc, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
	} else if (get.cs[[t]]$name %in% c("PINA, unweighted", "PINA, weighted")) {
		curr.cs <- community.similarity.corr.par(ot, S=S.phylo, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
	}
	
	#Correct for rounding errors
	curr.cs[curr.cs < 0] <- 0;
	#Export histogram
	curr.data <- curr.cs[upper.tri(curr.cs)];
	curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + xlim(0,1) + ggtitle(paste("Distribution of Pairwise", get.cs[[t]]$name, "Distances")) + ylab("Density");
	ggsave(curr.plot, width=10, height=10, filename=paste(PARAM$folder.output, "hmp_samples.similarity_hist.", get.cs[[t]]$name, ".pdf", sep = ""), useDingbats=F);
	#Store current similarities
	get.cs[[t]]$cs <- curr.cs;
	#Tidy up
	rm(curr.data, curr.cs);
	print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
	gc();
}
############################

#Save pairwise community similarities across methods
all.methods <- sapply(get.cs, function(i) {i$name});
all.calls <- sapply(get.cs, function(i) {i$call});
save(get.cs, file=paste(PARAM$folder.data, "hmp_samples.beta_div.get_cs.RData", sep = ""));
################################################################################
################################################################################

q()




