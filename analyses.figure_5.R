#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Prepare analyses shown in Figure 5
#=> down-sampling experiments
#
#=>	load HMP data (V35) for all samples, de novo 97% AL OTUs
#=> load community similarities (previously calculated in "prepare.community_similarity_data.R")
#=> perform downsampling per samples and per sequences (as shown in Figure 5A and 5B)
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
#Robustness to downsampling samples
#=> Figure 5A
#
#=> process body sites and sub-sites (within Oral, Urogenital and Skin) separately
#=> (randomly) select subsets of samples (n=5, 10, 20, etc per body (sub-)site) from the original dataset
#=> reduce the OTU table accordingly, in OTU and sample dimensions
#=> prune the phylogenetic tree accordingly
#=> re-calculate correlation & S matrix per subset
#=> iteratively downsample dataset, 10 times per step, and recalculate diversity indices
#=> restrict downsampling to include equal numbers of samples per body site (unlike full dataset)
#=> compare downsampled diversities to full diversities, in terms of 
#		(i)		distance (distribution) of per-sample diversity to baseline
#		(ii)	mean shift of per-sample diversity to baseline
#		(iii)	variance of shift of per-sample diversity to baseline
#		(iv)	F value of PERMANOVA between body sites, subsites, etc
#		(v)		R2 value of PERMANOVA between body sites, subsites, etc
################################################################################
############################
#Calculate baseline similarities for full dataset(s)
#=> process samples within Oral, Urogenital and Skin separately
#=> re-learn S matrices
#=> calculate diversities
############################
#Preallocate collector lists for baseline community similarities and S matrices
bs.cs <- list(); length(bs.cs) <- 3; names(bs.cs) <- c("Skin", "Urogenital_tract", "Oral");
bs.S <- list(); length(bs.S) <- 3; names(bs.S) <- c("Skin", "Urogenital_tract", "Oral");
bs.adonis <- list(); length(bs.adonis) <- 3; names(bs.adonis) <- c("Skin", "Urogenital_tract", "Oral");
#Iterate through body sites
for (bs in c("Skin", "Urogenital_tract", "Oral")) {
	############################
	#Sub-set OTU table
	bs.samples <- rownames(sample.data)[sample.data$Body_Site == bs];
	tmp.ot <- ot[, bs.samples]; bs.ot <- tmp.ot[rowSums(tmp.ot) > 0, ]; rm(tmp.ot);
	bs.sample.data <- sample.data[bs.samples, ];
	bs.otus <- rownames(bs.ot); bs.otu_sizes <- rowSums(bs.ot); bs.sample_sizes <- colSums(bs.ot);
	#Prune tree
	bs.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[ !rownames(ot) %in% bs.otus]);
	#Report
	cat(paste("Body site", bs, "has", ncol(bs.ot), "samples and", nrow(bs.ot), "OTUs =>", Sys.time()), sep="\n");
	############################
	
	############################
	#Calculate S matrices
	############################
	#SparCC
	if ("SparCC" == "SparCC") {
		########################
		#Set parameters
		size.thresh <- 1;
		pseudocount <- 10^-6;
		nblocks <- 400;
		use.cores <- 20;
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
		remove.otus <- rowSums(bs.ot) < size.thresh; keep.otus <- ! remove.otus;
		#o.t <- otu.table[1:1000, ] + pseudocount;
		o.t <- bs.ot[! remove.otus, ] + pseudocount;
		otus <- rownames(bs.ot);
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
		curr.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
		rownames(curr.sparcc) <- colnames(curr.sparcc) <- otus;
		cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
		########################
		
		########################
		#Tidy up
		rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
		########################
	}
	#Calculate derived OTU similarity matrix S
	cat(paste(bs, "calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	bs.S[[bs]]$S.sparcc <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S, curr.sparcc);
	############################
	#Phylogenetic distances
	#Calculate cophenetic phylogenetic distances
	curr.cophenetic_tree <- cophenetic(bs.tree);
	#Reorder to match order in OTU table
	curr.cophenetic_tree <- curr.cophenetic_tree[rownames(bs.ot), rownames(bs.ot)];
	#Calculate derived OTU similarity matrix S
	cat(paste(bs, "calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	bs.S[[bs]]$S.phylo <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S, curr.cophenetic_tree);
	############################
	
	############################
	#Calculate community similarities
	bs.cs[[bs]] <- list(); t <- 1;
	############################
	#Classical indices
	############################
	#Jaccard index, classical, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical"; bs.cs[[bs]][[t]]$call <- "jaccard"; t <- t + 1;
	#Jaccard index, classical, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, weighted, no fractions"; bs.cs[[bs]][[t]]$call <- "jaccard.abd"; t <- t + 1;
	#Jaccard index, classical, weighted, fraction-based
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
	#Jaccard index, classical, weighted, fraction-based, Chao's version
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, Chao fraction-based"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.frac_chao"; t <- t + 1;
	#Jaccard index, classical, weighted, alternative formula
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, alternative formula"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.alt"; t <- t + 1;
	#Jaccard index, classical, Chao' abundance correction
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, Chao"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
	#Bray-Curtis, classical
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Bray-Curtis"; bs.cs[[bs]][[t]]$call <- "bray_curtis"; t <- t + 1;
	#Morisita-Horn, classical
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Morisita-Horn"; bs.cs[[bs]][[t]]$call <- "morisita_horn"; t <- t + 1;
	#Iterate through indices and calculate pairwise sample similarities
	for (t in seq(1, length(bs.cs[[bs]]))) {
		print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
		#Calculate all pairwise distances
		my.cs <- community.similarity.par(bs.ot, distance= bs.cs[[bs]][[t]]$call, use.cores=PARAM$use.cores);
		#Manually correct for rounding errors
		my.cs[my.cs < 0] <- 0;
		#Store current similarities
		bs.cs[[bs]][[t]]$cs <- my.cs;
		#Report
		print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	}
	############################
	#UniFrac
	############################
	t <- length(bs.cs[[bs]]) + 1;
	#Unweighted UniFrac
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "UniFrac, unweighted"; bs.cs[[bs]][[t]]$call <- "UniFrac_unweighted";
	print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
	bs.cs[[bs]][[t]]$cs <- fasterUniFrac(bs.ot, bs.tree, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
	print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	t <- t + 1;
	#Weighted UniFrac
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "UniFrac, weighted"; bs.cs[[bs]][[t]]$call <- "UniFrac_weighted";
	print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
	bs.cs[[bs]][[t]]$cs <- fasterUniFrac(bs.ot, bs.tree, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
	print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	t <- t + 1;
	############################
	#Corrected indices
	############################
	t <- start.t <- length(bs.cs[[bs]]) + 1;
	#TINA, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "TINA, unweighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
	#TINA, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "TINA, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
	#PINA, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "PINA, unweighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
	#PINA, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "PINA, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
	############################
	#Iterate through (corrected) indices and calculate pairwise community similarities
	for (t in seq(start.t, length(bs.cs[[bs]]))) {
		print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
		#Calculate all pairwise similarities
		if (bs.cs[[bs]][[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
			curr.cs <- community.similarity.corr.par(bs.ot, S=bs.S[[bs]]$S.sparcc, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
		} else {
			curr.cs <- community.similarity.corr.par(bs.ot, S=bs.S[[bs]]$S.phylo, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
		}
		#Manually correct for rounding errors
		curr.cs[curr.cs < 0] <- 0;
		#Store current similarities
		bs.cs[[bs]][[t]]$cs <- curr.cs;
		#Report
		print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
		rm(curr.cs)
	}
	############################
	
	############################
	#Calculate PERMANOVA scores, in parallel
	#=> as a baseline for subsequent downsampling experiments
	############################
	registerDoMC(cores=length(bs.cs[[bs]]));
	bs.adonis[[bs]] <- foreach(t = seq(1, length(bs.cs[[bs]]))) %dopar% {curr.dist <- as.dist(bs.cs[[bs]][[t]]$cs); adonis(curr.dist ~ Body_Subsite, data=bs.sample.data)}
	
	############################
	#Save current progress
	save(bs.S, file=paste(PARAM$folder.data, "hmp_samples.data.body_subsites.S.RData", sep=""))
	save(bs.cs, bs.adonis, file=paste(PARAM$folder.data, "hmp_samples.data.body_subsites.RData", sep=""))
	############################
}
############################

############################
#Perform downsampling experiments on body sites
#=> randomly sub-select n samples per body site (n=5, 10, 20, etc.)
#=> iterate 10 times per downsampling experiment
#=> prune down dataset to include only these samples
#=> re-learn S matrices
#=> calculate community similarities
#=> compare community similarities to baseline
#=> calculate PERMANOVA stats
############################
#Set downsampling parameters
#=> downsampling "steps", i.e. number of samples per body (sub-)site to use
bs.steps <- c(5, 10, 20, 50);
#=> iterations per downsampling step
bs.iterations <- 10;
############################
#Preallocate
all.methods <- unlist(lapply(bs.cs[["Skin"]], function(t) {t$name}));
############################
#Results collector data.frame
bs_ds.res <- data.frame(
	Step=c(rep(ncol(ot), times=4*length(all.methods)), rep(bs.steps, each=length(all.methods) * 4 * bs.iterations)),
	Iteration=c(rep(1, times=4*length(all.methods)), rep(seq(1, bs.iterations), times=length(bs.steps), each=4*length(all.methods))),
	Method=c(rep(all.methods, each=4), rep(all.methods, times=length(bs.steps) * bs.iterations, each=4)),
	Set=c(rep(c("Body_Site", "Skin", "Urogenital_tract", "Oral"), times=length(all.methods)), rep(c("Body_Site", "Skin", "Urogenital_tract", "Oral"), times=length(all.methods) * bs.iterations * 4)),
	F=NA,
	R2=NA,
	P_Value=NA,
	Mean_Log2FC=NA,
	SD_Log2FC=NA
);
#Populate
for (t in seq(1, length(all.methods))) {for (bs in c("Body_Site", "Skin", "Urogenital_tract", "Oral")) {
	if (bs == "Body_Site") {curr.adonis <- get.cs[[t]]$permanova.body_site} else {curr.adonis <- bs.adonis[[bs]][[t]]}
	curr.idx <- which(bs_ds.res$Step == ncol(ot) & bs_ds.res$Iteration == 1 & bs_ds.res$Method == all.methods[t] & bs_ds.res$Set == bs);
	bs_ds.res[curr.idx, "F"] <- curr.adonis$aov.tab$F.Model[1];
	bs_ds.res[curr.idx, "R2"] <- curr.adonis$aov.tab$R2[1];
	bs_ds.res[curr.idx, "P_Value"] <- curr.adonis$aov.tab[1,"Pr(>F)"];
	bs_ds.res[curr.idx, "Mean_Log2FC"] <- 0;
	bs_ds.res[curr.idx, "SD_Log2FC"] <- 0;
}}
############################
#Iterate through downsampling steps
############################
#Body sites
bs <- "Body_Site";
for (step in bs.steps) {
	cat(paste(Sys.time(), "=> Start downsampling with", step, "samples per body site"), sep="\n");
	#Loop through iterations
	for (i in seq(1, bs.iterations)) {
		cat(paste(Sys.time(), "=> Iteration", i), sep="\n");
		############################
		#Downsample OTU table
		curr.samples <- unlist(lapply(levels(sample.data$Body_Site), function(s) {sample(rownames(sample.data)[which(sample.data$Body_Site == s)], size=step, replace=F)}));
		tmp.ot <- ot[,curr.samples]; curr.ot <- tmp.ot[rowSums(tmp.ot) > 0,];
		curr.sample.data <- sample.data[curr.samples, ];
		curr.otus <- rownames(curr.ot); curr.otu_sizes <- rowSums(curr.ot); curr.sample_sizes <- colSums(curr.ot);
		#Prune tree
		curr.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[ !rownames(ot) %in% curr.otus]);
		#Report
		cat(paste("In current iteration, body sites have", ncol(curr.ot), "samples and", nrow(curr.ot), "OTUs =>", Sys.time()), sep="\n");
		############################
		
		############################
		#Calculate S matrices
		############################
		#SparCC
		if ("SparCC" == "SparCC") {
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
			remove.otus <- rowSums(curr.ot) < size.thresh; keep.otus <- ! remove.otus;
			#o.t <- otu.table[1:1000, ] + pseudocount;
			o.t <- curr.ot[! remove.otus, ] + pseudocount;
			otus <- rownames(curr.ot);
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
			curr.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
			rownames(curr.sparcc) <- colnames(curr.sparcc) <- otus;
			cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
			########################
			
			########################
			#Tidy up
			rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
			########################
		}
		#Calculate derived OTU similarity matrix S
		cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
		tmp.S <- cor.par(curr.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
		cat(paste("Done =>", Sys.time()), sep="\n");
		curr.S.sparcc <- 0.5 * (tmp.S + 1);
		#Tidy up
		rm(tmp.S, curr.sparcc);
		############################
		#Phylogenetic distances
		#Calculate cophenetic phylogenetic distances
		curr.cophenetic_tree <- cophenetic(curr.tree);
		#Reorder to match order in OTU table
		curr.cophenetic_tree <- curr.cophenetic_tree[rownames(curr.ot), rownames(curr.ot)];
		#Calculate derived OTU similarity matrix S
		cat(paste("Calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
		tmp.S <- cor.par(curr.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
		cat(paste("Done =>", Sys.time()), sep="\n");
		curr.S.phylo <- 0.5 * (tmp.S + 1);
		#Tidy up
		rm(tmp.S, curr.cophenetic_tree);
		############################
		
		
		############################
		#Calculate community similarities
		#=> although most of the classical indices should be impartial to downsampling, they are still re-calculated, just to be sure
		############################
		############################
		#Calculate community similarities
		ds.cs <- list(); t <- 1;
		############################
		#Classical indices
		############################
		#Jaccard index, classical, unweighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical"; ds.cs[[t]]$call <- "jaccard"; t <- t + 1;
		#Jaccard index, classical, weighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, weighted, no fraction"; ds.cs[[t]]$call <- "jaccard.abd"; t <- t + 1;
		#Jaccard index, classical, weighted, fraction-based
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, weighted"; ds.cs[[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
		#Jaccard index, classical, weighted, fraction-based, Chao's version
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, Chao fraction-based"; ds.cs[[t]]$call <- "jaccard.abd.frac_chao"; t <- t + 1;
		#Jaccard index, classical, weighted, alternative formula
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, alternative formula"; ds.cs[[t]]$call <- "jaccard.abd.alt"; t <- t + 1;
		#Jaccard index, classical, Chao' abundance correction
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, Chao"; ds.cs[[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
		#Bray-Curtis, classical
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Bray-Curtis"; ds.cs[[t]]$call <- "bray_curtis"; t <- t + 1;
		#Morisita-Horn, classical
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Morisita-Horn"; ds.cs[[t]]$call <- "morisita_horn"; t <- t + 1;
		#Iterate through indices and calculate pairwise sample similarities
		for (t in seq(1, length(ds.cs))) {
			print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
			#Calculate all pairwise distances
			my.cs <- community.similarity.par(curr.ot, distance= ds.cs[[t]]$call, use.cores=PARAM$use.cores);
			#Manually correct for rounding errors
			my.cs[my.cs < 0] <- 0;
			#Store current similarities
			ds.cs[[t]]$cs <- my.cs;
			#Report
			print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
		}
		############################
		#UniFrac
		############################
		t <- length(ds.cs) + 1;
		#Unweighted UniFrac
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "UniFrac, unweighted"; ds.cs[[t]]$call <- "UniFrac_unweighted";
		print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
		ds.cs[[t]]$cs <- fasterUniFrac(curr.ot, curr.tree, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
		print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
		t <- t + 1;
		#Weighted UniFrac
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "UniFrac, weighted"; ds.cs[[t]]$call <- "UniFrac_weighted";
		print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
		ds.cs[[t]]$cs <- fasterUniFrac(curr.ot, curr.tree, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
		print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
		t <- t + 1;
		############################
		#Corrected indices
		############################
		t <- start.t <- length(ds.cs) + 1;
		#TINA, unweighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "TINA, unweighted"; ds.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
		#TINA, weighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "TINA, weighted"; ds.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
		#PINA, unweighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "PINA, unweighted"; ds.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
		#PINA, weighted
		ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "PINA, weighted"; ds.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
		############################
		#Iterate through (corrected) indices and calculate pairwise community similarities
		for (t in seq(start.t, length(ds.cs))) {
			print(paste(Sys.time(), bs, "starting", ds.cs[[t]]$name));
			#Calculate all pairwise similarities
			if (ds.cs[[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
				curr.cs <- community.similarity.corr.par(curr.ot, S=curr.S.sparcc, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
			} else {
				curr.cs <- community.similarity.corr.par(curr.ot, S=curr.S.phylo, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
			}
			#Manually correct for rounding errors
			curr.cs[curr.cs < 0] <- 0;
			#Store current similarities
			ds.cs[[t]]$cs <- curr.cs;
			#Report
			print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
			rm(curr.cs)
		}
		############################
		
		############################
		#Calculate PERMANOVA statistics
		############################
		registerDoMC(cores=length(ds.cs));
		curr.adonis <- foreach(t = seq(1, length(ds.cs))) %dopar% {curr.dist <- as.dist(ds.cs[[t]]$cs); adonis(curr.dist ~ Body_Site, data=curr.sample.data)}
		for (t in seq(1, length(all.methods))) {
			curr.idx <- which(bs_ds.res$Step == step & bs_ds.res$Iteration == i & bs_ds.res$Method == all.methods[t] & bs_ds.res$Set == "Body_Site");
			#Store current PERMANOVA stats
			bs_ds.res[curr.idx, "F"] <- curr.adonis[[t]]$aov.tab$F.Model[1];
			bs_ds.res[curr.idx, "R2"] <- curr.adonis[[t]]$aov.tab$R2[1];
			bs_ds.res[curr.idx, "P_Value"] <- curr.adonis[[t]]$aov.tab[1,"Pr(>F)"];
			#Store current (averaged) log2FC of pairwise community distances
			curr.ref_cs <- as.dist(get.cs[[t]]$cs[curr.samples, curr.samples]);
			curr.log2fc <- log2(as.dist(ds.cs[[t]]$cs) / curr.ref_cs);
			bs_ds.res[curr.idx, "Mean_Log2FC"] <- mean(curr.log2fc, na.rm=T);
			bs_ds.res[curr.idx, "SD_Log2FC"] <- sd(curr.log2fc, na.rm=T);
		}
		############################
	}
	#Save current data
	save(bs_ds.res, file=paste(PARAM$folder.data, "hmp_samples.downsampling_samples.RData", sep=""));
}
############################
#Body subsites
for (bs in c("Skin", "Urogenital_tract", "Oral")) {
	#Iterate through downsampling steps
	for (step in bs.steps) {
		cat(paste(Sys.time(), "=> Start downsampling with", step, "samples per", bs, "body subsites"), sep="\n");
		#Loop through iterations
		for (i in seq(1, bs.iterations)) {
			cat(paste(Sys.time(), "=> Iteration", i, bs), sep="\n");
			############################
			#Downsample OTU table
			curr.samples <- unlist(lapply(as.character(unique(sample.data$Body_Subsite[which(sample.data$Body_Site == bs)])), function(s) {sample(rownames(sample.data)[which(sample.data$Body_Subsite == s)], size=step, replace=F)}));
			tmp.ot <- ot[,curr.samples]; curr.ot <- tmp.ot[rowSums(tmp.ot) > 0,];
			curr.sample.data <- sample.data[curr.samples, ];
			curr.otus <- rownames(curr.ot); curr.otu_sizes <- rowSums(curr.ot); curr.sample_sizes <- colSums(curr.ot);
			#Prune tree
			curr.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[ !rownames(ot) %in% curr.otus]);
			#Report
			cat(paste("In current iteration, body sites have", ncol(curr.ot), "samples and", nrow(curr.ot), "OTUs =>", Sys.time()), sep="\n");
			############################
			
			############################
			#Calculate S matrices
			############################
			#SparCC
			if ("SparCC" == "SparCC") {
				########################
				#Set parameters
				size.thresh <- 1;
				pseudocount <- 10^-6;
				if (nrow(curr.ot) <= 800) {nblocks <- 100} else {nblocks <- 400}
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
				remove.otus <- rowSums(curr.ot) < size.thresh; keep.otus <- ! remove.otus;
				#o.t <- otu.table[1:1000, ] + pseudocount;
				o.t <- curr.ot[! remove.otus, ] + pseudocount;
				otus <- rownames(curr.ot);
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
				curr.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
				rownames(curr.sparcc) <- colnames(curr.sparcc) <- otus;
				cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
				########################
				
				########################
				#Tidy up
				rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
				########################
			}
			#Calculate derived OTU similarity matrix S
			cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
			tmp.S <- cor.par(curr.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
			cat(paste("Done =>", Sys.time()), sep="\n");
			curr.S.sparcc <- 0.5 * (tmp.S + 1);
			#Tidy up
			rm(tmp.S, curr.sparcc);
			############################
			#Phylogenetic distances
			#Calculate cophenetic phylogenetic distances
			curr.cophenetic_tree <- cophenetic(curr.tree);
			#Reorder to match order in OTU table
			curr.cophenetic_tree <- curr.cophenetic_tree[rownames(curr.ot), rownames(curr.ot)];
			#Calculate derived OTU similarity matrix S
			cat(paste("Calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
			tmp.S <- cor.par(curr.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
			cat(paste("Done =>", Sys.time()), sep="\n");
			curr.S.phylo <- 0.5 * (tmp.S + 1);
			#Tidy up
			rm(tmp.S, curr.cophenetic_tree);
			############################
			
			
			############################
			#Calculate community similarities
			#=> although most of the classical indices should be impartial to downsampling, they are still re-calculated, just to be sure
			############################
			############################
			#Calculate community similarities
			ds.cs <- list(); t <- 1;
			############################
			#Classical indices
			############################
			#Jaccard index, classical, unweighted
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical"; ds.cs[[t]]$call <- "jaccard"; t <- t + 1;
			#Jaccard index, classical, weighted
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, weighted, no fractions"; ds.cs[[t]]$call <- "jaccard.abd"; t <- t + 1;
			#Jaccard index, classical, weighted, fraction-based
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, weighted"; ds.cs[[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
			#Jaccard index, classical, weighted, fraction-based, Chao's version
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, Chao fraction-based"; ds.cs[[t]]$call <- "jaccard.abd.frac_chao"; t <- t + 1;
			#Jaccard index, classical, weighted, alternative formula
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, alternative formula"; ds.cs[[t]]$call <- "jaccard.abd.alt"; t <- t + 1;
			#Jaccard index, classical, Chao' abundance correction
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Jaccard, classical, Chao"; ds.cs[[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
			#Bray-Curtis, classical
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Bray-Curtis"; ds.cs[[t]]$call <- "bray_curtis"; t <- t + 1;
			#Morisita-Horn, classical
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "Morisita-Horn"; ds.cs[[t]]$call <- "morisita_horn"; t <- t + 1;
			#Iterate through indices and calculate pairwise sample similarities
			for (t in seq(1, length(ds.cs))) {
				print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
				#Calculate all pairwise distances
				my.cs <- community.similarity.par(curr.ot, distance= ds.cs[[t]]$call, use.cores=PARAM$use.cores);
				#Manually correct for rounding errors
				my.cs[my.cs < 0] <- 0;
				#Store current similarities
				ds.cs[[t]]$cs <- my.cs;
				#Report
				print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
			}
			############################
			#UniFrac
			############################
			t <- length(ds.cs) + 1;
			#Unweighted UniFrac
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "UniFrac, unweighted"; ds.cs[[t]]$call <- "UniFrac_unweighted";
			print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
			ds.cs[[t]]$cs <- fasterUniFrac(curr.ot, curr.tree, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
			print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
			t <- t + 1;
			#Weighted UniFrac
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "UniFrac, weighted"; ds.cs[[t]]$call <- "UniFrac_weighted";
			print(paste(Sys.time(), "Starting", ds.cs[[t]]$name));
			ds.cs[[t]]$cs <- fasterUniFrac(curr.ot, curr.tree, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
			print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
			t <- t + 1;
			############################
			#Corrected indices
			############################
			t <- start.t <- length(ds.cs) + 1;
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "TINA, unweighted"; ds.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "TINA, weighted"; ds.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "PINA, unweighted"; ds.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
			ds.cs[[t]] <- list(); ds.cs[[t]]$name <- "PINA, weighted"; ds.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
			############################
			#Iterate through (corrected) indices and calculate pairwise community similarities
			for (t in seq(start.t, length(ds.cs))) {
				print(paste(Sys.time(), bs, "starting", ds.cs[[t]]$name));
				#Calculate all pairwise similarities
				if (ds.cs[[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
					curr.cs <- community.similarity.corr.par(curr.ot, S=curr.S.sparcc, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
				} else {
					curr.cs <- community.similarity.corr.par(curr.ot, S=curr.S.phylo, distance=ds.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
				}
				#Manually correct for rounding errors
				curr.cs[curr.cs < 0] <- 0;
				#Store current similarities
				ds.cs[[t]]$cs <- curr.cs;
				#Report
				print(paste(Sys.time(), "Done with", ds.cs[[t]]$name));
				rm(curr.cs)
			}
			############################
			
			############################
			#Calculate PERMANOVA statistics
			############################
			registerDoMC(cores=length(ds.cs));
			curr.adonis <- foreach(t = seq(1, length(ds.cs))) %dopar% {curr.dist <- as.dist(ds.cs[[t]]$cs); adonis(curr.dist ~ Body_Subsite, data=curr.sample.data)}
			for (t in seq(1, length(all.methods))) {
				curr.idx <- which(bs_ds.res$Step == step & bs_ds.res$Iteration == i & bs_ds.res$Method == all.methods[t] & bs_ds.res$Set == bs);
				#Store current PERMANOVA stats
				bs_ds.res[curr.idx, "F"] <- curr.adonis[[t]]$aov.tab$F.Model[1];
				bs_ds.res[curr.idx, "R2"] <- curr.adonis[[t]]$aov.tab$R2[1];
				bs_ds.res[curr.idx, "P_Value"] <- curr.adonis[[t]]$aov.tab[1,"Pr(>F)"];
				#Store current (averaged) log2FC of pairwise community distances
				curr.ref_cs <- as.dist(get.cs[[t]]$cs[curr.samples, curr.samples]);
				curr.log2fc <- log2(as.dist(ds.cs[[t]]$cs) / curr.ref_cs);
				bs_ds.res[curr.idx, "Mean_Log2FC"] <- mean(curr.log2fc, na.rm=T);
				bs_ds.res[curr.idx, "SD_Log2FC"] <- sd(curr.log2fc, na.rm=T);
			}
			############################
		}
		#Save data
		save(bs_ds.res, file=paste(PARAM$folder.data, "hmp_samples.downsampling_samples.RData", sep=""));
	}
}
############################
################################################################################
################################################################################
#Analyze downsampling runs
bs_ds.res$Samples_per_Class <- as.factor(bs_ds.res$Step);
full.bs_ds.res <- bs_ds.res;
bs_ds.res <- subset(bs_ds.res, Method %in% c("Jaccard, classical", "Jaccard, classical, weighted", "Jaccard, classical, Chao", "Bray-Curtis", "Morisita-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted"));
bs_ds.res$Method <- factor(bs_ds.res$Method, c("Jaccard, classical", "Jaccard, classical, weighted", "Jaccard, classical, Chao", "Bray-Curtis", "Morisita-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted"));
for (bs in c("Body_Site", "Skin", "Urogenital_tract", "Oral")) {
	curr.res <- bs_ds.res[which(bs_ds.res$Set == bs & bs_ds.res$Step != 3849), ];
	#Plot log2(fold change)
	curr.plot <- ggplot(curr.res, aes(x=Samples_per_Class, y=Mean_Log2FC, color=Method, fill=Method)) + geom_boxplot(alpha=0.7, outlier.colour=NA) + geom_point(size=2, position=position_jitterdodge(), alpha=0.8) + scale_color_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + theme_bw();
	ggsave(plot=curr.plot,height=10,width=20,dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.downsampling_sample.", bs, ".log2FC.pdf", sep = ""), useDingbats=FALSE);
	#Plot PERMANOVA F value
	curr.plot <- ggplot(curr.res, aes(x=Samples_per_Class, y=F, color=Method, fill=Method)) + geom_boxplot(alpha=0.7, outlier.colour=NA) + geom_point(size=2, position=position_jitterdodge(), alpha=0.8) + scale_color_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + theme_bw();
	ggsave(plot=curr.plot,height=10,width=20,dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.downsampling_sample.", bs, ".F.pdf", sep = ""), useDingbats=FALSE);
	#Plot PERMANOVA R2 value
	curr.plot <- ggplot(curr.res, aes(x=Samples_per_Class, y=R2, color=Method, fill=Method)) + geom_boxplot(alpha=0.7, outlier.colour=NA) + geom_point(size=2, position=position_jitterdodge(), alpha=0.8) + scale_color_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + theme_bw();
	ggsave(plot=curr.plot,height=10,width=20,dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.downsampling_sample.", bs, ".R2.pdf", sep = ""), useDingbats=FALSE);
	#Plot PERMANOVA -log10(p value)
	curr.plot <- ggplot(curr.res, aes(x=Samples_per_Class, y=P_Value, color=Method, fill=Method)) + geom_boxplot(alpha=0.7, outlier.colour=NA) + geom_point(size=2, position=position_jitterdodge(), alpha=0.8) + scale_y_log10() + scale_color_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#08306b", brewer.pal(n = 10, name="Paired")[3:10])) + theme_bw();
	ggsave(plot=curr.plot,height=10,width=20,dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.downsampling_sample.", bs, ".P_Value.pdf", sep = ""), useDingbats=FALSE);
}
################################################################################
################################################################################




q()